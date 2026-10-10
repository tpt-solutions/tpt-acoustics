//! Room impulse response (RIR) construction and WAV export.
//!
//! The image-source echogram carries energy weights, not signed pressures.
//! [`rir_from_echogram`] turns each weight `w` into an impulse of amplitude
//! `sqrt(w)` with a pseudo-random sign, placed at the nearest sample. Random
//! signs keep the impulses incoherent, so the energy envelope is correct
//! while the waveform is only a statistical stand-in for a real RIR. Use it
//! for level and decay analysis and for auralization previews, not for phase.
//!
//! [`write_wav_f32`] serialises a response as a mono 32-bit IEEE float WAV
//! file (RIFF format tag 3), which every common audio tool reads.

use alloc::vec;
use alloc::vec::Vec;

use tpt_acou_core::error::AcousticError;

use crate::rng::Xoshiro256;

/// Builds a sampled impulse response from an echogram of
/// `(arrival_time_s, energy_weight)` pairs at sample rate `fs`.
///
/// The output runs from time zero to the latest arrival (rounded up to a
/// whole sample). Impulses landing on the same sample add. The signs are
/// drawn from a deterministic generator seeded with `seed`.
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] if `fs` is not positive and
/// finite, if `echogram` is empty, or if any arrival time is negative or
/// non-finite, or any weight is negative or non-finite.
pub fn rir_from_echogram(
    echogram: &[(f64, f64)],
    fs: f64,
    seed: u64,
) -> Result<Vec<f64>, AcousticError> {
    if !(fs.is_finite() && fs > 0.0) {
        return Err(AcousticError::InvalidArgument("fs must be finite and > 0"));
    }
    if echogram.is_empty() {
        return Err(AcousticError::InvalidArgument("echogram is empty"));
    }
    if echogram
        .iter()
        .any(|&(t, w)| !(t.is_finite() && t >= 0.0 && w.is_finite() && w >= 0.0))
    {
        return Err(AcousticError::InvalidArgument(
            "arrival times and weights must be finite and >= 0",
        ));
    }
    let t_max = echogram.iter().map(|&(t, _)| t).fold(0.0_f64, f64::max);
    let len = libm::ceil(t_max * fs) as usize + 1;
    let mut h = vec![0.0; len];
    let mut rng = Xoshiro256::new(seed);
    for &(t, w) in echogram {
        let n = libm::round(t * fs) as usize;
        let sign = if rng.next_f64() < 0.5 { -1.0 } else { 1.0 };
        h[n] += sign * libm::sqrt(w);
    }
    Ok(h)
}

/// Result of [`hybrid_rir`].
#[derive(Debug, Clone, PartialEq)]
pub struct HybridRir {
    /// The impulse response, sampled at the requested rate.
    pub samples: Vec<f64>,
    /// Gain applied to the ray tail. It equals the ISM energy in the
    /// calibration window divided by the ray energy in the same window.
    pub tail_gain: f64,
}

/// Hybrid impulse response: image-source early part, ray-traced late part.
///
/// The ISM echogram gives absolute energy weights, and the ray tracer gives
/// energies normalised to the injected total, so the two cannot be joined
/// without a common scale. This function uses the ISM energy as the
/// reference: it measures both energies within `±window_s` of
/// `crossover_s`, and scales the ray captures by the ratio of the two.
/// Before the crossover the response comes from the echogram. From the
/// crossover on it comes from the scaled ray captures. The energy level is
/// therefore continuous at the join. The calibration is an energy match, not
/// a physical derivation, and the gain is returned so callers can check it.
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] if `fs`, `crossover_s` or
/// `window_s` is not positive and finite, if the ray time and energy slices
/// differ in length, if any time or energy is negative or non-finite, or if
/// no ray capture falls in the calibration window (no scale can be set).
pub fn hybrid_rir(
    echogram: &[(f64, f64)],
    ray_times: &[f64],
    ray_energies: &[f64],
    fs: f64,
    crossover_s: f64,
    window_s: f64,
    seed: u64,
) -> Result<HybridRir, AcousticError> {
    if !(fs.is_finite() && fs > 0.0) {
        return Err(AcousticError::InvalidArgument("fs must be finite and > 0"));
    }
    if !(crossover_s.is_finite() && crossover_s > 0.0) {
        return Err(AcousticError::InvalidArgument(
            "crossover must be finite and > 0",
        ));
    }
    if !(window_s.is_finite() && window_s > 0.0) {
        return Err(AcousticError::InvalidArgument(
            "window must be finite and > 0",
        ));
    }
    if ray_times.len() != ray_energies.len() {
        return Err(AcousticError::InvalidArgument(
            "ray times and energies must have the same length",
        ));
    }
    if echogram
        .iter()
        .any(|&(t, w)| !(t.is_finite() && t >= 0.0 && w.is_finite() && w >= 0.0))
        || ray_times
            .iter()
            .zip(ray_energies)
            .any(|(&t, &e)| !(t.is_finite() && t >= 0.0 && e.is_finite() && e >= 0.0))
    {
        return Err(AcousticError::InvalidArgument(
            "times and energies must be finite and >= 0",
        ));
    }

    let in_window = |t: f64| (t - crossover_s).abs() <= window_s;
    let ism_energy: f64 = echogram
        .iter()
        .filter(|&&(t, _)| in_window(t))
        .map(|&(_, w)| w)
        .sum();
    let ray_energy: f64 = ray_times
        .iter()
        .zip(ray_energies)
        .filter(|&(&t, _)| in_window(t))
        .map(|(_, &e)| e)
        .sum();
    if ray_energy <= 0.0 {
        return Err(AcousticError::InvalidArgument(
            "no ray energy in the calibration window; widen it or trace more rays",
        ));
    }
    let tail_gain = ism_energy / ray_energy;

    // Early part: echogram arrivals before the crossover.
    let early: Vec<(f64, f64)> = echogram
        .iter()
        .copied()
        .filter(|&(t, _)| t < crossover_s)
        .collect();
    // Late part: scaled ray captures from the crossover on.
    let late: Vec<(f64, f64)> = ray_times
        .iter()
        .zip(ray_energies)
        .filter(|&(&t, _)| t >= crossover_s)
        .map(|(&t, &e)| (t, e * tail_gain))
        .collect();

    let mut combined = early;
    combined.extend(late);
    let samples = rir_from_echogram(&combined, fs, seed)?;
    Ok(HybridRir { samples, tail_gain })
}

/// Serialises `samples` (nominally in `[-1, 1]`) as a mono 32-bit float WAV
/// file at `sample_rate_hz`. Values are written as-is; the caller is
/// responsible for scaling.
#[must_use]
pub fn write_wav_f32(samples: &[f64], sample_rate_hz: u32) -> Vec<u8> {
    let data_bytes = samples.len() as u32 * 4;
    let mut out = Vec::with_capacity(44 + data_bytes as usize);
    out.extend_from_slice(b"RIFF");
    out.extend_from_slice(&(36 + data_bytes).to_le_bytes());
    out.extend_from_slice(b"WAVE");
    out.extend_from_slice(b"fmt ");
    out.extend_from_slice(&16_u32.to_le_bytes());
    out.extend_from_slice(&3_u16.to_le_bytes()); // IEEE float
    out.extend_from_slice(&1_u16.to_le_bytes()); // mono
    out.extend_from_slice(&sample_rate_hz.to_le_bytes());
    out.extend_from_slice(&(sample_rate_hz * 4).to_le_bytes()); // byte rate
    out.extend_from_slice(&4_u16.to_le_bytes()); // block align
    out.extend_from_slice(&32_u16.to_le_bytes()); // bits per sample
    out.extend_from_slice(b"data");
    out.extend_from_slice(&data_bytes.to_le_bytes());
    for &s in samples {
        out.extend_from_slice(&(s as f32).to_le_bytes());
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn echogram_energy_is_preserved_for_separated_arrivals() {
        let fs = 1000.0;
        // Arrivals on distinct samples, so no coherent addition.
        let echo = [(0.010, 0.5), (0.020, 0.25), (0.100, 0.125)];
        let h = rir_from_echogram(&echo, fs, 1).unwrap();
        let energy: f64 = h.iter().map(|x| x * x).sum();
        let expected: f64 = echo.iter().map(|&(_, w)| w).sum();
        assert!((energy - expected).abs() < 1e-12, "{energy} vs {expected}");
        // Each impulse sits on its nearest sample.
        assert!((h[10].abs() - libm::sqrt(0.5)).abs() < 1e-12);
        assert!((h[100].abs() - libm::sqrt(0.125)).abs() < 1e-12);
    }

    #[test]
    fn echogram_rir_is_deterministic_for_a_seed() {
        let echo = [(0.01, 0.5), (0.03, 0.2), (0.05, 0.1)];
        let a = rir_from_echogram(&echo, 2000.0, 42).unwrap();
        let b = rir_from_echogram(&echo, 2000.0, 42).unwrap();
        assert_eq!(a, b);
    }

    #[test]
    fn echogram_rejects_bad_input() {
        assert!(rir_from_echogram(&[], 1000.0, 0).is_err());
        assert!(rir_from_echogram(&[(0.1, 1.0)], 0.0, 0).is_err());
        assert!(rir_from_echogram(&[(-0.1, 1.0)], 1000.0, 0).is_err());
        assert!(rir_from_echogram(&[(0.1, -1.0)], 1000.0, 0).is_err());
    }

    #[test]
    fn hybrid_tail_is_calibrated_to_the_echogram_at_the_crossover() {
        let fs = 1000.0;
        // ISM window (50 ms ± 10 ms) holds 0.4; the ray window holds 0.3.
        let echo = [(0.010, 0.5), (0.045, 0.4)];
        let rays_t = [0.052, 0.058, 0.200];
        let rays_e = [0.1, 0.2, 0.05];
        let hybrid = hybrid_rir(&echo, &rays_t, &rays_e, fs, 0.05, 0.01, 1).unwrap();
        assert!(
            (hybrid.tail_gain - 0.4 / 0.3).abs() < 1e-12,
            "{}",
            hybrid.tail_gain
        );
        // Post-crossover energy inside the window equals the ISM reference.
        // Indices use the same rounding as `rir_from_echogram`.
        let i0 = libm::round(0.052 * fs) as usize;
        let i1 = libm::round(0.058 * fs) as usize;
        let window_energy: f64 = hybrid.samples[i0..=i1].iter().map(|h| h * h).sum();
        assert!((window_energy - 0.4).abs() < 1e-9, "{window_energy}");
        // The early echogram arrival is kept as-is, and the late ray arrival
        // is scaled by the same gain.
        assert!((hybrid.samples[10].abs() - libm::sqrt(0.5)).abs() < 1e-12);
        assert!((hybrid.samples[200].abs() - libm::sqrt(0.05 * hybrid.tail_gain)).abs() < 1e-12);
    }

    #[test]
    fn hybrid_rejects_bad_calibration_inputs() {
        let echo = [(0.010, 0.5)];
        // No ray energy in the window.
        assert!(hybrid_rir(&echo, &[0.3], &[0.1], 1000.0, 0.05, 0.01, 1).is_err());
        // Mismatched slices.
        assert!(hybrid_rir(&echo, &[0.05], &[], 1000.0, 0.05, 0.01, 1).is_err());
        // Non-positive window.
        assert!(hybrid_rir(&echo, &[0.05], &[0.1], 1000.0, 0.05, 0.0, 1).is_err());
    }

    #[test]
    fn wav_header_describes_float_mono_data() {
        let samples = [0.0, 0.5, -0.5, 1.0];
        let wav = write_wav_f32(&samples, 48000);
        assert_eq!(&wav[0..4], b"RIFF");
        assert_eq!(&wav[8..12], b"WAVE");
        assert_eq!(u16::from_le_bytes([wav[20], wav[21]]), 3); // IEEE float
        assert_eq!(u16::from_le_bytes([wav[22], wav[23]]), 1); // mono
        assert_eq!(
            u32::from_le_bytes([wav[24], wav[25], wav[26], wav[27]]),
            48000
        );
        assert_eq!(u16::from_le_bytes([wav[34], wav[35]]), 32); // bits
        assert_eq!(&wav[36..40], b"data");
        let data_len = u32::from_le_bytes([wav[40], wav[41], wav[42], wav[43]]) as usize;
        assert_eq!(data_len, samples.len() * 4);
        assert_eq!(wav.len(), 44 + data_len);
        let second = f32::from_le_bytes([wav[48], wav[49], wav[50], wav[51]]);
        assert_eq!(second, 0.5);
    }
}
