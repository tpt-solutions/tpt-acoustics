//! Broadband hybrid impulse response: modal below the Schroeder frequency,
//! image sources above it.
//!
//! The two transfer functions share the free-field Green's-function
//! normalisation, so their spectra join directly. The combined spectrum has
//! Hermitian symmetry, and the inverse FFT gives a real impulse response.
//!
//! * Below `f_S = 2000 √(T/V)` the response is modal, and the rigid-box modal
//!   sum ([`crate::modal::RigidBoxModes`]) is used, with the Sabine damping of
//!   the room at 1 kHz.
//! * Above `f_S` the image-source transfer function is used, with reflection
//!   coefficients from the wall absorption.
//!
//! **Limits.** The high band uses image sources only. The ray-traced late
//! reverberation is not merged here, because the ray tracer is energy-only and
//! has no phase. The crossover is a hard cut at `f_S`. The result is a
//! frequency-domain construction of the Green's function, with no separate
//! level calibration between the two bands.

use alloc::vec;
use alloc::vec::Vec;

use tpt_acou_core::error::AcousticError;
use tpt_acou_core::medium::Medium;
use tpt_math_linalg_complex::Complex;

use crate::fft::{fft_in_place, next_pow2};
use crate::ism::ShoeboxRoom;
use crate::modal::RigidBoxModes;
use crate::vec3::Vec3;

/// Broadband impulse response of `room` between `src` and `rcv`.
///
/// `fs` is the sample rate, and the output has `n_samples` samples, rounded up
/// to a power of two internally and truncated back. The image-source order is
/// `max_order`.
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] for non-positive `fs` or
/// `n_samples`, and propagates geometry errors for positions outside the room.
pub fn hybrid_broadband_ir(
    room: &ShoeboxRoom,
    src: Vec3,
    rcv: Vec3,
    fs: f64,
    n_samples: usize,
    max_order: u32,
    medium: &Medium,
) -> Result<Vec<f64>, AcousticError> {
    if !(fs.is_finite() && fs > 0.0) || n_samples == 0 {
        return Err(AcousticError::InvalidArgument(
            "fs must be positive and n_samples at least 1",
        ));
    }
    let dims = room.dims();
    let volume = dims[0] * dims[1] * dims[2];
    let t60 = room.sabine_rt60(1000.0);
    let f_schroeder = crate::iso3382::schroeder_frequency_hz(t60, volume)
        .ok_or(AcousticError::InvalidArgument(
            "room has no finite reverberation time",
        ))?
        .min(fs / 2.0);

    let n = next_pow2(n_samples);
    let half = n / 2;
    let modes = RigidBoxModes::new(
        dims,
        [src.x, src.y, src.z],
        [rcv.x, rcv.y, rcv.z],
        medium.speed_of_sound,
        2.0 * f_schroeder,
        t60,
    )?;

    // Spectrum on the FFT grid, bins 0..=N/2.
    let mut spectrum = vec![Complex::new(0.0, 0.0); n];
    for k in 1..=half {
        let f = k as f64 * fs / n as f64;
        let value = if f < f_schroeder {
            modes.green_at_hz(f)
        } else {
            room.transfer_function(src, rcv, &[f], max_order, medium)?[0]
        };
        spectrum[k] = value;
        if k != half {
            spectrum[n - k] = value.conj();
        }
    }
    // Inverse transform: h[n] = fs · IFFT(H)[n].
    fft_in_place(&mut spectrum, true);
    let out: Vec<f64> = spectrum.iter().take(n_samples).map(|c| c.re * fs).collect();
    Ok(out)
}

/// Result of [`three_way_ir`].
#[derive(Debug, Clone, PartialEq)]
pub struct ThreeWayIr {
    /// The impulse response, `n_samples` long.
    pub samples: Vec<f64>,
    /// Gain applied to the ray tail: the broadband energy in the calibration
    /// window divided by the ray energy in the same window.
    pub tail_gain: f64,
}

/// Three-way impulse response: modal low frequencies and image-source high
/// frequencies up to the crossover, then a ray-traced energy tail.
///
/// Before `crossover_s` the samples are the signed broadband response
/// ([`hybrid_broadband_ir`]), so modal low-frequency structure and the
/// image-source early reflections both appear. From the crossover on, the
/// samples come from the ray captures scaled by `tail_gain`, placed with the
/// energy-only convention of [`crate::rir::rir_from_echogram`].
///
/// The gain matches the broadband energy in `±window_s` of the crossover to
/// the ray energy in the same window, the same energy calibration as
/// [`crate::rir::hybrid_rir`]. The tail is therefore broadband: it has no
/// modal low-frequency structure. Long low-frequency modal decay that runs
/// past `crossover_s` is lost, so choose a crossover late enough for the
/// modal part to have decayed, or accept the cut.
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] for invalid sampling, times or
/// energies, mismatched ray slices, a crossover outside the response, or no ray
/// energy (or no broadband energy) in the calibration window.
#[allow(clippy::too_many_arguments)]
pub fn three_way_ir(
    room: &ShoeboxRoom,
    src: Vec3,
    rcv: Vec3,
    fs: f64,
    n_samples: usize,
    max_order: u32,
    medium: &Medium,
    ray_times: &[f64],
    ray_energies: &[f64],
    crossover_s: f64,
    window_s: f64,
    seed: u64,
) -> Result<ThreeWayIr, AcousticError> {
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
    if ray_times
        .iter()
        .zip(ray_energies)
        .any(|(&t, &e)| !(t.is_finite() && t >= 0.0 && e.is_finite() && e >= 0.0))
    {
        return Err(AcousticError::InvalidArgument(
            "ray times and energies must be finite and >= 0",
        ));
    }
    let mut samples = hybrid_broadband_ir(room, src, rcv, fs, n_samples, max_order, medium)?;
    let split = libm::round(crossover_s * fs) as usize;
    if split >= n_samples {
        return Err(AcousticError::InvalidArgument(
            "crossover lies beyond the end of the response",
        ));
    }

    // Calibration: broadband energy and ray energy within ±window_s.
    let in_window = |t: f64| (t - crossover_s).abs() <= window_s;
    let broadband_energy: f64 = samples
        .iter()
        .enumerate()
        .filter(|&(i, _)| in_window(i as f64 / fs))
        .map(|(_, &h)| h * h)
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
    if broadband_energy <= 0.0 {
        return Err(AcousticError::InvalidArgument(
            "no broadband energy in the calibration window",
        ));
    }
    let tail_gain = broadband_energy / ray_energy;

    let late: Vec<(f64, f64)> = ray_times
        .iter()
        .zip(ray_energies)
        .filter(|&(&t, _)| t >= crossover_s)
        .map(|(&t, &e)| (t, e * tail_gain))
        .collect();
    samples.truncate(split);
    samples.resize(n_samples, 0.0);
    if !late.is_empty() {
        let tail = crate::rir::rir_from_echogram(&late, fs, seed)?;
        for (slot, &h) in samples[split..]
            .iter_mut()
            .zip(&tail[split.min(tail.len())..])
        {
            *slot = h;
        }
    }
    Ok(ThreeWayIr { samples, tail_gain })
}

/// The spectrum on the same grid as [`hybrid_broadband_ir`] (bins 0..=N/2 for
/// the padded length), for checks and plots. It uses the same crossover.
///
/// # Errors
/// As [`hybrid_broadband_ir`].
pub fn hybrid_broadband_spectrum(
    room: &ShoeboxRoom,
    src: Vec3,
    rcv: Vec3,
    fs: f64,
    n_samples: usize,
    max_order: u32,
    medium: &Medium,
) -> Result<Vec<Complex>, AcousticError> {
    let ir = hybrid_broadband_ir(room, src, rcv, fs, n_samples, max_order, medium)?;
    let n = next_pow2(n_samples);
    let mut padded: Vec<Complex> = vec![Complex::new(0.0, 0.0); n];
    for (slot, &h) in padded.iter_mut().zip(&ir) {
        *slot = Complex::new(h, 0.0);
    }
    fft_in_place(&mut padded, false);
    Ok(padded.iter().map(|c| *c * (1.0 / fs)).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::BoundaryMaterial;

    fn room() -> ShoeboxRoom {
        ShoeboxRoom::uniform([5.0, 4.0, 3.0], BoundaryMaterial::concrete_painted()).unwrap()
    }

    #[test]
    fn output_is_real_and_has_the_requested_length() {
        let air = Medium::air_20c();
        let ir = hybrid_broadband_ir(
            &room(),
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(4.0, 3.0, 2.0),
            2_000.0,
            1_000,
            2,
            &air,
        )
        .unwrap();
        assert_eq!(ir.len(), 1_000);
        assert!(ir.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn spectrum_matches_the_two_transfer_functions_in_their_bands() {
        let air = Medium::air_20c();
        let r = room();
        let src = Vec3::new(1.0, 1.0, 1.0);
        let rcv = Vec3::new(4.0, 3.0, 2.0);
        let fs = 2_000.0;
        let n_samples = 1_024;
        let spec = hybrid_broadband_spectrum(&r, src, rcv, fs, n_samples, 3, &air).unwrap();
        let n = next_pow2(n_samples);
        let t60 = r.sabine_rt60(1000.0);
        let f_s = crate::iso3382::schroeder_frequency_hz(t60, 60.0)
            .unwrap()
            .min(fs / 2.0);
        let modes = RigidBoxModes::new(
            r.dims(),
            [src.x, src.y, src.z],
            [rcv.x, rcv.y, rcv.z],
            air.speed_of_sound,
            2.0 * f_s,
            t60,
        )
        .unwrap();
        // Bins well inside each band (away from the crossover bin).
        for k in [1_usize, 3, 7] {
            let f = k as f64 * fs / n as f64;
            assert!(f < f_s, "test bin {k} must lie below f_S = {f_s}");
            let want = modes.green_at_hz(f);
            assert!(
                (spec[k] - want).norm() < 1e-9 * want.norm().max(1e-30),
                "bin {k}"
            );
        }
        for k in [n / 3, n / 2 - 10] {
            let f = k as f64 * fs / n as f64;
            assert!(f > f_s, "test bin {k} must lie above f_S");
            let want = r.transfer_function(src, rcv, &[f], 3, &air).unwrap()[0];
            assert!(
                (spec[k] - want).norm() < 1e-9 * want.norm().max(1e-30),
                "bin {k}"
            );
        }
    }

    #[test]
    fn rejects_invalid_sampling() {
        let air = Medium::air_20c();
        assert!(hybrid_broadband_ir(
            &room(),
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(4.0, 3.0, 2.0),
            0.0,
            10,
            1,
            &air
        )
        .is_err());
        assert!(hybrid_broadband_ir(
            &room(),
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(4.0, 3.0, 2.0),
            1000.0,
            0,
            1,
            &air
        )
        .is_err());
    }

    /// Ray captures 1 ms apart with unit energy, so no two impulses share a
    /// sample and the tail's energy is exactly the gain times the ray energy.
    fn synthetic_rays(start_s: f64, count: usize) -> (Vec<f64>, Vec<f64>) {
        let times = (0..count).map(|i| start_s + i as f64 * 1e-3).collect();
        (times, vec![1.0; count])
    }

    #[test]
    fn three_way_matches_broadband_before_the_crossover() {
        let air = Medium::air_20c();
        let r = room();
        let src = Vec3::new(1.0, 1.0, 1.0);
        let rcv = Vec3::new(4.0, 3.0, 2.0);
        let (fs, n, max_order) = (2_000.0, 1_024, 2);
        let (times, energies) = synthetic_rays(0.2, 200);
        let crossover = 0.2;
        let out = three_way_ir(
            &r, src, rcv, fs, n, max_order, &air, &times, &energies, crossover, 0.05, 7,
        )
        .unwrap();
        let broadband = hybrid_broadband_ir(&r, src, rcv, fs, n, max_order, &air).unwrap();
        let split = (crossover * fs).round() as usize;
        assert_eq!(out.samples.len(), n);
        assert!(out.samples.iter().all(|v| v.is_finite()));
        assert_eq!(&out.samples[..split], &broadband[..split]);
        assert!(out.tail_gain > 0.0 && out.tail_gain.is_finite());
    }

    #[test]
    fn three_way_tail_energy_equals_the_calibrated_ray_energy() {
        let air = Medium::air_20c();
        let r = room();
        let src = Vec3::new(1.0, 1.0, 1.0);
        let rcv = Vec3::new(4.0, 3.0, 2.0);
        let (fs, n) = (2_000.0, 1_024);
        let (times, energies) = synthetic_rays(0.2, 200);
        let out = three_way_ir(
            &r, src, rcv, fs, n, 2, &air, &times, &energies, 0.2, 0.05, 7,
        )
        .unwrap();
        let split = (0.2 * fs).round() as usize;
        let tail_energy: f64 = out.samples[split..].iter().map(|h| h * h).sum();
        let ray_energy: f64 = energies
            .iter()
            .zip(&times)
            .filter(|&(_, &t)| t >= 0.2)
            .map(|(&e, _)| e)
            .sum();
        assert!(
            (tail_energy - out.tail_gain * ray_energy).abs() < 1e-9 * ray_energy,
            "tail {tail_energy} vs {}",
            out.tail_gain * ray_energy
        );
    }

    #[test]
    fn three_way_rejects_bad_inputs() {
        let air = Medium::air_20c();
        let r = room();
        let src = Vec3::new(1.0, 1.0, 1.0);
        let rcv = Vec3::new(4.0, 3.0, 2.0);
        let (times, energies) = synthetic_rays(0.2, 200);
        // Crossover past the end of a 1024-sample response at 2 kHz (0.512 s).
        assert!(three_way_ir(
            &r, src, rcv, 2_000.0, 1_024, 2, &air, &times, &energies, 0.6, 0.05, 1
        )
        .is_err());
        // No ray energy near the crossover: the last capture is at 0.399 s.
        assert!(three_way_ir(
            &r, src, rcv, 2_000.0, 1_024, 2, &air, &times, &energies, 0.45, 0.01, 1
        )
        .is_err());
        // Mismatched ray slices.
        assert!(three_way_ir(
            &r,
            src,
            rcv,
            2_000.0,
            1_024,
            2,
            &air,
            &times,
            &energies[1..],
            0.2,
            0.05,
            1
        )
        .is_err());
    }
}
