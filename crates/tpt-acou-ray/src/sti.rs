//! Speech Transmission Index (STI) from a room impulse response.
//!
//! The method follows the structure of IEC 60268-16: the impulse response is
//! split into seven octave bands (125 Hz to 8 kHz), the modulation transfer
//! is measured at 14 modulation frequencies per band (0.63 to 12.5 Hz), each
//! value becomes a transmission index `TI = (SNR + 15) / 30` with the apparent
//! SNR clipped to ±15 dB, and the band indices are combined with the
//! weighting factors `α_k` and the redundancy factors `β_k`:
//!
//! ```text
//! STI = Σ α_k · MTI_k − Σ β_k · MTI_k · MTI_{k+1}
//! ```
//!
//! Both weight sets satisfy `Σα − Σβ = 1`, so a perfect channel gives STI = 1.
//! The weights are a parameter, and the male set is the default.
//!
//! **Approximations.** The octave bands are rectangular in the frequency
//! domain (zero-phase, through the FFT), not the standard's band-pass
//! filters. The envelope is the squared band-limited impulse response. The
//! speech-spectrum (intensity-level) weighting of the modulation is not
//! applied. The result is an RIR-based STI estimate. It is not a measured STI
//! and is not certified to IEC 60268-16.

use alloc::vec;
use alloc::vec::Vec;

use tpt_acou_core::error::AcousticError;
use tpt_math_linalg_complex::Complex;

use crate::fft::{fft_in_place, next_pow2};

/// Octave band centres used by STI, Hz.
pub const STI_OCTAVE_HZ: [f64; 7] = [125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0];

/// Modulation frequencies used per band, Hz.
pub const STI_MOD_HZ: [f64; 14] = [
    0.63, 0.80, 1.00, 1.25, 1.60, 2.00, 2.50, 3.15, 4.00, 5.00, 6.30, 8.00, 10.00, 12.50,
];

/// Octave-band weighting factors `α_k` and redundancy factors `β_k`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StiWeights {
    /// Band weights, 125 Hz to 8 kHz.
    pub alpha: [f64; 7],
    /// Adjacent-band redundancy factors, 125/250 Hz to 4/8 kHz.
    pub beta: [f64; 6],
}

/// Male-talker weighting factors (IEC 60268-16 male set, as tabulated in
/// the secondary source documented in `todo.md`). `Σα − Σβ = 1`.
pub const STI_WEIGHTS_MALE: StiWeights = StiWeights {
    alpha: [0.085, 0.127, 0.230, 0.233, 0.309, 0.224, 0.173],
    beta: [0.085, 0.078, 0.065, 0.011, 0.047, 0.095],
};

/// Female-talker weighting factors. The female set is not in IEC 60268-16:2020
/// (it was removed in that edition), so these values come from the older
/// edition via the same secondary source as the male set. Unverified.
/// `Σα − Σβ = 1`.
pub const STI_WEIGHTS_FEMALE: StiWeights = StiWeights {
    alpha: [0.000, 0.117, 0.223, 0.216, 0.328, 0.250, 0.194],
    beta: [0.000, 0.099, 0.066, 0.062, 0.025, 0.076],
};

/// Octave-band impulse responses of `rir`, one per STI band.
fn octave_bands(rir: &[f64], fs: f64) -> Vec<Vec<f64>> {
    // Pad to twice the length so the zero-phase band-pass tail (which extends
    // before time zero) wraps past the end and is discarded by truncation.
    let n = next_pow2(2 * rir.len());
    let mut spectrum = vec![Complex::new(0.0, 0.0); n];
    for (slot, &h) in spectrum.iter_mut().zip(rir) {
        *slot = Complex::new(h, 0.0);
    }
    fft_in_place(&mut spectrum, false);
    STI_OCTAVE_HZ
        .iter()
        .map(|&fc| {
            let lo = fc / core::f64::consts::SQRT_2;
            let hi = fc * core::f64::consts::SQRT_2;
            let mut band = vec![Complex::new(0.0, 0.0); n];
            for k in 0..n {
                // Bin k has frequency k·fs/n, mirrored above n/2.
                let bin = if k <= n / 2 { k } else { n - k };
                let f = bin as f64 * fs / n as f64;
                if f >= lo && f <= hi {
                    band[k] = spectrum[k];
                }
            }
            fft_in_place(&mut band, true);
            band.iter().take(rir.len()).map(|c| c.re).collect()
        })
        .collect()
}

/// Modulation transfer `m(F)` of an energy envelope at modulation frequency
/// `f_mod`: `|Σ e(t) e^{−2πi F t}| / Σ e(t)`.
fn modulation_transfer(envelope: &[f64], fs: f64, f_mod: f64) -> f64 {
    let total: f64 = envelope.iter().sum();
    if total <= 0.0 {
        return 0.0;
    }
    let (mut re, mut im) = (0.0, 0.0);
    for (i, &e) in envelope.iter().enumerate() {
        let phase = -core::f64::consts::TAU * f_mod * i as f64 / fs;
        re += e * libm::cos(phase);
        im += e * libm::sin(phase);
    }
    libm::sqrt(re * re + im * im) / total
}

/// STI of a room impulse response, sampled at `fs`, with the given weights.
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] if `fs` is not positive and
/// finite, if the response is empty or has no energy, or if it contains a
/// non-finite sample.
pub fn sti_from_rir(rir: &[f64], fs: f64, weights: &StiWeights) -> Result<f64, AcousticError> {
    // The top band (8 kHz, edges to 11.3 kHz) must lie below Nyquist.
    if fs.is_finite() && fs < 2.0 * 8_000.0 * core::f64::consts::SQRT_2 {
        return Err(AcousticError::InvalidArgument(
            "STI needs fs above 22.6 kHz so the 8 kHz band lies below Nyquist",
        ));
    }
    if !(fs.is_finite() && fs > 0.0) {
        return Err(AcousticError::InvalidArgument("fs must be finite and > 0"));
    }
    if rir.is_empty() || rir.iter().any(|h| !h.is_finite()) {
        return Err(AcousticError::InvalidArgument(
            "impulse response must be non-empty and finite",
        ));
    }
    if rir.iter().map(|h| h * h).sum::<f64>() <= 0.0 {
        return Err(AcousticError::InvalidArgument(
            "impulse response has no energy",
        ));
    }

    let bands = octave_bands(rir, fs);
    let mut mti = [0.0_f64; 7];
    for (k, band) in bands.iter().enumerate() {
        let envelope: Vec<f64> = band.iter().map(|h| h * h).collect();
        let mut sum = 0.0;
        for &f_mod in &STI_MOD_HZ {
            let m = modulation_transfer(&envelope, fs, f_mod).clamp(0.0, 0.99999);
            // Apparent SNR in dB, clipped to ±15, then the transmission index.
            let snr = (10.0 * libm::log10(m / (1.0 - m))).clamp(-15.0, 15.0);
            sum += (snr + 15.0) / 30.0;
        }
        mti[k] = sum / STI_MOD_HZ.len() as f64;
    }

    let mut sti = 0.0;
    for k in 0..7 {
        sti += weights.alpha[k] * mti[k];
    }
    for k in 0..6 {
        sti -= weights.beta[k] * mti[k] * mti[k + 1];
    }
    Ok(sti.clamp(0.0, 1.0))
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    #[test]
    fn weight_sets_sum_to_one() {
        for w in [STI_WEIGHTS_MALE, STI_WEIGHTS_FEMALE] {
            let a: f64 = w.alpha.iter().sum();
            let b: f64 = w.beta.iter().sum();
            assert!((a - b - 1.0).abs() < 1e-3, "sum alpha {a} beta {b}");
        }
    }

    #[test]
    fn a_perfect_impulse_gives_sti_one() {
        let mut rir = vec![0.0; 4096];
        // Impulse after the band-pass pre-ringing (~180 samples at 16 kHz).
        rir[1600] = 1.0;
        let sti = sti_from_rir(&rir, 48_000.0, &STI_WEIGHTS_MALE).unwrap();
        // Band-limited impulses have finite spread, so the modulation transfer
        // falls slightly below one at 12.5 Hz; the index stays very close to 1.
        assert!(sti > 0.95 && sti <= 1.0, "{sti}");
    }

    #[test]
    fn longer_reverberation_lowers_sti() {
        // Reverberant tails are noise under an exponential envelope. A
        // deterministic exponential with a sharp onset would instead measure
        // the onset transient, which band-pass filtering turns into a short
        // burst that dominates the envelope.
        let fs = 48_000.0;
        let n = 288_000;
        let decay = |t60: f64| -> Vec<f64> {
            let mut rng = crate::rng::Xoshiro256::new(11);
            (0..n)
                .map(|i| {
                    let t = i as f64 / fs;
                    rng.next_normal() * libm::pow(10.0, -3.0 * t / t60)
                })
                .collect()
        };
        let dry = sti_from_rir(&decay(0.4), fs, &STI_WEIGHTS_MALE).unwrap();
        let wet = sti_from_rir(&decay(4.0), fs, &STI_WEIGHTS_MALE).unwrap();
        assert!(dry > wet, "dry {dry} wet {wet}");
        assert!((0.0..=1.0).contains(&dry) && (0.0..=1.0).contains(&wet));
    }

    #[test]
    fn noise_tail_matches_the_closed_form_index() {
        // For an exponential energy envelope with decay rate k = 6 ln10 / T,
        // m(F) = k / √(k² + Ω²), Ω = 2πF. The same m applies in every band,
        // so the STI has a closed form to compare against.
        let fs = 48_000.0;
        let t60 = 1.2;
        let n = 288_000;
        let mut rng = crate::rng::Xoshiro256::new(5);
        let rir: Vec<f64> = (0..n)
            .map(|i| rng.next_normal() * libm::pow(10.0, -3.0 * (i as f64 / fs) / t60))
            .collect();
        let numeric = sti_from_rir(&rir, fs, &STI_WEIGHTS_MALE).unwrap();

        let k = 6.0 * core::f64::consts::LN_10 / t60;
        let mut ti_sum = 0.0;
        for &f in &STI_MOD_HZ {
            let omega = core::f64::consts::TAU * f;
            let m = (k / libm::sqrt(k * k + omega * omega)).clamp(0.0, 0.99999);
            let snr = (10.0 * libm::log10(m / (1.0 - m))).clamp(-15.0, 15.0);
            ti_sum += (snr + 15.0) / 30.0;
        }
        let mti = ti_sum / STI_MOD_HZ.len() as f64;
        let analytic = mti
            * (STI_WEIGHTS_MALE.alpha.iter().sum::<f64>()
                - STI_WEIGHTS_MALE.beta.iter().sum::<f64>() * mti);
        assert!(
            (numeric - analytic).abs() < 0.03,
            "numeric {numeric:.4} vs closed form {analytic:.4}"
        );
    }

    #[test]
    fn rejects_bad_responses() {
        assert!(sti_from_rir(&[], 48_000.0, &STI_WEIGHTS_MALE).is_err());
        assert!(sti_from_rir(&[0.0, 0.0], 48_000.0, &STI_WEIGHTS_MALE).is_err());
        assert!(sti_from_rir(&[1.0], 0.0, &STI_WEIGHTS_MALE).is_err());
        assert!(sti_from_rir(&[1.0], 16_000.0, &STI_WEIGHTS_MALE).is_err());
        assert!(sti_from_rir(&[f64::NAN], 48_000.0, &STI_WEIGHTS_MALE).is_err());
    }
}
