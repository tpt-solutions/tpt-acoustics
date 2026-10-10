//! ISO 3382-1 room acoustic metrics computed from a room impulse response.
//!
//! The input is a sampled impulse response `h[n]` at sample rate `fs`. The
//! energy envelope is `h²`, so:
//!
//! * Clarity `C_t = 10 log10( ∫₀ᵗ h² / ∫ₜ^∞ h² )` (C50 for speech,
//!   C80 for music), in dB.
//! * Definition `D50 = ∫₀^50ms h² / ∫₀^∞ h²`, a ratio in `[0, 1]`.
//! * Reverberation times EDT, T20 and T30 from the Schroeder backward
//!   integral, each extrapolated to a 60 dB decay (see
//!   [`crate::ray::decay_time_from_curve`]).
//!
//! The time windows use the sample grid, so results are exact for the
//! sampled energy and the usual ±1 sample resolution.

use alloc::vec;
use alloc::vec::Vec;

use crate::ray::decay_time_from_curve;

/// Reverberation times in seconds from one impulse response.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReverbTimes {
    /// Early decay time (0 → -10 dB), extrapolated to 60 dB.
    pub edt: f64,
    /// T20 (-5 → -25 dB), extrapolated to 60 dB.
    pub t20: f64,
    /// T30 (-5 → -35 dB), extrapolated to 60 dB.
    pub t30: f64,
}

/// Sum of `h²` up to but excluding sample `split`, and total energy.
fn early_and_total(rir: &[f64], split: usize) -> (f64, f64) {
    let mut early = 0.0;
    let mut total = 0.0;
    for (n, &h) in rir.iter().enumerate() {
        let e = h * h;
        total += e;
        if n < split {
            early += e;
        }
    }
    (early, total)
}

/// Sample index corresponding to `t_ms` milliseconds at `fs` Hz.
fn split_index(t_ms: f64, fs: f64) -> usize {
    libm::round(t_ms * 1e-3 * fs) as usize
}

/// Clarity `C_t` in dB for an early window of `t_ms` milliseconds.
///
/// Returns `None` if the impulse response is empty, if `fs` or `t_ms` is not
/// positive, or if there is no energy after the window (the ratio would be
/// undefined).
#[must_use]
pub fn clarity_db(rir: &[f64], fs: f64, t_ms: f64) -> Option<f64> {
    if rir.is_empty() || fs <= 0.0 || t_ms <= 0.0 {
        return None;
    }
    let (early, total) = early_and_total(rir, split_index(t_ms, fs));
    let late = total - early;
    if early <= 0.0 || late <= 0.0 {
        return None;
    }
    Some(10.0 * libm::log10(early / late))
}

/// Definition `D` for an early window of `t_ms` milliseconds, in `[0, 1]`.
///
/// Returns `None` for an empty or all-zero impulse response, or invalid
/// `fs`/`t_ms`.
#[must_use]
pub fn definition(rir: &[f64], fs: f64, t_ms: f64) -> Option<f64> {
    if rir.is_empty() || fs <= 0.0 || t_ms <= 0.0 {
        return None;
    }
    let (early, total) = early_and_total(rir, split_index(t_ms, fs));
    if total <= 0.0 {
        return None;
    }
    Some(early / total)
}

/// Schroeder backward-integrated energy curve of an impulse response:
/// `curve[i] = Σ_{j ≥ i} h[j]²`.
#[must_use]
pub fn schroeder_from_rir(rir: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; rir.len()];
    let mut acc = 0.0;
    for (i, &h) in rir.iter().enumerate().rev() {
        acc += h * h;
        out[i] = acc;
    }
    out
}

/// EDT, T20 and T30 of an impulse response. A field is `NaN` if the response
/// does not decay far enough to fit its window (or `fs` is not positive).
#[must_use]
pub fn reverberation_times(rir: &[f64], fs: f64) -> ReverbTimes {
    let curve = schroeder_from_rir(rir);
    let bin = if fs > 0.0 { 1.0 / fs } else { 0.0 };
    let fit = |upper, lower| {
        if fs > 0.0 {
            decay_time_from_curve(&curve, bin, upper, lower).unwrap_or(f64::NAN)
        } else {
            f64::NAN
        }
    };
    ReverbTimes {
        edt: fit(0.0, -10.0),
        t20: fit(-5.0, -25.0),
        t30: fit(-5.0, -35.0),
    }
}

/// Schroeder frequency in Hz, above which a room's modal overlap is high
/// enough for statistical (ray and diffuse-field) methods to apply:
/// `f_S = 2000 √(T / V)` with `T` the reverberation time in seconds and `V`
/// the room volume in cubic metres. Hybrid ISM/ray/wave crossovers use it.
///
/// Returns `None` for non-positive or non-finite `t60_s` or `volume_m3`.
#[must_use]
pub fn schroeder_frequency_hz(t60_s: f64, volume_m3: f64) -> Option<f64> {
    if !(t60_s.is_finite() && t60_s > 0.0 && volume_m3.is_finite() && volume_m3 > 0.0) {
        return None;
    }
    Some(2000.0 * libm::sqrt(t60_s / volume_m3))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn schroeder_frequency_worked_example() {
        // A 200 m³ room with T = 0.5 s: f_S = 2000 √(0.5/200) = 100 Hz.
        let f = schroeder_frequency_hz(0.5, 200.0).unwrap();
        assert!((f - 100.0).abs() < 1e-9, "{f}");
        // Scales with √T and 1/√V.
        let f4 = schroeder_frequency_hz(2.0, 200.0).unwrap();
        assert!((f4 / f - 2.0).abs() < 1e-12);
        assert!(schroeder_frequency_hz(0.0, 200.0).is_none());
        assert!(schroeder_frequency_hz(0.5, -1.0).is_none());
    }

    /// Deterministic exponential impulse response with a 60 dB energy decay
    /// in `t60` seconds: amplitude `10^(-3 t / t60)`, energy `10^(-6 t / t60)`.
    fn exponential_rir(t60: f64, fs: f64, seconds: f64) -> Vec<f64> {
        let n = libm::round(seconds * fs) as usize;
        (0..n)
            .map(|i| libm::pow(10.0, -3.0 * (i as f64 / fs) / t60))
            .collect()
    }

    #[test]
    fn decay_times_recover_exact_exponential_t60() {
        let fs = 8000.0;
        let rir = exponential_rir(1.0, fs, 2.5);
        let t = reverberation_times(&rir, fs);
        assert!((t.t30 - 1.0).abs() < 0.02, "T30 {}", t.t30);
        assert!((t.t20 - 1.0).abs() < 0.02, "T20 {}", t.t20);
        assert!((t.edt - 1.0).abs() < 0.02, "EDT {}", t.edt);
    }

    #[test]
    fn clarity_matches_closed_form_for_exponential_decay() {
        // For energy e^(-k t) with k = 6 ln10 / T: C_t = 10 log10(e^(k t) - 1).
        let fs = 16000.0;
        let t60 = 1.0;
        let rir = exponential_rir(t60, fs, 4.0);
        let k = 6.0 * core::f64::consts::LN_10 / t60;
        let c50 = clarity_db(&rir, fs, 50.0).unwrap();
        let expected = 10.0 * libm::log10(libm::exp(k * 0.050) - 1.0);
        assert!((c50 - expected).abs() < 0.05, "C50 {c50} vs {expected}");
    }

    #[test]
    fn definition_matches_closed_form_for_exponential_decay() {
        // D50 = 1 - e^(-k · 0.05) for the same energy decay.
        let fs = 16000.0;
        let rir = exponential_rir(1.0, fs, 4.0);
        let k = 6.0 * core::f64::consts::LN_10;
        let d50 = definition(&rir, fs, 50.0).unwrap();
        let expected = 1.0 - libm::exp(-k * 0.050);
        assert!((d50 - expected).abs() < 1e-3, "D50 {d50} vs {expected}");
        assert!((0.0..=1.0).contains(&d50));
    }

    #[test]
    fn c80_exceeds_c50_for_the_same_response() {
        let fs = 16000.0;
        let rir = exponential_rir(1.0, fs, 3.0);
        let c50 = clarity_db(&rir, fs, 50.0).unwrap();
        let c80 = clarity_db(&rir, fs, 80.0).unwrap();
        assert!(c80 > c50);
    }

    #[test]
    fn degenerate_inputs_return_none() {
        assert!(clarity_db(&[], 1000.0, 50.0).is_none());
        assert!(clarity_db(&[1.0, 0.5], 0.0, 50.0).is_none());
        assert!(definition(&[0.0, 0.0], 1000.0, 50.0).is_none());
        // Entire energy inside the early window: no late energy.
        assert!(clarity_db(&[1.0, 0.0], 1000.0, 50.0).is_none());
    }
}
