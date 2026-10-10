//! Moore–Glasberg loudness model for tonal sounds.
//!
//! The structure follows Glasberg & Moore (1990) and Moore, Glasberg & Baer
//! (1997):
//!
//! * Auditory filters are rounded-exponential, with width set by the
//!   equivalent rectangular bandwidth `ERB(f) = 24.7 (4.37 f/1000 + 1)` Hz.
//!   A tone at `f` produces excitation `P·(1 + p g) e^(−p g)` at a channel
//!   centred on `fc`, where `g = |f − fc|/fc` and `p = 4 fc / ERB(fc)`.
//! * The absolute threshold is `T(f) = 3.64 (f/1000)^−0.8 − 6.5 e^(−0.6 (f/1000 − 3.3)²)`
//!   `+ 10⁻³ (f/1000)⁴` dB SPL. A channel's excitation threshold is its threshold level
//!   converted to a power.
//! * Specific loudness is compressive, `N'(fc) = C [(E + E_thr)^α − E_thr^α]`
//!   with `α = 0.2` and `C = 0.046`, and total loudness is the integral over
//!   channels spaced by 0.1 ERB-number.
//!
//! **Calibration and limits.** The model is scaled so that a 1 kHz tone at
//! 40 dB SPL gives exactly one sone, which is the definition of the phon
//! scale. Other levels and frequencies follow from the model's filters and
//! constants. The outer- and middle-ear filter, the binaural and temporal
//! effects, and partial loudness are not modelled. This is a tonal model: it
//! takes a list of pure tones, not a spectrum. The absolute values have not
//! been checked against the ISO 532-2 reference data, which are licensed.

use tpt_acou_core::error::AcousticError;

/// Specific-loudness compression exponent.
const ALPHA: f64 = 0.2;
/// Specific-loudness scale constant.
const C_SCALE: f64 = 0.046;
/// Channel spacing in ERB-number units.
const ERB_STEP: f64 = 0.1;

/// Equivalent rectangular bandwidth of the auditory filter at `f` Hz.
#[must_use]
pub fn erb_hz(f: f64) -> f64 {
    24.7 * (4.37 * f / 1000.0 + 1.0)
}

/// ERB-number scale: `21.4 log10(4.37 f/1000 + 1)`.
#[must_use]
pub fn erb_number(f: f64) -> f64 {
    21.4 * libm::log10(4.37 * f / 1000.0 + 1.0)
}

/// Inverse of [`erb_number`]: the frequency in Hz with ERB-number `x`.
#[must_use]
pub fn erb_number_inverse(x: f64) -> f64 {
    (libm::pow(10.0, x / 21.4) - 1.0) * 1000.0 / 4.37
}

/// Absolute threshold of hearing in dB SPL at `f` Hz (Moore et al., 1997).
#[must_use]
pub fn absolute_threshold_db(f: f64) -> f64 {
    let x = f / 1000.0;
    3.64 * libm::pow(x, -0.8) - 6.5 * libm::exp(-0.6 * (x - 3.3) * (x - 3.3))
        + 1e-3 * libm::pow(x, 4.0)
}

/// Excitation (power) that a tone of power `power` at `f` produces at the
/// auditory-filter channel centred on `fc`.
fn excitation(power: f64, f: f64, fc: f64) -> f64 {
    let p = 4.0 * fc / erb_hz(fc);
    let g = (f - fc).abs() / fc;
    power * (1.0 + p * g) * libm::exp(-p * g)
}

/// Unscaled loudness of a tone list: the integral of the specific loudness.
fn raw_loudness(tones: &[(f64, f64)]) -> f64 {
    let lo = erb_number(50.0);
    let hi = erb_number(15_000.0);
    let channels = libm::floor((hi - lo) / ERB_STEP) as usize;
    let mut total = 0.0;
    for i in 0..=channels {
        let fc = erb_number_inverse(lo + i as f64 * ERB_STEP);
        let mut e = 0.0;
        for &(f, level) in tones {
            e += excitation(libm::pow(10.0, level / 10.0), f, fc);
        }
        let e_thr = libm::pow(10.0, absolute_threshold_db(fc) / 10.0);
        let specific = C_SCALE * (libm::pow(e + e_thr, ALPHA) - libm::pow(e_thr, ALPHA));
        total += specific * ERB_STEP;
    }
    total
}

/// Loudness in sones of a list of pure tones `(frequency_hz, level_db_spl)`.
///
/// The result is scaled so that one 1 kHz tone at 40 dB SPL gives 1 sone.
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] for a non-positive or non-finite
/// frequency, or a non-finite level.
pub fn moore_loudness_tones(tones: &[(f64, f64)]) -> Result<f64, AcousticError> {
    if tones
        .iter()
        .any(|&(f, l)| !(f.is_finite() && f > 0.0 && l.is_finite()))
    {
        return Err(AcousticError::InvalidArgument(
            "tone frequencies must be positive and levels finite",
        ));
    }
    let anchor = raw_loudness(&[(1000.0, 40.0)]);
    let raw = raw_loudness(tones);
    Ok(raw / anchor)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn erb_scale_round_trips() {
        for f in [50.0, 500.0, 1000.0, 8000.0] {
            let back = erb_number_inverse(erb_number(f));
            assert!((back - f).abs() < 1e-6 * f, "{f} -> {back}");
        }
    }

    #[test]
    fn threshold_is_most_sensitive_near_three_kilohertz() {
        assert!(absolute_threshold_db(3500.0) < absolute_threshold_db(1000.0));
        assert!(absolute_threshold_db(100.0) > absolute_threshold_db(1000.0));
    }

    #[test]
    fn anchor_is_one_sone() {
        let n = moore_loudness_tones(&[(1000.0, 40.0)]).unwrap();
        assert!((n - 1.0).abs() < 1e-12, "{n}");
    }

    #[test]
    fn loudness_grows_with_level() {
        let at = |db: f64| moore_loudness_tones(&[(1000.0, db)]).unwrap();
        assert!(at(30.0) < at(40.0));
        assert!(at(40.0) < at(50.0));
        assert!(at(50.0) < at(60.0));
    }

    #[test]
    fn nearby_tones_mask_each_other() {
        // Two tones inside one filter share excitation, so the pair is louder
        // than one tone but not twice as loud.
        let one = moore_loudness_tones(&[(1000.0, 40.0)]).unwrap();
        let two = moore_loudness_tones(&[(1000.0, 40.0), (1020.0, 40.0)]).unwrap();
        assert!(two > one && two < 2.0 * one, "one {one} two {two}");
    }

    #[test]
    fn silence_and_bad_input() {
        assert_eq!(moore_loudness_tones(&[]).unwrap(), 0.0);
        assert!(moore_loudness_tones(&[(0.0, 40.0)]).is_err());
        assert!(moore_loudness_tones(&[(1000.0, f64::NAN)]).is_err());
    }
}
