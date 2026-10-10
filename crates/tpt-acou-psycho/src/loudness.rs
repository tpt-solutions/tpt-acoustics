//! Loudness: phon/sone conversion and a simplified Bark-domain stationary
//! model.
//!
//! * [`phon_to_sone`] / [`sone_to_phon`] — the ISO 532-1 mapping between
//!   loudness level (phon, anchored to a 1 kHz tone) and loudness (sone):
//!   `N = 2^((Ln - 40)/10)`.
//! * [`zwicker_stationary`] — a *simplified* Zwicker-style model: band input
//!   levels are mapped onto the 24 Zwicker critical bands, threshold-in-quiet
//!   comes from `tpt-dsp-core`'s Terhardt formula, and specific loudness
//!   follows the standard power law `N' = 0.08 (E/E_thr)^0.23` (sone/Bark)
//!   summed over bands.
//!
//! This is an engineering approximation, *not* a certified ISO 532-1
//! implementation (the standard's exact threshold and spreading tables are
//! licensed content). It has no spectral spreading yet, so a single-band tone
//! rates ~0.57 sone at 1 kHz / 40 dB rather than the 1 sone of the standard;
//! the test suite pins that value and the model is monotone in level.

/// Threshold in quiet (dB SPL) of a pure tone, delegated to
/// `tpt_dsp_core::psychoacoustic::threshold_in_quiet` (Terhardt's analytic
/// approximation: minimum ~2-5 kHz, rising sharply toward the band edges).
#[must_use]
pub fn threshold_in_quiet_db(center_hz: f64) -> f64 {
    tpt_dsp_core::threshold_in_quiet(center_hz)
}

/// Converts loudness level (phon) to loudness (sone): `N = 2^((Ln-40)/10)`.
#[must_use]
pub fn phon_to_sone(phon: f64) -> f64 {
    libm::pow(2.0, (phon - 40.0) / 10.0)
}

/// Converts loudness (sone) to loudness level (phon): `Ln = 40 + 10 log2 N`.
#[must_use]
pub fn sone_to_phon(sone: f64) -> f64 {
    40.0 + 10.0 * libm::log2(sone.max(1e-12))
}

/// The 24 Zwicker critical-band center frequencies (Hz).
pub const BARK_CENTERS_HZ: [f64; 24] = [
    50.0, 150.0, 250.0, 350.0, 450.0, 570.0, 700.0, 840.0, 1000.0, 1170.0, 1370.0, 1600.0, 1850.0,
    2150.0, 2500.0, 2900.0, 3400.0, 4000.0, 4800.0, 5800.0, 7000.0, 8500.0, 10500.0, 13500.0,
];

/// Simplified Zwicker-style stationary loudness (sone) from band sound
/// pressure levels.
///
/// `band_levels_db` holds SPL values (dB re 20 µPa) for each of the 24
/// [`BARK_CENTERS_HZ`]; missing bands (slice shorter than 24) contribute
/// nothing. A band at or below its threshold contributes nothing, so use
/// `f64::NEG_INFINITY` (not 0 dB) to mark a band as silent. Specific loudness
/// per band:
///
/// ```text
/// N'(z) = 0.08 * max(E_z/E_thr(z) - 1, 0)^0.23   [sone/Bark]
/// ```
///
/// with `E = 10^(L/10)` the band excitation proxy and `E_thr` the
/// Terhardt threshold-in-quiet excitation at the band center. Total loudness
/// `N = sum N'(z)` (dz = 1 Bark each). A 1 kHz tone at 40 dB SPL rates about
/// 0.57 sone under this simplified model (no spreading).
///
/// Negative SPL values are accepted; they simply fall below threshold.
///
/// # Errors
/// Currently infallible; the `Result` is kept for API stability.
pub fn zwicker_stationary(band_levels_db: &[f64]) -> Result<f64, tpt_acou_core::AcousticError> {
    let mut total = 0.0;
    for (z, &level) in band_levels_db.iter().enumerate() {
        if z >= BARK_CENTERS_HZ.len() {
            break;
        }
        let lthr = threshold_in_quiet_db(BARK_CENTERS_HZ[z]);
        if level <= lthr {
            continue;
        }
        let e = libm::pow(10.0, level / 10.0);
        let ethr = libm::pow(10.0, lthr / 10.0);
        total += 0.08 * libm::pow(e / ethr - 1.0, 0.23);
    }
    Ok(total)
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn phon_sone_anchor_and_doubling() {
        assert!(approx(phon_to_sone(40.0), 1.0, 1e-12));
        assert!(approx(phon_to_sone(50.0), 2.0, 1e-12));
        assert!(approx(phon_to_sone(30.0), 0.5, 1e-12));
        // Round trips.
        for &p in &[20.0, 40.0, 63.0, 90.0] {
            assert!(approx(sone_to_phon(phon_to_sone(p)), p, 1e-9));
        }
    }

    #[test]
    fn zwicker_one_khz_tone_near_one_sone() {
        // A 1 kHz tone at 40 dB SPL: E/E_thr ~ 5000 gives 0.08 * 5000^0.23
        // ~ 0.57 sone. The ISO 532-1 value is 1 sone; the gap is the missing
        // spreading, so this pins the model's own value rather than the
        // standard's.
        // Other bands are set to -inf (true silence): 0 dB SPL is NOT silence,
        // since several bands have threshold-in-quiet below 0 dB.
        let mut levels = vec![f64::NEG_INFINITY; 24];
        levels[8] = 40.0; // band 9 centers at 1 kHz
        let n = zwicker_stationary(&levels).unwrap();
        assert!(n > 0.5 && n < 0.6, "N = {n:.3} sone");
    }

    #[test]
    fn zwicker_monotone_in_level() {
        let quiet = vec![30.0_f64; 24];
        let loud = vec![70.0_f64; 24];
        let nq = zwicker_stationary(&quiet).unwrap();
        let nl = zwicker_stationary(&loud).unwrap();
        assert!(nl > nq, "quiet {nq} loud {nl}");
        // Doubling per +10 dB holds approximately across the whole spectrum.
        assert!((nl / nq).abs() > 1.5, "loudness ratio {}", nl / nq);
    }

    #[test]
    fn threshold_in_quiet_shape() {
        // Minimum near 3-4 kHz (most sensitive), rising sharply at LF/HF.
        let l_low = threshold_in_quiet_db(50.0);
        let l_mid = threshold_in_quiet_db(3500.0);
        let l_high = threshold_in_quiet_db(13500.0);
        assert!(l_low > l_mid + 5.0);
        assert!(l_high > l_mid + 5.0);
    }

    #[test]
    fn accepts_negative_levels() {
        // Negative SPL is physically valid (very quiet signals); bands below
        // threshold contribute nothing, so an all-negative input is silent.
        let levels = vec![-10.0_f64; 24];
        assert_eq!(zwicker_stationary(&levels).unwrap(), 0.0);
    }
}
