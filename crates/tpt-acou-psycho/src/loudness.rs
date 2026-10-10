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
//! licensed content). Spreading uses Zwicker's published slopes; a single
//! 1 kHz / 40 dB tone rates about 1.02 sone against 1 sone in the standard.
//! The test suite pins that value and checks that loudness grows with level.

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

/// Attenuation in dB of a band of level `level_db` and centre `center_hz`
/// when its energy is heard `dz` Barks away (positive `dz` is higher in
/// frequency). Zwicker's slopes: 27 dB/Bark toward lower frequencies, and
/// `24 + 230/f − 0.2·L` dB/Bark toward higher ones, which gets shallower at
/// high level. The upper slope is clamped at zero.
fn spreading_attenuation_db(level_db: f64, center_hz: f64, dz: f64) -> f64 {
    if dz == 0.0 {
        0.0
    } else if dz < 0.0 {
        -27.0 * dz
    } else {
        let slope = (24.0 + 230.0 / center_hz - 0.2 * level_db).max(0.0);
        slope * dz
    }
}

/// Simplified Zwicker-style stationary loudness (sone) from band sound
/// pressure levels, with spectral spreading.
///
/// `band_levels_db` holds SPL values (dB re 20 µPa) for each of the 24
/// [`BARK_CENTERS_HZ`]. Missing bands (slice shorter than 24) and non-finite
/// entries (use `f64::NEG_INFINITY` for true silence, not 0 dB) contribute
/// nothing. Each band's excitation is the spread sum of every band's energy:
///
/// ```text
/// E_z = Σ_j 10^((L_j − A(z, j)) / 10)
/// N'(z) = 0.08 · max(E_z/E_thr(z) − 1, 0)^0.23   [sone/Bark]
/// ```
///
/// with `A` the spreading attenuation and `E_thr` the Terhardt
/// threshold-in-quiet excitation at the band centre. Total loudness
/// `N = Σ_z N'(z)` (one Bark per band). A 1 kHz tone at 40 dB SPL rates
/// about 1.02 sone, against 1 sone in ISO 532-1.
///
/// Negative SPL values are accepted; they fall below threshold.
///
/// # Errors
/// Currently infallible; the `Result` is kept for API stability.
pub fn zwicker_stationary(band_levels_db: &[f64]) -> Result<f64, tpt_acou_core::AcousticError> {
    let n = band_levels_db.len().min(BARK_CENTERS_HZ.len());
    let mut total = 0.0;
    for z in 0..BARK_CENTERS_HZ.len() {
        let lthr = threshold_in_quiet_db(BARK_CENTERS_HZ[z]);
        let ethr = libm::pow(10.0, lthr / 10.0);
        let mut excitation = 0.0;
        for j in 0..n {
            let lj = band_levels_db[j];
            if !lj.is_finite() {
                continue;
            }
            let dz = z as f64 - j as f64;
            let attenuation = spreading_attenuation_db(lj, BARK_CENTERS_HZ[j], dz);
            excitation += libm::pow(10.0, (lj - attenuation) / 10.0);
        }
        if excitation > ethr {
            total += 0.08 * libm::pow(excitation / ethr - 1.0, 0.23);
        }
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
        // A 1 kHz tone at 40 dB SPL. With spreading the excitation reaches
        // the neighbouring bands, so the model gives ~1.02 sone against the
        // standard's 1 sone. Tolerance is 5 %.
        // Other bands are set to -inf (true silence): 0 dB SPL is NOT silence,
        // since several bands have threshold-in-quiet below 0 dB.
        let mut levels = vec![f64::NEG_INFINITY; 24];
        levels[8] = 40.0; // band 9 centers at 1 kHz
        let n = zwicker_stationary(&levels).unwrap();
        assert!((n - 1.0).abs() < 0.05, "N = {n:.3} sone");
    }

    #[test]
    fn spreading_slopes_match_zwicker() {
        // 27 dB/Bark toward lower frequencies, level-dependent toward higher.
        assert_eq!(spreading_attenuation_db(40.0, 1000.0, 0.0), 0.0);
        assert!((spreading_attenuation_db(40.0, 1000.0, -1.0) - 27.0).abs() < 1e-12);
        let upper = spreading_attenuation_db(40.0, 1000.0, 1.0);
        assert!(
            (upper - (24.0 + 0.23 - 8.0)).abs() < 1e-9,
            "upper slope {upper}"
        );
        // Louder bands spread less steeply upward, and the slope never goes negative.
        assert!(spreading_attenuation_db(80.0, 1000.0, 1.0) < upper);
        assert!(spreading_attenuation_db(400.0, 1000.0, 1.0) >= 0.0);
    }

    #[test]
    fn spreading_makes_neighbours_audible() {
        // Two adjacent equal bands are louder than one band alone, but less
        // than twice as loud: the second band's energy is partly shared.
        let mut one = vec![f64::NEG_INFINITY; 24];
        one[8] = 40.0;
        let mut two = one.clone();
        two[9] = 40.0;
        let n1 = zwicker_stationary(&one).unwrap();
        let n2 = zwicker_stationary(&two).unwrap();
        assert!(n2 > n1 && n2 < 2.0 * n1, "one {n1} two {n2}");
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
