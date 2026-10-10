//! Sound transmission loss and the ASTM E413 Sound Transmission Class.
//!
//! * Normal-incidence mass law: `tau = 1/(1 + beta^2)`,
//!   `beta = m*omega / (2 rho c)`.
//! * Diffuse-field mass law: exact integral of the angular mass law over the
//!   hemisphere, `tau_diff = ln(1 + beta^2)/beta^2` — equivalent to the
//!   familiar field-incidence approximation `TL = 20 log10(m f) - 47.5 dB`
//!   (SI units) in the mass-law range.
//! * STC (ASTM E413): 1/3-octave TL spectrum from 125 Hz to 4 kHz compared
//!   against the reference contour (-16, -13, -10, -7, -4, -1, 0 dB relative
//!   to 500 Hz at 125–500 Hz, then +1, +2, +3, +4 at 630–1250 Hz, then flat
//!   +4 dB from 1600–4000 Hz), shifted in 1 dB steps to the highest value
//!   with total deficit ≤ 32 dB and no single deficit > 8 dB.

use tpt_acou_core::error::AcousticError;
use tpt_acou_core::medium::Medium;

/// The ASTM E413 reference contour in dB relative to its 500 Hz value, at the
/// 16 STC 1/3-octave bands in [`STC_BAND_HZ`] (125 Hz … 4000 Hz). The STC
/// rating is the absolute contour value at 500 Hz (index 6).
pub const STC_CONTOUR_REL_DB: [f64; 16] = [
    -16.0, -13.0, -10.0, -7.0, -4.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 4.0, 4.0, 4.0, 4.0, 4.0,
];

/// The STC band center frequencies in Hz, matching [`STC_CONTOUR_REL_DB`].
pub const STC_BAND_HZ: [f64; 16] = [
    125.0, 160.0, 200.0, 250.0, 315.0, 400.0, 500.0, 630.0, 800.0, 1000.0, 1250.0, 1600.0, 2000.0,
    2500.0, 3150.0, 4000.0,
];

/// Normal-incidence mass-law transmission loss of a single leaf of surface
/// mass `m_kg_per_m2` at frequency `f` into `medium`:
/// `TL = 10 log10(1 + beta^2)`, `beta = m omega / (2 rho c)`.
#[must_use]
pub fn mass_law_normal_tl(m_kg_per_m2: f64, f: f64, medium: &Medium) -> f64 {
    let beta = beta_of(m_kg_per_m2, f, medium);
    10.0 * libm::log10(1.0 + beta * beta)
}

/// Exact diffuse-field (random-incidence) mass-law transmission loss:
/// `TL = 10 log10(beta^2 / ln(1 + beta^2))` — the analytic average of the
/// angular mass law over the hemisphere. Limit: 0 dB as `beta -> 0`.
#[must_use]
pub fn mass_law_diffuse_tl(m_kg_per_m2: f64, f: f64, medium: &Medium) -> f64 {
    let beta = beta_of(m_kg_per_m2, f, medium);
    if beta < 1e-9 {
        return 0.0;
    }
    10.0 * libm::log10(beta * beta / libm::log(1.0 + beta * beta))
}

fn beta_of(m_kg_per_m2: f64, f: f64, medium: &Medium) -> f64 {
    let omega = core::f64::consts::TAU * f;
    m_kg_per_m2 * omega / (2.0 * medium.characteristic_impedance())
}

/// ASTM E413 Sound Transmission Class of a 1/3-octave TL spectrum.
///
/// `spectrum` is a slice of `(frequency_hz, tl_db)` pairs; the 16 STC bands
/// (125 Hz … 4 kHz) must all be present (frequencies matched to a 1%
/// tolerance). Returns the highest integer STC value `R` such that the
/// reference contour, positioned so its 500 Hz value is `R`, has every
/// deficit (`contour - TL`, positive parts only) ≤ 8 dB and total deficit
/// ≤ 32 dB.
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] if any STC band is missing, or
/// if no rating in the search range satisfies the criteria.
pub fn sound_transmission_class(spectrum: &[(f64, f64)]) -> Result<i32, AcousticError> {
    // Extract the TL at each STC band.
    let mut tl = [0.0_f64; 16];
    for (band, &f) in STC_BAND_HZ.iter().enumerate() {
        let got = spectrum
            .iter()
            .find(|&&(ff, _)| (ff - f).abs() < 0.01 * f)
            .map(|&(_, tl)| tl)
            .ok_or(AcousticError::InvalidArgument("missing STC band"))?;
        tl[band] = got;
    }
    // Deficits grow monotonically with the rating, so the feasible set is
    // downward-closed and a scan from the top finds the maximum.
    for rating in (-60..=150).rev() {
        let mut total = 0.0;
        let mut max_single = 0.0_f64;
        for b in 0..16 {
            let contour = STC_CONTOUR_REL_DB[b] + f64::from(rating);
            let deficit = contour - tl[b];
            if deficit > 0.0 {
                total += deficit;
                max_single = max_single.max(deficit);
            }
        }
        if total <= 32.0 && max_single <= 8.0 {
            return Ok(rating);
        }
    }
    Err(AcousticError::InvalidArgument(
        "spectrum fails the STC contour criteria at every rating",
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    const AIR: Medium = Medium {
        density: 1.2041,
        speed_of_sound: 343.21,
    };

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn diffuse_mass_law_near_field_law_empirical() {
        // 5 kg/m^2 at 500 Hz: the exact diffuse integral sits ~2.5 dB below
        // the classical field-incidence empirical (20log10(mf) - 47.5),
        // which is calibrated to near-45-degree incidence rather than a
        // true hemisphere average.
        let tl = mass_law_diffuse_tl(5.0, 500.0, &AIR);
        let empirical = 20.0 * libm::log10(5.0 * 500.0) - 47.5;
        assert!(
            (tl - empirical).abs() < 3.5,
            "diffuse {tl:.2} vs empirical {empirical:.2}"
        );
        assert!(tl > 15.0 && tl < 19.0, "diffuse {tl:.2}");
    }

    #[test]
    fn diffuse_is_below_normal_incidence() {
        for f in [200.0_f64, 630.0, 2000.0] {
            let n = mass_law_normal_tl(8.0, f, &AIR);
            let d = mass_law_diffuse_tl(8.0, f, &AIR);
            assert!(d < n && n - d < 15.0, "f={f}: normal {n} diffuse {d}");
        }
    }

    #[test]
    fn mass_law_doubles_per_surface_mass() {
        let a = mass_law_normal_tl(4.0, 1000.0, &AIR);
        let b = mass_law_normal_tl(8.0, 1000.0, &AIR);
        assert!(approx(b - a, 6.0, 0.5), "{a} vs {b}");
    }

    #[test]
    fn stc_of_the_reference_contour_is_52() {
        // Spectrum = the STC-50 contour (500 Hz value 50). At rating R the
        // deficit is (R - 50) in every band: 16 x (R-50) <= 32 gives R = 52.
        let spectrum: Vec<(f64, f64)> = STC_BAND_HZ
            .iter()
            .zip(STC_CONTOUR_REL_DB.iter())
            .map(|(&f, &rel)| (f, 50.0 + rel))
            .collect();
        assert_eq!(sound_transmission_class(&spectrum).unwrap(), 52);
    }

    #[test]
    fn stc_of_uniform_60_db_spectrum() {
        // Flat 60 dB: at rating R the deficits are max(0, rel + R - 60).
        // R = 60 gives 1+2+3+6*4 = 30 dB total (<= 32, max 4 dB); R = 61
        // adds 1 to every band with rel >= 0 for 40 dB total -> STC = 60.
        let spectrum: Vec<(f64, f64)> = STC_BAND_HZ.iter().map(|&f| (f, 60.0)).collect();
        assert_eq!(sound_transmission_class(&spectrum).unwrap(), 60);
    }

    #[test]
    fn stc_contour_matches_astm_e413_shape() {
        // Relative contour: -16 at 125 Hz, 0 at 500 Hz, flat +4 from 1600 Hz.
        assert_eq!(STC_CONTOUR_REL_DB[0], -16.0);
        assert_eq!(STC_CONTOUR_REL_DB[6], 0.0);
        assert_eq!(STC_CONTOUR_REL_DB[10], 4.0);
        assert!(STC_CONTOUR_REL_DB[11..].iter().all(|&v| v == 4.0));
    }

    #[test]
    fn stc_rejects_missing_bands() {
        assert!(sound_transmission_class(&[(125.0, 30.0)]).is_err());
    }
}
