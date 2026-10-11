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

/// A single-number rating defined by a reference contour and two deficit
/// limits. The rating is the highest integer `R` such that the contour, shifted
/// so its reference value is `R`, has total deficit (`contour − TL`, positive
/// parts only) at most `max_total_deficit_db` and no single deficit above
/// `max_single_deficit_db`.
///
/// The contour values are inputs. [`RatingContour::astm_e413_shape`] is a
/// built-in contour with the E413 shape and limits. Its values are not checked
/// against the published standard. If you have the standard, pass its table
/// with [`RatingContour::new`], and label the result with the table's source.
#[derive(Debug, Clone, PartialEq)]
pub struct RatingContour {
    /// Band centre frequencies in Hz, strictly increasing.
    pub bands_hz: Vec<f64>,
    /// Contour in dB relative to the rating, one value per band.
    pub relative_db: Vec<f64>,
    /// Largest allowed sum of deficits, in dB.
    pub max_total_deficit_db: f64,
    /// Largest allowed deficit in any one band, in dB.
    pub max_single_deficit_db: f64,
    /// Where the contour values came from. Reports should quote this.
    pub source: String,
}

impl RatingContour {
    /// Builds a contour from user-supplied values.
    ///
    /// # Errors
    /// Returns [`AcousticError::InvalidArgument`] if the lists are empty or
    /// differ in length, a frequency is not positive and strictly increasing, a
    /// value or limit is not finite, a limit is not positive, or `source` is
    /// empty.
    pub fn new(
        bands_hz: Vec<f64>,
        relative_db: Vec<f64>,
        max_total_deficit_db: f64,
        max_single_deficit_db: f64,
        source: &str,
    ) -> Result<Self, AcousticError> {
        if bands_hz.is_empty() || bands_hz.len() != relative_db.len() {
            return Err(AcousticError::InvalidArgument(
                "contour needs one value per band, and at least one band",
            ));
        }
        if bands_hz.iter().any(|f| !(f.is_finite() && *f > 0.0))
            || bands_hz.windows(2).any(|w| w[0] >= w[1])
        {
            return Err(AcousticError::InvalidArgument(
                "contour bands must be positive and strictly increasing",
            ));
        }
        if relative_db.iter().any(|v| !v.is_finite()) {
            return Err(AcousticError::InvalidArgument(
                "contour values must be finite",
            ));
        }
        if !(max_total_deficit_db.is_finite()
            && max_total_deficit_db > 0.0
            && max_single_deficit_db.is_finite()
            && max_single_deficit_db > 0.0)
        {
            return Err(AcousticError::InvalidArgument(
                "deficit limits must be positive and finite",
            ));
        }
        if source.trim().is_empty() {
            return Err(AcousticError::InvalidArgument(
                "a contour needs a source description",
            ));
        }
        Ok(Self {
            bands_hz,
            relative_db,
            max_total_deficit_db,
            max_single_deficit_db,
            source: source.to_string(),
        })
    }

    /// The built-in E413-shaped contour at the 16 STC bands, with the E413
    /// limits (32 dB total, 8 dB single). See the type docs for its status.
    #[must_use]
    pub fn astm_e413_shape() -> Self {
        Self {
            bands_hz: STC_BAND_HZ.to_vec(),
            relative_db: STC_CONTOUR_REL_DB.to_vec(),
            max_total_deficit_db: 32.0,
            max_single_deficit_db: 8.0,
            source: "built-in E413-shaped contour (values not checked against the standard)"
                .to_string(),
        }
    }
}

/// Rating of a TL spectrum against a contour. `spectrum` holds
/// `(frequency_hz, tl_db)` pairs; each contour band must be present, matched to
/// a 1 % frequency tolerance. The search runs over integer ratings from −60 to
/// 150.
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] if a band is missing, a TL value
/// is not finite, or no rating in the search range meets the limits.
pub fn rating_against_contour(
    spectrum: &[(f64, f64)],
    contour: &RatingContour,
) -> Result<i32, AcousticError> {
    let mut tl = Vec::with_capacity(contour.bands_hz.len());
    for &f in &contour.bands_hz {
        let got = spectrum
            .iter()
            .find(|&&(ff, _)| (ff - f).abs() < 0.01 * f)
            .map(|&(_, tl)| tl)
            .ok_or(AcousticError::InvalidArgument("missing contour band"))?;
        if !got.is_finite() {
            return Err(AcousticError::InvalidArgument("TL values must be finite"));
        }
        tl.push(got);
    }
    // Deficits grow monotonically with the rating, so the feasible set is
    // downward-closed and a scan from the top finds the maximum.
    for rating in (-60..=150).rev() {
        let mut total = 0.0;
        let mut max_single = 0.0_f64;
        for (b, &rel) in contour.relative_db.iter().enumerate() {
            let deficit = rel + f64::from(rating) - tl[b];
            if deficit > 0.0 {
                total += deficit;
                max_single = max_single.max(deficit);
            }
        }
        if total <= contour.max_total_deficit_db && max_single <= contour.max_single_deficit_db {
            return Ok(rating);
        }
    }
    Err(AcousticError::InvalidArgument(
        "spectrum fails the contour criteria at every rating",
    ))
}

/// Slope of the TII reference contour in dB per octave. It is a design choice
/// of this project, close to the diffuse mass-law slope (about 5 dB per
/// octave above coincidence). It is not taken from any standard.
pub const TII_SLOPE_DB_PER_OCTAVE: f64 = 5.0;

/// Band centre at which the TII contour equals its rating. This is the
/// reference frequency of the rating, as 500 Hz is for STC.
pub const TII_REFERENCE_HZ: f64 = 500.0;

/// TPT Isolation Index v1: an open, single-number rating defined in
/// `docs/open-rating.md`. It is not ASTM E413 or any other standard, and its
/// results must not be reported as such.
///
/// The reference contour is a straight line on a log-frequency axis, equal to
/// the rating at [`TII_REFERENCE_HZ`] and rising at [`TII_SLOPE_DB_PER_OCTAVE`]
/// dB per octave, at the 16 STC bands. The deficit limits are 32 dB total and
/// 8 dB single, the same rule as the rating procedure above.
///
/// # Errors
/// As [`rating_against_contour`].
pub fn tii_rating(spectrum: &[(f64, f64)]) -> Result<i32, AcousticError> {
    tii_rating_on_bands(spectrum, &STC_BAND_HZ)
}

/// Bands of TII version 2: the 16 STC bands extended to 63 Hz, 80 Hz, 100 Hz
/// and 5 kHz, 20 in all. The version 2 limits scale with the band count (see
/// [`tii_v2_rating`]).
pub const TII_V2_BANDS_HZ: [f64; 20] = [
    63.0, 80.0, 100.0, 125.0, 160.0, 200.0, 250.0, 315.0, 400.0, 500.0, 630.0, 800.0, 1000.0,
    1250.0, 1600.0, 2000.0, 2500.0, 3150.0, 4000.0, 5000.0,
];

/// TPT Isolation Index v2: the same line contour and deficit rule as v1, on the
/// 20 bands [`TII_V2_BANDS_HZ`]. The total-deficit limit is 2 dB per band
/// (40 dB), so a line spectrum still gives rating `level + 2`, as in v1.
///
/// Below 125 Hz the diffuse mass law, and the 5 dB/octave line, are weak
/// models of real partitions, which have stiffness and panel resonances there.
/// v2 is for comparison with v1 and for showing low-frequency sensitivity. It
/// is not a replacement for v1. See `docs/open-rating.md`.
///
/// # Errors
/// As [`rating_against_contour`].
pub fn tii_v2_rating(spectrum: &[(f64, f64)]) -> Result<i32, AcousticError> {
    tii_rating_on_bands(spectrum, &TII_V2_BANDS_HZ)
}

/// The TII line contour on `bands`, with the total-deficit limit at 2 dB per
/// band and the single-deficit limit at 8 dB.
fn tii_rating_on_bands(spectrum: &[(f64, f64)], bands: &[f64]) -> Result<i32, AcousticError> {
    let relative_db = bands
        .iter()
        .map(|&f| TII_SLOPE_DB_PER_OCTAVE * libm::log2(f / TII_REFERENCE_HZ))
        .collect();
    let contour = RatingContour::new(
        bands.to_vec(),
        relative_db,
        2.0 * bands.len() as f64,
        8.0,
        "TPT Isolation Index (project definition, not a standard)",
    )?;
    rating_against_contour(spectrum, &contour)
}

/// Sound Transmission Class of a 1/3-octave TL spectrum with the built-in
/// E413-shaped contour. See [`RatingContour::astm_e413_shape`] for the status
/// of its values, and [`rating_against_contour`] to use another contour.
///
/// # Errors
/// As [`rating_against_contour`].
pub fn sound_transmission_class(spectrum: &[(f64, f64)]) -> Result<i32, AcousticError> {
    rating_against_contour(spectrum, &RatingContour::astm_e413_shape())
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

    fn spectrum_from(tl: &[f64]) -> Vec<(f64, f64)> {
        STC_BAND_HZ
            .iter()
            .copied()
            .zip(tl.iter().copied())
            .collect()
    }

    #[test]
    fn user_contour_equal_to_builtin_gives_the_builtin_rating() {
        let contour = RatingContour::new(
            STC_BAND_HZ.to_vec(),
            STC_CONTOUR_REL_DB.to_vec(),
            32.0,
            8.0,
            "test copy of the built-in values",
        )
        .unwrap();
        for tl in [
            [30.0; 16],
            [45.0; 16],
            [
                20.0, 22.0, 25.0, 28.0, 31.0, 34.0, 38.0, 41.0, 44.0, 46.0, 48.0, 50.0, 51.0, 52.0,
                52.0, 52.0,
            ],
        ] {
            let spec = spectrum_from(&tl);
            assert_eq!(
                rating_against_contour(&spec, &contour).unwrap(),
                sound_transmission_class(&spec).unwrap()
            );
        }
    }

    #[test]
    fn shifting_the_contour_shifts_the_rating_by_the_same_amount() {
        let spec = spectrum_from(&[
            35.0, 36.0, 38.0, 40.0, 42.0, 44.0, 46.0, 47.0, 48.0, 49.0, 50.0, 50.0, 51.0, 52.0,
            52.0, 52.0,
        ]);
        let base = RatingContour::astm_e413_shape();
        let mut raised = base.clone();
        for v in raised.relative_db.iter_mut() {
            *v += 3.0;
        }
        let r0 = rating_against_contour(&spec, &base).unwrap();
        let r1 = rating_against_contour(&spec, &raised).unwrap();
        // A contour 3 dB higher everywhere is 3 dB harder to meet at the same rating.
        assert_eq!(r1, r0 - 3);
    }

    #[test]
    fn looser_limits_never_lower_the_rating() {
        let spec = spectrum_from(&[
            20.0, 22.0, 25.0, 28.0, 31.0, 34.0, 38.0, 41.0, 44.0, 46.0, 48.0, 50.0, 51.0, 52.0,
            52.0, 52.0,
        ]);
        let strict = RatingContour::astm_e413_shape();
        let mut loose = strict.clone();
        loose.max_total_deficit_db = 64.0;
        loose.max_single_deficit_db = 16.0;
        let a = rating_against_contour(&spec, &strict).unwrap();
        let b = rating_against_contour(&spec, &loose).unwrap();
        assert!(b >= a, "{b} < {a}");
        assert!(b > a, "the looser limits should admit a higher rating here");
    }

    #[test]
    fn tii_of_a_line_spectrum_is_its_level_plus_the_allowed_excess() {
        // TL on the TII line itself, level `a` at 500 Hz. Each band can sit
        // 2 dB above the contour before the total deficit (16 bands) passes
        // 32 dB, so the rating is a + 2 exactly.
        for a in [20.0, 30.0, 41.0] {
            let spec: Vec<(f64, f64)> = STC_BAND_HZ
                .iter()
                .map(|&f| {
                    (
                        f,
                        a + TII_SLOPE_DB_PER_OCTAVE * libm::log2(f / TII_REFERENCE_HZ),
                    )
                })
                .collect();
            assert_eq!(tii_rating(&spec).unwrap(), a as i32 + 2);
        }
    }

    #[test]
    fn tii_follows_the_diffuse_mass_law() {
        let air = Medium::air_20c();
        let spectrum = |m: f64| -> Vec<(f64, f64)> {
            STC_BAND_HZ
                .iter()
                .map(|&f| (f, mass_law_diffuse_tl(m, f, &air)))
                .collect()
        };
        let light = tii_rating(&spectrum(5.0)).unwrap();
        let heavy = tii_rating(&spectrum(10.0)).unwrap();
        // The mass law adds about 6 dB per doubling of mass; the TII slope is
        // 5 dB per octave, so the rating step is close to but not exactly 6.
        let step = heavy - light;
        assert!(
            (5..=7).contains(&step),
            "step {step} from {light} to {heavy}"
        );
        assert!(heavy > light);
    }

    #[test]
    fn tii_v2_of_a_line_spectrum_is_its_level_plus_two() {
        for a in [20.0, 30.0, 41.0] {
            let spec: Vec<(f64, f64)> = TII_V2_BANDS_HZ
                .iter()
                .map(|&f| {
                    (
                        f,
                        a + TII_SLOPE_DB_PER_OCTAVE * libm::log2(f / TII_REFERENCE_HZ),
                    )
                })
                .collect();
            assert_eq!(tii_v2_rating(&spec).unwrap(), a as i32 + 2);
        }
    }

    #[test]
    fn tii_v2_follows_the_mass_law_like_v1() {
        let air = Medium::air_20c();
        let spectrum = |m: f64| -> Vec<(f64, f64)> {
            TII_V2_BANDS_HZ
                .iter()
                .map(|&f| (f, mass_law_diffuse_tl(m, f, &air)))
                .collect()
        };
        let step = tii_v2_rating(&spectrum(10.0)).unwrap() - tii_v2_rating(&spectrum(5.0)).unwrap();
        assert!((5..=7).contains(&step), "step {step}");
    }

    #[test]
    fn tii_v2_needs_every_band() {
        assert!(tii_v2_rating(&[(125.0, 30.0)]).is_err());
    }

    #[test]
    fn tii_shifts_one_for_one_with_the_level() {
        let base = spectrum_from(&[
            35.0, 36.0, 38.0, 40.0, 42.0, 44.0, 46.0, 47.0, 48.0, 49.0, 50.0, 50.0, 51.0, 52.0,
            52.0, 52.0,
        ]);
        let raised: Vec<(f64, f64)> = base.iter().map(|&(f, tl)| (f, tl + 3.0)).collect();
        assert_eq!(tii_rating(&raised).unwrap(), tii_rating(&base).unwrap() + 3);
    }

    #[test]
    fn tii_needs_every_band() {
        assert!(tii_rating(&[(125.0, 30.0)]).is_err());
    }

    #[test]
    fn user_contour_rejects_bad_input() {
        let ok = |bands: Vec<f64>, rel: Vec<f64>, t: f64, s: f64, src: &str| {
            RatingContour::new(bands, rel, t, s, src).is_err()
        };
        assert!(ok(vec![], vec![], 32.0, 8.0, "x"));
        assert!(ok(vec![125.0, 250.0], vec![0.0], 32.0, 8.0, "x"));
        assert!(ok(vec![250.0, 125.0], vec![0.0, 0.0], 32.0, 8.0, "x"));
        assert!(ok(vec![125.0], vec![f64::NAN], 32.0, 8.0, "x"));
        assert!(ok(vec![125.0], vec![0.0], 0.0, 8.0, "x"));
        assert!(ok(vec![125.0], vec![0.0], 32.0, 8.0, "  "));
        assert!(RatingContour::new(vec![125.0], vec![0.0], 32.0, 8.0, "provided").is_ok());
    }

    #[test]
    fn user_contour_reports_missing_bands_and_unreachable_ratings() {
        let contour =
            RatingContour::new(vec![125.0, 1000.0], vec![0.0, 0.0], 32.0, 8.0, "t").unwrap();
        assert!(rating_against_contour(&[(125.0, 30.0)], &contour).is_err());
        // With every TL far below the contour the total deficit can never be met.
        let hopeless =
            RatingContour::new(vec![125.0, 1000.0], vec![0.0, 0.0], 0.5, 0.5, "t").unwrap();
        assert!(rating_against_contour(&[(125.0, -200.0), (1000.0, -200.0)], &hopeless).is_err());
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
