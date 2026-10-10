//! Simultaneous masking thresholds.
//!
//! The spreading function is the classic psychoacoustic-model form
//! (Johnston/Zwicker family), evaluated in Bark:
//!
//! ```text
//! SF(dz) = 15.81 + 7.5 (dz + 0.474) - 17.5 sqrt(1 + (dz + 0.474)^2)   [dB]
//! ```
//!
//! with `dz` the Bark distance from the masker (negative below). A masker of
//! level `Lm` masks down to `Lm + SF(dz) - offset`, where the offset depends
//! on the masker type (noise maskers are ~2 dB more masking than tonal ones
//! in the classic data; both offsets are parameters here).

/// Johnston/Zwicker spreading function value (dB) at Bark distance `dz`
/// (negative = below the masker). Peaks at `dz ≈ -0.474` Bark; delegated to
/// `tpt_dsp_core::psychoacoustic::spreading_function_db` (same constants).
#[must_use]
pub fn spreading_db(dz: f64) -> f64 {
    tpt_dsp_core::spreading_function_db(dz)
}

/// Bark↔Hz mapping, delegated to
/// `tpt_dsp_core::bark::hz_to_bark` (Traunmüller 1990 formula with the
/// standard low/high-frequency corrections).
#[must_use]
pub fn hz_to_bark(f: f64) -> f64 {
    tpt_dsp_core::hz_to_bark(f)
}

/// Masker description for [`masking_threshold_db`].
#[derive(Debug, Clone, Copy)]
pub struct Masker {
    /// Center frequency (Hz).
    pub frequency_hz: f64,
    /// Sound pressure level (dB SPL).
    pub level_db: f64,
    /// Masking offset (dB): how far below its own level the masker masks
    /// (≈2 dB for narrowband noise, ≈6 dB for tonal maskers).
    pub offset_db: f64,
}

/// Simultaneous masking threshold (dB SPL) at frequency `f` from a set of
/// maskers: the maximum over maskers of `level + SF(dz) - offset`, floored at
/// 0 dB SPL (absolute-hearing floor handled by callers via threshold-in-quiet
/// if desired).
#[must_use]
pub fn masking_threshold_db(maskers: &[Masker], f: f64) -> f64 {
    let z = hz_to_bark(f);
    let mut thr = 0.0_f64;
    for m in maskers {
        let dz = z - hz_to_bark(m.frequency_hz);
        thr = thr.max(m.level_db + spreading_db(dz) - m.offset_db);
    }
    thr
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn spreading_function_shape() {
        // Peaks ~0 dB near -0.474 Bark, decays both ways, asymmetric
        // (falls faster above than below at moderate distances).
        // Peak value is -1.69 dB at dz = -0.474 (Johnston's constants).
        assert!(spreading_db(-0.474) > -3.0 && spreading_db(-0.474) <= 0.0);
        // Masking spreads further ABOVE the masker frequency (gentler upper
        // slope): at equal Bark distance the upper side sits higher.
        assert!(spreading_db(1.0) > spreading_db(-1.0));
        assert!(spreading_db(-4.0) < -60.0);
        assert!(spreading_db(4.0) < -30.0);
    }

    #[test]
    fn threshold_decays_with_distance() {
        let masker = Masker {
            frequency_hz: 1000.0,
            level_db: 80.0,
            offset_db: 2.0,
        };
        let near = masking_threshold_db(&[masker], 1100.0);
        let far = masking_threshold_db(&[masker], 4000.0);
        assert!(near > far, "near {near:.1} far {far:.1}");
        assert!(near > 50.0 && far < 60.0);
    }

    #[test]
    fn stronger_masker_raises_threshold() {
        let weak = Masker {
            frequency_hz: 1000.0,
            level_db: 40.0,
            offset_db: 2.0,
        };
        let strong = Masker {
            frequency_hz: 1000.0,
            level_db: 80.0,
            offset_db: 2.0,
        };
        let f = 1500.0;
        assert!(masking_threshold_db(&[strong], f) > masking_threshold_db(&[weak], f));
    }

    #[test]
    fn multiple_maskers_take_maximum() {
        let m1 = Masker {
            frequency_hz: 1000.0,
            level_db: 60.0,
            offset_db: 2.0,
        };
        let m2 = Masker {
            frequency_hz: 4000.0,
            level_db: 90.0,
            offset_db: 2.0,
        };
        let f = 3000.0;
        let both = masking_threshold_db(&[m1, m2], f);
        let only_m1 = masking_threshold_db(&[m1], f);
        let only_m2 = masking_threshold_db(&[m2], f);
        assert!(approx(both, only_m1.max(only_m2), 1e-12));
        assert!(both > 60.0);
    }

    #[test]
    fn bark_mapping_monotone_and_sane() {
        // Traunmüller maps 0 Hz to -0.15 Bark (the low-frequency correction
        // does not fully clamp it) — accepted; the scale is monotone and
        // lands in the canonical ranges at 1 kHz and 16 kHz.
        assert!(hz_to_bark(0.0) > -0.5 && hz_to_bark(0.0) < 0.0);
        assert!(hz_to_bark(1000.0) > 8.0 && hz_to_bark(1000.0) < 9.5);
        assert!(hz_to_bark(16000.0) > 22.0 && hz_to_bark(16000.0) < 25.0);
        assert!(hz_to_bark(4000.0) > hz_to_bark(2000.0));
    }
}
