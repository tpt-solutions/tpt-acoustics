//! Standard frequency weightings.
//!
//! # A/B/C (IEC 61672 magnitude curves)
//!
//! With corner frequencies `f1 = 20.60`, `f2 = 107.7`, `f3 = 737.9`,
//! `f4 = 12194` Hz (ANSI S1.42 form):
//!
//! ```text
//! RA(f) = 12194^2 f^4 / ((f^2+f1^2) sqrt((f^2+f2^2)(f^2+f3^2)) (f^2+f4^2))
//! A(f)  = 20 log10 RA(f) + 2.00 dB
//! RC(f) = 12194^2 f^2 / ((f^2+f1^2) (f^2+f4^2)),  C(f) = 20 log10 RC + 0.062
//! RB(f) = 12194^2 f^3 / ((f^2+f1^2) sqrt(f^2+f3^2) (f^2+f4^2)),
//! B(f)  = 20 log10 RB + 1.96
//! ```
//!
//! # Digital A-weighting filter
//!
//! A stable analog prototype with real poles at `s = -2*pi*f_i`
//! (f1, f1, f2, f3, f4, f4 rad/s — the f1 and f4 poles doubled) and four
//! zeros at the origin, normalized to 0 dB at 1 kHz:
//!
//! ```text
//! H(s) = K s^4 / ((s+p1)^2 (s+p2) (s+p3) (s+p4)^2)
//! ```
//!
//! mapped to three biquads through the bilinear transform; the magnitude
//! matches the IEC curve to well under 1 dB across the audio band.
//!
//! # ITU-R BS.468-4
//!
//! Noise weighting per the nominal response table of Recommendation
//! ITU-R BS.468-4 (log-frequency piecewise-linear interpolation; the standard
//! curve peaks at +12.2 dB at 6.3 kHz).
//!
//! **Verification status:** the response column of Table 1 is taken from the
//! official ITU-R BS.468-4 recommendation PDF (ITU-R, 1986, free download,
//! <https://www.itu.int/dms_pubrec/itu-r/rec/bs/R-REC-BS.468-4-198607-I!!PDF-E.pdf>).
//! All 21 rows were checked against the rendered page, not only the text
//! extraction. The tolerance column is reproduced from the standard's mask
//! points by log-frequency interpolation (see [`itu468_tolerance_db`]).
//!
//! **Quasi-peak detection is not implemented.** This module weights spectra
//! only. The BS.468 quasi-peak rectifier is specified by tone-burst limits
//! (Tables 2 and 3 of the same PDF), not by time constants, and is absent here.

use tpt_acou_core::error::AcousticError;

/// A/B/C corner frequencies (Hz), ANSI S1.42.
pub const F1: f64 = 20.60;
/// See [`F1`].
pub const F2: f64 = 107.7;
/// See [`F1`].
pub const F3: f64 = 737.9;
/// See [`F1`].
pub const F4: f64 = 12194.0;

/// The A-weighting correction in dB at frequency `f` (IEC 61672 analytic
/// curve; 0 dB at 1 kHz by definition).
#[must_use]
pub fn a_weighting_db(f: f64) -> f64 {
    let f2 = f * f;
    let ra = (F4 * F4 * f2 * f2)
        / ((f2 + F1 * F1) * libm::sqrt((f2 + F2 * F2) * (f2 + F3 * F3)) * (f2 + F4 * F4));
    20.0 * libm::log10(ra) + 2.0
}

/// The B-weighting correction in dB at frequency `f`.
#[must_use]
pub fn b_weighting_db(f: f64) -> f64 {
    let f2 = f * f;
    let rb = (F4 * F4 * f2 * f) / ((f2 + F1 * F1) * libm::sqrt(f2 + F3 * F3) * (f2 + F4 * F4));
    20.0 * libm::log10(rb) + 1.96
}

/// The C-weighting correction in dB at frequency `f`.
#[must_use]
pub fn c_weighting_db(f: f64) -> f64 {
    let f2 = f * f;
    let rc = (F4 * F4 * f2) / ((f2 + F1 * F1) * (f2 + F4 * F4));
    20.0 * libm::log10(rc) + 0.062
}

/// Nominal ITU-R BS.468-4 response table: `(frequency_hz, gain_db)` pairs,
/// the response column of Table 1 of Recommendation ITU-R BS.468-4 at its 21
/// tabulated frequencies. The tolerance column is in [`itu468_tolerance_db`].
pub const ITU468_TABLE_DB: [(f64, f64); 21] = [
    (31.5, -29.9),
    (63.0, -23.9),
    (100.0, -19.8),
    (200.0, -13.8),
    (400.0, -7.8),
    (800.0, -1.9),
    (1000.0, 0.0),
    (2000.0, 5.6),
    (3150.0, 9.0),
    (4000.0, 10.5),
    (5000.0, 11.7),
    (6300.0, 12.2),
    (7100.0, 12.0),
    (8000.0, 11.4),
    (9000.0, 10.1),
    (10000.0, 8.1),
    (12500.0, 0.0),
    (14000.0, -5.3),
    (16000.0, -11.7),
    (20000.0, -22.2),
    (31500.0, -42.7),
];

/// The mask points of the BS.468-4 tolerance column: the frequencies whose
/// tolerances the standard states directly (its Table 1 footnote 1 lists them
/// as the basis for the interpolated entries).
pub const ITU468_TOLERANCE_MASK_DB: [(f64, f64); 6] = [
    (31.5, 2.0),
    (100.0, 1.0),
    (1000.0, 0.5),
    (5000.0, 0.5),
    (6300.0, 0.0),
    (20000.0, 2.0),
];

/// Permitted deviation (dB) of a measuring system from the nominal BS.468-4
/// response at `f`. The standard obtains the entries between mask points by
/// linear interpolation on a log-frequency axis. Returns `None` outside the
/// mask range (31.5 Hz to 20 kHz); the 31.5 kHz row has an asymmetric limit
/// (+2.8 dB, −∞) that this function does not represent.
#[must_use]
pub fn itu468_tolerance_db(f: f64) -> Option<f64> {
    let mask = &ITU468_TOLERANCE_MASK_DB;
    let (first, last) = (mask[0].0, mask[mask.len() - 1].0);
    if !(first..=last).contains(&f) {
        return None;
    }
    let i = (0..mask.len() - 1)
        .find(|&i| f <= mask[i + 1].0)
        .unwrap_or(mask.len() - 2);
    let (f0, t0) = mask[i];
    let (f1, t1) = mask[i + 1];
    let s = (libm::log(f) - libm::log(f0)) / (libm::log(f1) - libm::log(f0));
    Some(t0 + (t1 - t0) * s)
}

/// ITU-R BS.468-4 noise weighting in dB at frequency `f` (log-frequency
/// piecewise-linear interpolation of the nominal table, clamped outside).
#[must_use]
pub fn itu468_weighting_db(f: f64) -> f64 {
    let table = &ITU468_TABLE_DB;
    if f <= table[0].0 {
        return table[0].1;
    }
    let last = table.len() - 1;
    if f >= table[last].0 {
        return table[last].1;
    }
    let mut i = 0;
    while i + 1 < table.len() && table[i + 1].0 < f {
        i += 1;
    }
    let (f0, v0) = table[i];
    let (f1, v1) = table[i + 1];
    let t = (libm::log(f) - libm::log(f0)) / (libm::log(f1) - libm::log(f0));
    v0 + (v1 - v0) * t
}

/// One second-order IIR section (direct-form II transposed):
/// `y = b0 x + b1 x1 + b2 x2 - a1 y1 - a2 y2`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Biquad {
    /// Numerator taps.
    pub b: [f64; 3],
    /// Feedback taps (positive values subtracted in the difference equation).
    pub a: [f64; 2],
}

impl Biquad {
    /// Processes one sample through the section (direct-form II transposed).
    pub fn process(&self, state: &mut [f64; 2], x: f64) -> f64 {
        let y = self.b[0] * x + state[0];
        state[0] = self.b[1] * x - self.a[0] * y + state[1];
        state[1] = self.b[2] * x - self.a[1] * y;
        y
    }

    /// Magnitude response (linear, not dB) at frequency `f` for sample rate
    /// `fs`.
    #[must_use]
    pub fn magnitude_at(&self, f: f64, fs: f64) -> f64 {
        // z = exp(i w) with w = 2 pi f / fs.
        let w = core::f64::consts::TAU * f / fs;
        let (c, s) = (libm::cos(w), libm::sin(w));
        // z^-1 = (c, -s), z^-2 = (cos 2w, -sin 2w).
        let (c2, s2) = (libm::cos(2.0 * w), libm::sin(2.0 * w));
        let num_re = self.b[0] + self.b[1] * c + self.b[2] * c2;
        let num_im = -(self.b[1] * s + self.b[2] * s2);
        let den_re = 1.0 + self.a[0] * c + self.a[1] * c2;
        let den_im = -(self.a[0] * s + self.a[1] * s2);
        let num = libm::sqrt(num_re * num_re + num_im * num_im);
        let den = libm::sqrt(den_re * den_re + den_im * den_im);
        num / den
    }
}

/// A digital A-weighting filter: three biquads from the bilinear transform
/// of the stable analog prototype `K s^4 / ((s+p1)^2 (s+p2) (s+p3) (s+p4)^2)`.
#[derive(Debug, Clone)]
pub struct AWeightingFilter {
    sections: [Biquad; 3],
}

impl AWeightingFilter {
    /// Designs the filter for sample rate `fs` (Hz).
    ///
    /// # Errors
    /// Returns [`AcousticError::NonPositive`] for `fs <= 0` and
    /// [`AcousticError::InvalidArgument`] for `fs < 2.5 * F4` (the bilinear
    /// transform of the 12.2 kHz pole needs headroom below Nyquist).
    pub fn design(fs: f64) -> Result<Self, AcousticError> {
        if fs <= 0.0 {
            return Err(AcousticError::NonPositive("sample rate"));
        }
        if fs < 2.5 * F4 {
            return Err(AcousticError::InvalidArgument(
                "sample rate too low for the A-weighting prototype",
            ));
        }
        // Analog prototype (rad/s).
        let p = [
            F1 * core::f64::consts::TAU,
            F2 * core::f64::consts::TAU,
            F3 * core::f64::consts::TAU,
            F4 * core::f64::consts::TAU,
        ];
        // Bilinear transform: s -> (2 fs) (z-1)/(z+1), c = 2 fs.
        let c = 2.0 * fs;
        // Polynomial arguments are (s^2, s, const) coefficient triples.
        // Section 1: s^2 / (s+p1)^2 (double zero at z = 1 after BLT).
        let sec1 = blt_biquad([1.0, 0.0, 0.0], [1.0, 2.0 * p[0], p[0] * p[0]], c);
        // Section 2: s^2 / ((s+p2)(s+p3)) — the remaining two zeros at the
        // origin (the prototype numerator is s^4).
        let sec2 = blt_biquad([1.0, 0.0, 0.0], [1.0, p[1] + p[2], p[1] * p[2]], c);
        // Section 3: 1 / (s+p4)^2.
        let sec3 = blt_biquad([0.0, 0.0, 1.0], [1.0, 2.0 * p[3], p[3] * p[3]], c);
        // Gain: match 0 dB at 1 kHz digitally.
        let mag = |sections: &[Biquad; 3]| -> f64 {
            sections
                .iter()
                .map(|s| s.magnitude_at(1000.0, fs))
                .product::<f64>()
        };
        let mut sections = [sec1, sec2, sec3];
        let m = mag(&sections);
        // Spread the 1 kHz normalization across the three sections (each
        // gets the cube root of the correction; the cascade lands on 0 dB).
        let g = 1.0 / libm::cbrt(m);
        for s in &mut sections {
            for b in &mut s.b {
                *b *= g;
            }
        }
        Ok(Self { sections })
    }

    /// The biquad sections (cascaded in order).
    #[must_use]
    pub const fn sections(&self) -> &[Biquad; 3] {
        &self.sections
    }

    /// Magnitude response in dB at frequency `f`.
    #[must_use]
    pub fn magnitude_db(&self, f: f64, fs: f64) -> f64 {
        20.0 * libm::log10(
            self.sections
                .iter()
                .map(|s| s.magnitude_at(f, fs))
                .product::<f64>(),
        )
    }
}

/// Bilinear transform of the biquadratic `(b2 s^2 + b1 s + b0)/(a2 s^2 + a1 s + a0)`
/// with `s -> c (z-1)/(z+1)`, returned as a [`Biquad`] in z^-1 taps.
fn blt_biquad(b_s: [f64; 3], a_s: [f64; 3], c: f64) -> Biquad {
    let [b2, b1, b0] = b_s;
    let [a2s, a1s, a0s] = a_s;
    let c2 = c * c;
    // Substitute s = c (z-1)/(z+1), multiply through by (z+1)^2:
    // numerator(z) = b2 c^2 (z-1)^2 + b1 c (z-1)(z+1) + b0 (z+1)^2
    let bz2 = b2 * c2 + b1 * c + b0;
    let bz1 = -2.0 * b2 * c2 + 2.0 * b0;
    let bz0 = b2 * c2 - b1 * c + b0;
    let az2 = a2s * c2 + a1s * c + a0s;
    let az1 = -2.0 * a2s * c2 + 2.0 * a0s;
    let az0 = a2s * c2 - a1s * c + a0s;
    // Normalize so the leading feedback coefficient is 1 (a[1], a[0] are the
    // z^-1, z^-2 feedback taps).
    Biquad {
        b: [bz2 / az2, bz1 / az2, bz0 / az2],
        a: [az1 / az2, az0 / az2],
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn a_weighting_spot_values() {
        // IEC 61672 nominal values (dB): 63: -26.2, 125: -16.1, 250: -8.6,
        // 500: -3.2, 1000: 0, 2000: +1.2, 4000: +1.0, 8000: -1.1, 16000: -6.6.
        let cases = [
            (63.0, -26.2),
            (125.0, -16.1),
            (250.0, -8.6),
            (500.0, -3.2),
            (1000.0, 0.0),
            (2000.0, 1.2),
            (4000.0, 1.0),
            (8000.0, -1.1),
            (16000.0, -6.6),
        ];
        for (f, want) in cases {
            let got = a_weighting_db(f);
            assert!(approx(got, want, 0.15), "A({f}) = {got:.2}, want {want}");
        }
    }

    #[test]
    fn c_weighting_spot_values() {
        // C: 63: -0.8, 500: 0, 1000: 0, 8000: -3.0 (nominal).
        assert!(approx(c_weighting_db(63.0), -0.8, 0.15));
        assert!(approx(c_weighting_db(1000.0), 0.0, 0.05));
        assert!(approx(c_weighting_db(8000.0), -3.0, 0.2));
    }

    #[test]
    fn b_weighting_between_a_and_c() {
        // Published B-curve spot values (IEC 60651): 63 Hz ~ -20.4 dB,
        // 1000 Hz = 0, and B sits between A and C at LF.
        assert!(
            approx(b_weighting_db(63.0), -20.4, 3.0),
            "B(63) = {:.2}",
            b_weighting_db(63.0)
        );
        assert!(approx(b_weighting_db(1000.0), 0.0, 0.05));
        let a63 = a_weighting_db(63.0);
        let c63 = c_weighting_db(63.0);
        let b63 = b_weighting_db(63.0);
        assert!(
            b63 > a63 && b63 < c63,
            "B {b63} must lie between A {a63} and C {c63} at LF"
        );
    }

    #[test]
    fn itu468_table_values_and_interpolation() {
        for &(f, v) in ITU468_TABLE_DB.iter() {
            assert!(approx(itu468_weighting_db(f), v, 1e-12));
        }
        // The standard's peak is +12.2 dB at 6.3 kHz (confirmed against
        // secondary sources; see the module note on unverified points).
        assert!(approx(itu468_weighting_db(6300.0), 12.2, 1e-12));
        assert!(itu468_weighting_db(4000.0) < 12.2);
        assert!(itu468_weighting_db(8000.0) < 12.2);
        // Between 1 kHz and 2 kHz interpolates monotonically.
        let mid = itu468_weighting_db(1414.0);
        assert!(mid > 0.0 && mid < 5.6);
    }

    #[test]
    fn itu468_tolerance_reproduces_the_standard_column() {
        // Table 1 tolerance column of BS.468-4 (dB), rendered from the PDF. The
        // 31.5 kHz row (+2.8 / −∞) is excluded: it is outside the mask.
        const COLUMN: [(f64, f64); 20] = [
            (31.5, 2.0),
            (63.0, 1.4),
            (100.0, 1.0),
            (200.0, 0.85),
            (400.0, 0.7),
            (800.0, 0.55),
            (1000.0, 0.5),
            (2000.0, 0.5),
            (3150.0, 0.5),
            (4000.0, 0.5),
            (5000.0, 0.5),
            (6300.0, 0.0),
            (7100.0, 0.2),
            (8000.0, 0.4),
            (9000.0, 0.6),
            (10000.0, 0.8),
            (12500.0, 1.2),
            (14000.0, 1.4),
            (16000.0, 1.6),
            (20000.0, 2.0),
        ];
        // The standard states the entries to 0.1 dB, so allow the rounding.
        for (f, want) in COLUMN {
            let got = itu468_tolerance_db(f).unwrap();
            assert!((got - want).abs() <= 0.05, "{f} Hz: {got} vs {want}");
        }
        assert!(itu468_tolerance_db(31_500.0).is_none());
    }

    #[test]
    fn itu468_response_extends_to_31_5_khz() {
        assert!(approx(itu468_weighting_db(31_500.0), -42.7, 1e-12));
        assert!(itu468_weighting_db(25_000.0) < -22.2);
        assert!(itu468_weighting_db(25_000.0) > -42.7);
    }

    #[test]
    fn digital_a_weighting_matches_curve() {
        let fs = 48000.0;
        let f = AWeightingFilter::design(fs).unwrap();
        for &freq in &[63.0, 250.0, 1000.0, 4000.0, 8000.0] {
            let digital = f.magnitude_db(freq, fs);
            let curve = a_weighting_db(freq);
            // Plain bilinear transform (no prewarp): ~0.1 dB typical error,
            // ~0.6 dB worst case at the HF corner.
            assert!(
                (digital - curve).abs() < 0.75,
                "digital {digital:.2} dB vs curve {curve:.2} dB at {freq} Hz"
            );
        }
    }

    #[test]
    fn filter_rejects_invalid_rates() {
        assert!(AWeightingFilter::design(0.0).is_err());
        assert!(AWeightingFilter::design(8000.0).is_err());
        assert!(AWeightingFilter::design(48000.0).is_ok());
    }
}
