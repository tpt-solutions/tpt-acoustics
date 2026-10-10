//! Octave and third-octave frequency band tables.
//!
//! Center frequencies follow ANSI S1.11 (base-10 system, reference 1 kHz):
//! the exact third-octave centers are `1000·10^{n/10}` Hz for integer `n`,
//! with the nominal (rounded) values used for reporting in the tables below.

/// Nominal third-octave band center frequencies in Hz, 20 Hz – 20 kHz
/// (ANSI/IEC nominal values).
pub const THIRD_OCTAVE_NOMINAL_HZ: [f64; 31] = [
    20.0, 25.0, 31.5, 40.0, 50.0, 63.0, 80.0, 100.0, 125.0, 160.0, 200.0, 250.0, 315.0, 400.0,
    500.0, 630.0, 800.0, 1000.0, 1250.0, 1600.0, 2000.0, 2500.0, 3150.0, 4000.0, 5000.0, 6300.0,
    8000.0, 10_000.0, 12_500.0, 16_000.0, 20_000.0,
];

/// Nominal octave band center frequencies in Hz, 31.5 Hz – 16 kHz.
pub const OCTAVE_NOMINAL_HZ: [f64; 10] = [
    31.5, 63.0, 125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0, 8000.0, 16_000.0,
];

/// Exact base-10 third-octave center frequency for band index `n`
/// (`n = 0` at 1 kHz): `1000·10^{n/10}` Hz.
#[must_use]
pub fn third_octave_exact_hz(n: i32) -> f64 {
    1000.0 * libm::pow(10.0, n as f64 / 10.0)
}

/// Exact base-10 octave center frequency for band index `n`
/// (`n = 0` at 1 kHz): `1000·10^{3n/10}` Hz.
#[must_use]
pub fn octave_exact_hz(n: i32) -> f64 {
    1000.0 * libm::pow(10.0, 3.0 * n as f64 / 10.0)
}

/// Exact band-edge frequencies of a base-10 third-octave band: the pair
/// `(f_c·10^{−1/20}, f_c·10^{+1/20})` — bandwidth ratio `10^{0.1} ≈ 2^{1/3}`.
#[must_use]
pub fn third_octave_edges_hz(center_hz: f64) -> (f64, f64) {
    let half = libm::pow(10.0, 0.05);
    (center_hz / half, center_hz * half)
}

/// Exact band-edge frequencies of a base-10 octave band: the pair
/// `(f_c·10^{−3/20}, f_c·10^{+3/20})` — bandwidth ratio `10^{0.3} ≈ 2`.
#[must_use]
pub fn octave_edges_hz(center_hz: f64) -> (f64, f64) {
    let half = libm::pow(10.0, 0.15);
    (center_hz / half, center_hz * half)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn tables_are_monotonic_and_contain_1khz() {
        for table in [
            THIRD_OCTAVE_NOMINAL_HZ.as_slice(),
            OCTAVE_NOMINAL_HZ.as_slice(),
        ] {
            for w in table.windows(2) {
                assert!(w[0] < w[1], "non-monotonic at {}", w[0]);
            }
            assert!(table.contains(&1000.0));
        }
        assert_eq!(THIRD_OCTAVE_NOMINAL_HZ.len(), 31);
    }

    #[test]
    fn exact_centers_match_formula() {
        assert!(approx(third_octave_exact_hz(0), 1000.0, 1e-9));
        // n = -17 -> 19.95 Hz (nominal 20 Hz) in the base-10 system.
        assert!(approx(third_octave_exact_hz(-17), 19.95, 0.01));
        // Octave: n = 1 -> 1.995 kHz (nominal 2 kHz).
        assert!(approx(octave_exact_hz(1), 1995.26, 0.01));
    }

    #[test]
    fn third_octave_bandwidth_ratio() {
        let (lo, hi) = third_octave_edges_hz(1000.0);
        assert!(approx(hi / lo, libm::pow(10.0, 0.1), 1e-12));
        // 1 kHz third-octave band: ~891 – 1122 Hz.
        assert!(approx(lo, 890.9, 0.5) && approx(hi, 1122.5, 0.5));
    }

    #[test]
    fn octave_bandwidth_ratio() {
        let (lo, hi) = octave_edges_hz(1000.0);
        assert!(approx(hi / lo, libm::pow(10.0, 0.3), 1e-12));
        assert!(approx(lo, 707.9, 0.1) && approx(hi, 1412.5, 0.1));
    }
}
