//! Deterministic pseudo-random number generation (no external `rand` dep).
//!
//! xoshiro256** seeded through SplitMix64: small, fast, statistically sound
//! enough for Monte Carlo room acoustics, fully deterministic across platforms
//! for reproducible auralization and property-based tests.

/// xoshiro256** generator state (4 × u64).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Xoshiro256 {
    s: [u64; 4],
}

impl Xoshiro256 {
    /// Seeds from a single u64 via SplitMix64 (avoids the all-zero state).
    #[must_use]
    pub fn new(seed: u64) -> Self {
        let mut sm = seed;
        let mut next = || {
            sm = sm.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut z = sm;
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        };
        Self {
            s: [next(), next(), next(), next()],
        }
    }

    /// Next raw u64 (xoshiro256**).
    #[must_use]
    pub fn next_u64(&mut self) -> u64 {
        let result = self.s[1].wrapping_mul(5).rotate_left(7).wrapping_mul(9);
        let t = self.s[1] << 17;
        self.s[2] ^= self.s[0];
        self.s[3] ^= self.s[1];
        self.s[1] ^= self.s[2];
        self.s[0] ^= self.s[3];
        self.s[2] ^= t;
        self.s[3] = self.s[3].rotate_left(45);
        result
    }

    /// Uniform f64 in `[0, 1)`.
    #[must_use]
    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / 9_007_199_254_740_992.0)
    }

    /// Uniform f64 in `[a, b)`.
    #[must_use]
    pub fn next_range(&mut self, a: f64, b: f64) -> f64 {
        a + (b - a) * self.next_f64()
    }

    /// Standard normal N(0, 1) via Box–Muller (cos/sin branch).
    #[must_use]
    pub fn next_normal(&mut self) -> f64 {
        let u1 = self.next_f64().max(f64::MIN_POSITIVE);
        let u2 = self.next_f64();
        let r = libm::sqrt(-2.0 * libm::log(u1));
        r * libm::cos(core::f64::consts::TAU * u2)
    }

    /// Uniform direction on the unit sphere (Marsaglia ziggurat-free method:
    /// `z = 2u−1`, `phi = 2πv` — exact up to floating point).
    #[must_use]
    pub fn next_unit_direction(&mut self) -> crate::vec3::Vec3 {
        let z = 2.0 * self.next_f64() - 1.0;
        let phi = core::f64::consts::TAU * self.next_f64();
        let r = libm::sqrt((1.0 - z * z).max(0.0));
        crate::vec3::Vec3::new(r * libm::cos(phi), r * libm::sin(phi), z)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deterministic_for_same_seed() {
        let mut a = Xoshiro256::new(42);
        let mut b = Xoshiro256::new(42);
        for _ in 0..100 {
            assert_eq!(a.next_u64(), b.next_u64());
        }
    }

    #[test]
    fn uniforms_in_range_with_reasonable_mean() {
        let mut rng = Xoshiro256::new(7);
        let n = 20_000;
        let mut sum = 0.0;
        for _ in 0..n {
            let u = rng.next_f64();
            assert!((0.0..1.0).contains(&u));
            sum += u;
        }
        assert!((sum / n as f64 - 0.5).abs() < 0.02);
    }

    #[test]
    fn normal_moments_are_sane() {
        let mut rng = Xoshiro256::new(123);
        let n = 40_000;
        let (mut s1, mut s2) = (0.0, 0.0);
        for _ in 0..n {
            let x = rng.next_normal();
            s1 += x;
            s2 += x * x;
        }
        let mean = s1 / n as f64;
        let var = s2 / n as f64 - mean * mean;
        assert!(mean.abs() < 0.05, "mean {mean}");
        assert!((var - 1.0).abs() < 0.1, "var {var}");
    }

    #[test]
    fn unit_directions_are_unit() {
        let mut rng = Xoshiro256::new(9);
        for _ in 0..1000 {
            let d = rng.next_unit_direction();
            assert!((d.norm() - 1.0).abs() < 1e-12);
        }
    }
}
