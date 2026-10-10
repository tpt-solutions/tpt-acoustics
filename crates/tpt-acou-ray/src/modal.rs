//! Modal transfer function of a rectangular room with rigid walls.
//!
//! The Helmholtz Green's function `G` with `∇²G + k²G = −δ(r − s)` has the
//! modal expansion
//!
//! ```text
//! G(r, s; ω) = c² Σ_n ψ_n(r) ψ_n(s) / ( N_n (ω_n² − ω² − 2 j δ ω) )
//! ```
//!
//! with Neumann eigenfunctions `ψ_n = ∏ cos(n_i π x_i / L_i)`, norms
//! `N_n = ∏ N_i` (`N_i = L_i/2` for `n_i > 0`, `L_i` for `n_i = 0`), mode
//! frequencies `ω_n = c π √Σ(n_i/L_i)²`, and a common damping rate
//! `δ = 3 ln10 / T` from the Sabine time `T`. This is the same free-field
//! normalisation as the image-source transfer function, so the two can be
//! combined in one spectrum.
//!
//! Below the Schroeder frequency the response is modal, and the image-source
//! sum (which is accurate at high frequency) does not resolve individual
//! modes. [`crate::rir::hybrid_broadband_rir`] uses the modal sum below the
//! Schroeder frequency and the image-source sum above it.

use alloc::vec::Vec;

use tpt_acou_core::error::AcousticError;
use tpt_math_linalg_complex::Complex;

/// Modal sum for a rectangular rigid room.
#[derive(Debug, Clone, PartialEq)]
pub struct RigidBoxModes {
    dims: [f64; 3],
    c: f64,
    damping: f64,
    /// `(ω_n, ψ_n(r)ψ_n(s)/N_n)` per retained mode.
    modes: Vec<(f64, f64)>,
}

impl RigidBoxModes {
    /// Builds the mode set for source `s` and receiver `r`, keeping every mode
    /// with frequency at or below `f_max_hz`.
    ///
    /// `t60_s` sets the common damping `δ = 3 ln10 / T`. Use a large value
    /// (or `f64::INFINITY`) for a lossless room.
    ///
    /// # Errors
    /// Returns [`AcousticError::InvalidArgument`] for non-positive dimensions,
    /// speed, `f_max_hz` or `t60_s`, or for positions outside the room.
    pub fn new(
        dims: [f64; 3],
        s: [f64; 3],
        r: [f64; 3],
        c: f64,
        f_max_hz: f64,
        t60_s: f64,
    ) -> Result<Self, AcousticError> {
        if dims.iter().any(|d| !(d.is_finite() && *d > 0.0)) {
            return Err(AcousticError::InvalidArgument(
                "room dimensions must be positive",
            ));
        }
        if !(c.is_finite() && c > 0.0) || !(f_max_hz.is_finite() && f_max_hz > 0.0) {
            return Err(AcousticError::InvalidArgument(
                "speed and f_max must be positive",
            ));
        }
        if t60_s.is_nan() || t60_s <= 0.0 {
            return Err(AcousticError::InvalidArgument("T60 must be positive"));
        }
        for p in [s, r] {
            for i in 0..3 {
                if !(p[i] > 0.0 && p[i] < dims[i]) {
                    return Err(AcousticError::InvalidArgument(
                        "source and receiver must lie strictly inside the room",
                    ));
                }
            }
        }
        let damping = if t60_s.is_finite() {
            3.0 * core::f64::consts::LN_10 / t60_s
        } else {
            0.0
        };
        let omega_max = core::f64::consts::TAU * f_max_hz;
        // Largest mode index per axis that can still lie below omega_max.
        let nmax: [i64; 3] = [0, 1, 2]
            .map(|i| libm::floor(omega_max * dims[i] / (core::f64::consts::PI * c)) as i64 + 1);
        let mut modes = Vec::new();
        for nx in 0..=nmax[0] {
            for ny in 0..=nmax[1] {
                for nz in 0..=nmax[2] {
                    if nx == 0 && ny == 0 && nz == 0 {
                        continue;
                    }
                    let n = [nx, ny, nz];
                    let kx = nx as f64 * core::f64::consts::PI / dims[0];
                    let ky = ny as f64 * core::f64::consts::PI / dims[1];
                    let kz = nz as f64 * core::f64::consts::PI / dims[2];
                    let omega_n = c * libm::sqrt(kx * kx + ky * ky + kz * kz);
                    if omega_n > omega_max {
                        continue;
                    }
                    let mut psi = 1.0;
                    let mut norm = 1.0;
                    for i in 0..3 {
                        let arg_s = n[i] as f64 * core::f64::consts::PI * s[i] / dims[i];
                        let arg_r = n[i] as f64 * core::f64::consts::PI * r[i] / dims[i];
                        psi *= libm::cos(arg_s) * libm::cos(arg_r);
                        norm *= if n[i] == 0 { dims[i] } else { dims[i] / 2.0 };
                    }
                    modes.push((omega_n, psi / norm));
                }
            }
        }
        Ok(Self {
            dims,
            c,
            damping,
            modes,
        })
    }

    /// Room dimensions.
    #[must_use]
    pub fn dims(&self) -> [f64; 3] {
        self.dims
    }

    /// Number of retained modes (the `(0,0,0)` mode is excluded).
    #[must_use]
    pub fn len(&self) -> usize {
        self.modes.len()
    }

    /// Whether no modes were retained.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.modes.is_empty()
    }

    /// Green's function value at angular frequency `ω` (rad/s).
    #[must_use]
    pub fn green_at_omega(&self, omega: f64) -> Complex {
        let mut sum = Complex::new(0.0, 0.0);
        for &(omega_n, weight) in &self.modes {
            let denom = Complex::new(
                omega_n * omega_n - omega * omega,
                -2.0 * self.damping * omega,
            );
            sum = sum + Complex::new(weight, 0.0) * denom.inv().unwrap_or(Complex::new(0.0, 0.0));
        }
        sum * (self.c * self.c)
    }

    /// Green's function at frequency `f` in hertz.
    #[must_use]
    pub fn green_at_hz(&self, f: f64) -> Complex {
        self.green_at_omega(core::f64::consts::TAU * f)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const DIMS: [f64; 3] = [5.0, 4.0, 3.0];
    const C: f64 = 343.0;

    #[test]
    fn reciprocity_is_exact() {
        let s = [1.2, 1.1, 1.0];
        let r = [3.7, 2.6, 1.8];
        let a = RigidBoxModes::new(DIMS, s, r, C, 300.0, 1.0).unwrap();
        let b = RigidBoxModes::new(DIMS, r, s, C, 300.0, 1.0).unwrap();
        for f in [30.0, 80.0, 150.0] {
            let ga = a.green_at_hz(f);
            let gb = b.green_at_hz(f);
            assert!((ga - gb).norm() < 1e-12 * ga.norm().max(1e-30), "f={f}");
        }
    }

    #[test]
    fn first_axial_mode_sets_the_resonance() {
        // Lossless room: the response peaks at the first axial mode c/(2Lx).
        let s = [1.0, 1.0, 1.0];
        let r = [4.0, 3.0, 2.0];
        let m = RigidBoxModes::new(DIMS, s, r, C, 400.0, f64::INFINITY).unwrap();
        let f1 = C / (2.0 * DIMS[0]);
        let below = m.green_at_hz(f1 * 0.9).norm();
        let at = m.green_at_hz(f1 * 1.0005).norm();
        assert!(at > below, "resonance not at {f1} Hz");
    }

    #[test]
    fn damping_limits_the_peak() {
        let s = [1.0, 1.0, 1.0];
        let r = [4.0, 3.0, 2.0];
        let lossless = RigidBoxModes::new(DIMS, s, r, C, 400.0, f64::INFINITY).unwrap();
        let damped = RigidBoxModes::new(DIMS, s, r, C, 400.0, 0.3).unwrap();
        let f1 = C / (2.0 * DIMS[0]);
        assert!(damped.green_at_hz(f1).norm() < lossless.green_at_hz(f1 * 1.0005).norm());
    }

    #[test]
    fn invalid_geometry_is_rejected() {
        assert!(RigidBoxModes::new(DIMS, [0.0, 1.0, 1.0], [1.0, 1.0, 1.0], C, 100.0, 1.0).is_err());
        assert!(RigidBoxModes::new(DIMS, [1.0; 3], [1.0; 3], -1.0, 100.0, 1.0).is_err());
        assert!(RigidBoxModes::new(DIMS, [1.0; 3], [1.0; 3], C, 100.0, 0.0).is_err());
    }
}
