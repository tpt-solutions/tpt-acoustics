//! Curle's extension of Lighthill's analogy for compact solid boundaries.
//!
//! Curle (1955): unsteady surface forces on a compact body radiate a dipole
//! whose far field is
//!
//! ```text
//! p'(x, t) = 1 / (4 pi c0 a) * x_hat_i * dF_i/dt (t - a/c0)
//! ```
//!
//! with `F_i(t)` the total unsteady force exerted by the body on the fluid
//! (N) and `a = |x|` the observer distance. The radiation scales as `U^6`
//! (via dF/dt ~ rho U^3 L for vortex shedding), the classic dipole law.

use crate::source::UnsteadyForce;

/// Acoustic environment.
#[derive(Debug, Clone, Copy)]
pub struct Ambient {
    /// Mean density `rho0` (kg/m³).
    pub density: f64,
    /// Mean speed of sound `c0` (m/s).
    pub c0: f64,
}

impl Ambient {
    /// Standard air at 20 °C.
    #[must_use]
    pub fn air_20c() -> Self {
        Self {
            density: 1.2041,
            c0: 343.21,
        }
    }
}

/// Far-field pressure `p'(x, t)` (Pa) of the compact dipole from an
/// unsteady force history.
///
/// The force derivative is evaluated with linear interpolation of the
/// sampled histories (forward difference).
///
/// # Panics
/// Panics if any force component's sample rate is non-positive (propagated
/// from [`crate::source::TimeSignal::at`] as zero output instead — the
/// sample-rate error surfaces at construction).
#[must_use]
pub fn far_field_pressure(
    ambient: &Ambient,
    force: &UnsteadyForce,
    observer: tpt_acou_core::Vec3,
    time: f64,
) -> f64 {
    let a = observer.norm();
    if a < 1e-9 {
        return 0.0;
    }
    let x_hat = observer / a;
    let t_ret = time - a / ambient.c0;
    let dfdt = |s: &crate::source::TimeSignal| {
        let dt = 1.0 / s.sample_rate;
        (s.at(t_ret + dt) - s.at(t_ret)) / dt
    };
    let f_dot = x_hat.x * dfdt(&force.fx) + x_hat.y * dfdt(&force.fy) + x_hat.z * dfdt(&force.fz);
    f_dot / (4.0 * core::f64::consts::PI * ambient.c0 * a)
}

/// Dipole directivity: the pressure amplitude at fixed distance scales with
/// `x_hat · f_hat` (cosine); returned as the fraction of the on-axis
/// amplitude for a force along `f_hat` observed from `x_hat`.
#[must_use]
pub fn dipole_directivity(f_hat: tpt_acou_core::Vec3, x_hat: tpt_acou_core::Vec3) -> f64 {
    f_hat.normalized().dot(x_hat.normalized())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source::TimeSignal;
    use alloc::vec;
    use alloc::vec::Vec;
    use tpt_acou_core::Vec3;

    #[test]
    fn dipole_amplitude_and_delay() {
        let amb = Ambient::air_20c();
        // Fz(t) = F0 sin(w t), observer on +z at 343.21 m (1 s delay).
        let f0 = 1.0;
        let f = 100.0;
        let w = core::f64::consts::TAU * f;
        let n = 8192;
        let sr = 16384.0_f64;
        let samples: Vec<f64> = (0..n).map(|i| f0 * libm::sin(w * i as f64 / sr)).collect();
        let force = UnsteadyForce {
            fx: TimeSignal::new(vec![0.0; n], sr).unwrap(),
            fy: TimeSignal::new(vec![0.0; n], sr).unwrap(),
            fz: TimeSignal::new(samples, sr).unwrap(),
            sample_rate: sr,
        };
        let obs = Vec3::new(0.0, 0.0, 343.21);
        // Analytic peak: (1/(4 pi c0 a)) * F0 * w.
        let peak = f0 * w / (4.0 * core::f64::consts::PI * amb.c0 * 343.21);
        let mut max_got = 0.0_f64;
        for k in 0..800 {
            let t = 1.0 + k as f64 / sr;
            max_got = max_got.max(far_field_pressure(&amb, &force, obs, t).abs());
        }
        assert!(
            (max_got - peak).abs() / peak < 0.05,
            "peak {max_got:.3e} vs {peak:.3e}"
        );
    }

    #[test]
    fn dipole_directivity_cosine() {
        let f_hat = Vec3::new(0.0, 0.0, 1.0);
        assert!((dipole_directivity(f_hat, Vec3::new(0.0, 0.0, 2.0)) - 1.0).abs() < 1e-12);
        assert!(dipole_directivity(f_hat, Vec3::new(3.0, 0.0, 0.0)).abs() < 1e-12);
        assert!(
            (dipole_directivity(f_hat, Vec3::new(1.0, 0.0, 1.0))
                - core::f64::consts::FRAC_1_SQRT_2)
                .abs()
                < 1e-12
        );
    }
}
