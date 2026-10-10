//! Lighthill's acoustic analogy.
//!
//! Lighthill (1952) rewrites the exact equations of motion as a wave equation
//! in a fictitious acoustic medium with the quadrupole source
//!
//! ```text
//! T_ij = rho u_i u_j + (p - c0^2 rho') delta_ij - tau_ij
//! ```
//!
//! (Reynolds-stress + entropy + viscous contributions; for high-Reynolds
//! cold flow the Reynolds-stress term `rho u_i u_j` dominates). In the far
//! field of a compact source region the density fluctuation at observer `x`
//! is the double-divergence of the retarded source distribution, which for a
//! discretized volume reduces to the standard quadrupole sum implemented
//! here (`1/(4 pi c0^2) * d^2/dt^2 integral T_rr / a` in the compact limit,
//! with `T_rr` the source strength projected on the radiation direction).

use crate::source::QuadrupoleCell;

/// Acoustic environment (density and speed of sound of the quiescent medium).
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

/// Far-field acoustic density `rho'(x, t)` (kg/m³) from a discretized
/// quadrupole volume: each cell contributes its compact-source term
///
/// `rho'_cell = volume / (4 pi c0^2 a) * d^2 T_rr / dt^2 (t - a/c0)`
///
/// evaluated with the retarded time and a linear `d²/dt²` stencil.
///
/// # Panics
/// Panics if `cells` is empty (the acoustic density of no source is
/// undefined at this level; callers should treat that as zero).
#[must_use]
pub fn far_field_density(
    ambient: &Ambient,
    cells: &[QuadrupoleCell],
    observer: tpt_acou_core::Vec3,
    time: f64,
) -> f64 {
    let mut acc = 0.0;
    for cell in cells {
        let d = observer.distance(cell.centroid);
        if d < 1e-9 {
            continue;
        }
        let t_ret = time - d / ambient.c0;
        let ddtt = second_derivative_at(cell, t_ret);
        acc += cell.volume * ddtt / (4.0 * core::f64::consts::PI * ambient.c0 * ambient.c0 * d);
    }
    acc
}

/// `d²/dt²` of a uniformly sampled signal at time `t`, via the centered
/// finite-difference on the interpolated array (returns 0 outside the
/// sampled window).
fn second_derivative_at(cell: &QuadrupoleCell, t: f64) -> f64 {
    let dt = 1.0 / cell.t_rr.sample_rate;
    let s = &cell.t_rr.samples;
    if s.len() < 3 {
        return 0.0;
    }
    let h = 2.0 * dt;
    let f_plus = cell.t_rr.at(t + dt);
    let f_mid = cell.t_rr.at(t);
    let f_minus = cell.t_rr.at(t - dt);
    let _ = h;
    (f_plus - 2.0 * f_mid + f_minus) / (dt * dt)
}

/// Far-field mean-square pressure from the Lighthill integral in the
/// acoustic analogy's compact limit, given the *time series* of the
/// far-field density at the observer: `p_rms = c0^2 * rho'_rms`.
#[must_use]
pub fn rms_pressure(density_series: &[f64], c0: f64) -> f64 {
    let n = density_series.len();
    if n == 0 {
        return 0.0;
    }
    let mean = density_series.iter().sum::<f64>() / n as f64;
    let ms = density_series
        .iter()
        .map(|d| (d - mean) * (d - mean))
        .sum::<f64>()
        / n as f64;
    c0 * c0 * libm::sqrt(ms)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source::TimeSignal;
    use alloc::vec::Vec;

    #[test]
    fn quadrupole_amplitude_and_retarded_time() {
        let amb = Ambient::air_20c();
        // One cell of volume 1 m^3 with T_rr = sin(2*pi*f*t), f = 100 Hz.
        let f = 100.0_f64;
        let n = 4096;
        let sr = 8192.0_f64;
        let samples: Vec<f64> = (0..n)
            .map(|i| libm::sin(core::f64::consts::TAU * f * i as f64 / sr))
            .collect();
        let cell = QuadrupoleCell {
            centroid: tpt_acou_core::Vec3::new(0.0, 0.0, 0.0),
            volume: 1.0,
            t_rr: TimeSignal::new(samples, sr).unwrap(),
        };
        let observer = tpt_acou_core::Vec3::new(0.0, 0.0, 343.21); // 1 s travel
                                                                   // Analytic: d2/dt2 sin(w t) = -w^2 sin(w t) -> rho' =
                                                                   // 1/(4 pi c0^2 a) * (-w^2) sin(w(t - a/c0)); peak value:
        let w = core::f64::consts::TAU * f;
        let peak = w * w / (4.0 * core::f64::consts::PI * amb.c0 * amb.c0 * 343.21);
        let mut max_got = 0.0_f64;
        for k in 0..400 {
            let t = 1.0 + k as f64 / sr; // after the wavefront arrives
            let got = far_field_density(&amb, &alloc::vec![cell.clone()], observer, t).abs();
            max_got = max_got.max(got);
        }
        assert!(
            (max_got - peak).abs() / peak < 0.05,
            "peak {max_got:.3e} vs analytic {peak:.3e}"
        );
    }

    #[test]
    fn rms_pressure_from_density_series() {
        let s: Vec<f64> = (0..1000).map(|i| libm::sin(i as f64)).collect();
        let rms = rms_pressure(&s, 343.21);
        let want = 343.21 * 343.21 * libm::sqrt(0.5); // p_rms = c0^2 * rho'_rms
        assert!((rms - want).abs() / want < 0.05);
    }
}
