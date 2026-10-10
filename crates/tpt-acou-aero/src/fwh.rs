//! Ffowcs Williams–Hawkings far-field forms for moving sources.
//!
//! The FW-H equation embeds a body surface in the fluid; its solution splits
//! into a thickness (monopole, `F1`-type) term from the surface normal
//! velocity and a loading (dipole, `F2`-type) term from the surface pressure.
//! For *compact* moving sources the far field reduces to
//!
//! ```text
//! p'_T(x, t) = 1/(4 pi)   d/dt [ Q(t_r) / (a (1 - M_r)) ]        (monopole)
//! p'_L(x, t) = 1/(4 pi c0) d/dt [ F(t_r) . x_hat / (a (1 - M_r)) ] (dipole)
//! ```
//!
//! with the source-emission time `t_r` solving the retarded-time equation
//! `t - t_r = |x - y(t_r)|/c0`, source strength `Q` (volume displacement
//! rate), force `F`, and `M_r` the source Mach component toward the observer.
//! Time differentiation is carried through the implicit `t_r(t)` dependence
//! with the standard factor `1/(1 - M_r)` — the convective amplification
//! that makes approaching sources louder.
//!
//! Rotating machinery: `rotating_sources` evaluates a ring of blades with
//! per-blade monopole/dipole strength histories, producing the
//! blade-passing-frequency (BPF = B·Ω/2π) harmonic ladder.

use alloc::vec;
use alloc::vec::Vec;

use tpt_acou_core::Vec3;

use crate::source::CompactSource;

/// Acoustic environment.
#[derive(Debug, Clone, Copy)]
pub struct Ambient {
    /// Mean density (kg/m³).
    pub density: f64,
    /// Mean speed of sound c0 (m/s).
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

/// Solves the source-retarded time `t_r` for a source moving at constant
/// velocity: `t - t_r = |x - y(t_r)|/c0` via fixed-point iteration (converges
/// quadratically for subsonic source speeds).
#[must_use]
pub fn retarded_time(x: Vec3, y0: Vec3, source_velocity: Vec3, c0: f64, t: f64) -> f64 {
    let mut tr = t;
    for _ in 0..32 {
        let d = x - (y0 + source_velocity * tr);
        let a = d.norm();
        let tr_new = t - a / c0;
        if (tr_new - tr).abs() < 1e-12 {
            return tr_new;
        }
        tr = tr_new;
    }
    tr
}

/// Far-field pressure (Pa) of a single moving monopole (thickness) source
/// with volume rate history `Q(t)` (m³/s, `TimeSignal` sampled in the
/// source's own time frame).
///
/// `p'_T = 1/(4 pi) * d/dt[ Q(t_r) / (a (1 - M_r)) ]`, differentiated
/// numerically in observer time with a central difference (step = one
/// sample).
///
/// Returns `NaN` where the far field is undefined: observer at the source
/// position, or `1 - M_r <= 0` (no subsonic retarded time).
#[must_use]
pub fn monopole_pressure(
    ambient: &Ambient,
    source: &CompactSource,
    source_velocity: Vec3,
    observer: Vec3,
    time: f64,
) -> f64 {
    let h = 1.0 / source.strength.sample_rate;
    let f = |t: f64| {
        let tr = retarded_time(observer, source.position, source_velocity, ambient.c0, t);
        let d = observer - (source.position + source_velocity * tr);
        let a = d.norm();
        if a <= 0.0 {
            return f64::NAN;
        }
        let mr = source_velocity.dot(d / a) / ambient.c0;
        if 1.0 - mr <= 0.0 {
            return f64::NAN;
        }
        source.strength.at(tr) / (a * (1.0 - mr))
    };
    (f(time + h) - f(time - h)) / (4.0 * core::f64::consts::PI * 2.0 * h)
}

/// Far-field pressure (Pa) of a single moving dipole (loading) source with
/// force history `F(t)` (N) along `source.direction`.
///
/// `p'_L = 1/(4 pi c0) * d/dt[ F(t_r) . x_hat / (a (1 - M_r)) ]`, central
/// difference in observer time. Returns `NaN` where the far field is
/// undefined (see [`monopole_pressure`]).
#[must_use]
pub fn dipole_pressure(
    ambient: &Ambient,
    source: &CompactSource,
    source_velocity: Vec3,
    observer: Vec3,
    time: f64,
) -> f64 {
    let h = 1.0 / source.strength.sample_rate;
    let f = |t: f64| {
        let tr = retarded_time(observer, source.position, source_velocity, ambient.c0, t);
        let d = observer - (source.position + source_velocity * tr);
        let a = d.norm();
        if a <= 0.0 {
            return f64::NAN;
        }
        let mr = source_velocity.dot(d / a) / ambient.c0;
        if 1.0 - mr <= 0.0 {
            return f64::NAN;
        }
        let x_hat = d / a;
        source.direction.dot(x_hat) * source.strength.at(tr) / (a * (1.0 - mr))
    };
    (f(time + h) - f(time - h)) / (4.0 * core::f64::consts::PI * ambient.c0 * 2.0 * h)
}

/// Builds the `B` compact monopole sources of a ring of radius `radius` in
/// its t = 0 layout (blades evenly spaced about the z axis), each with
/// constant strength history `q` (m³/s) sampled at `sample_rate` for
/// `n_samples`.
///
/// This is a static snapshot: it carries no rotation. The rotating dynamics,
/// including the emission-time blade positions and velocities, live in
/// [`rotating_ring_pressure`]. Use this builder when the sources are needed
/// as explicit [`CompactSource`] values, e.g. for the moving-monopole
/// functions.
///
/// # Errors
/// Propagates [`crate::source::TimeSignal::new`] errors.
pub fn rotating_monopole_ring(
    blade_count: u32,
    radius: f64,
    q: f64,
    sample_rate: f64,
    n_samples: usize,
) -> Result<Vec<CompactSource>, tpt_acou_core::AcousticError> {
    let mut sources = Vec::new();
    for b in 0..blade_count {
        let theta = core::f64::consts::TAU * b as f64 / f64::from(blade_count);
        sources.push(CompactSource {
            position: Vec3::new(radius * libm::cos(theta), radius * libm::sin(theta), 0.0),
            direction: Vec3::ZERO,
            strength: crate::source::TimeSignal::new(vec![q; n_samples], sample_rate)?,
        });
    }
    Ok(sources)
}

/// Total far-field thickness pressure (Pa) at `observer` and observer time
/// `time` of `blade_count` blades, each a monopole of constant strength `q`
/// (m³/s) on a ring of radius `radius` rotating at `omega` (rad/s) about the
/// z axis.
///
/// Each blade is a moving monopole evaluated at its emission time: position
/// `R (cos θ, sin θ, 0)` and velocity `R ω (-sin θ, cos θ, 0)` with
/// `θ = θ_b + ω t_r`. The far-field form `1/(4π) d/dt[ q / (a (1 - M_r)) ]`
/// includes the convective factor `1/(1 - M_r)`, with `M_r` the emission-time
/// Mach number toward the observer.
///
/// Returns `0.0` for `omega = 0` (a fixed ring of constant strength radiates
/// nothing), and `NaN` where the far field is undefined (observer at a blade,
/// or `1 - M_r <= 0`).
#[must_use]
pub fn rotating_ring_pressure(
    ambient: &Ambient,
    blade_count: u32,
    radius: f64,
    omega: f64,
    q: f64,
    observer: Vec3,
    time: f64,
) -> f64 {
    if omega == 0.0 {
        return 0.0;
    }
    // Central-difference step: 1/32 of a revolution.
    let h = core::f64::consts::TAU / (32.0 * omega.abs());
    let mut total = 0.0;
    for b in 0..blade_count {
        let theta0 = core::f64::consts::TAU * b as f64 / f64::from(blade_count);
        let f = |t: f64| {
            let tr = ring_retarded_time(observer, radius, omega, theta0, ambient.c0, t);
            let angle = theta0 + omega * tr;
            let (s, c) = (libm::sin(angle), libm::cos(angle));
            let y = Vec3::new(radius * c, radius * s, 0.0);
            let vel = Vec3::new(-radius * omega * s, radius * omega * c, 0.0);
            let d = observer - y;
            let a = d.norm();
            if a <= 0.0 {
                return f64::NAN;
            }
            let mr = vel.dot(d / a) / ambient.c0;
            if 1.0 - mr <= 0.0 {
                return f64::NAN;
            }
            q / (a * (1.0 - mr))
        };
        let p = (f(time + h) - f(time - h)) / (2.0 * h);
        total += p / (4.0 * core::f64::consts::PI);
    }
    total
}

/// Fixed-point retarded time for a blade on a rotating ring.
fn ring_retarded_time(
    observer: Vec3,
    radius: f64,
    omega: f64,
    theta0: f64,
    c0: f64,
    t: f64,
) -> f64 {
    let mut tr = t;
    for _ in 0..40 {
        let angle = theta0 + omega * tr;
        let y = Vec3::new(radius * libm::cos(angle), radius * libm::sin(angle), 0.0);
        let tr_new = t - observer.distance(y) / c0;
        if (tr_new - tr).abs() < 1e-12 {
            return tr_new;
        }
        tr = tr_new;
    }
    tr
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::source::TimeSignal;

    #[test]
    fn stationary_monopole_matches_curle_like_scaling() {
        let amb = Ambient::air_20c();
        // Q(t) = Q0 sin(w t): p = 1/(4 pi a) dQ/dt(t - a/c0).
        let q0 = 1.0;
        let f = 100.0;
        let w = core::f64::consts::TAU * f;
        let n = 8192;
        let sr = 16384.0_f64;
        let samples: Vec<f64> = (0..n).map(|i| q0 * libm::sin(w * i as f64 / sr)).collect();
        let source = CompactSource {
            position: Vec3::ZERO,
            direction: Vec3::ZERO,
            strength: TimeSignal::new(samples, sr).unwrap(),
        };
        let obs = Vec3::new(0.0, 0.0, 343.21);
        let peak = q0 * w / (4.0 * core::f64::consts::PI * 343.21);
        let mut max_got = 0.0_f64;
        for k in 0..800 {
            let t = 1.0 + k as f64 / sr;
            max_got = max_got.max(monopole_pressure(&amb, &source, Vec3::ZERO, obs, t).abs());
        }
        assert!(
            (max_got - peak).abs() / peak < 0.05,
            "peak {max_got:.3e} vs {peak:.3e}"
        );
    }

    #[test]
    fn moving_monopole_is_amplified_upstream() {
        // A source moving TOWARD the observer radiates a louder signal than
        // the stationary equivalent (convective amplification 1/(1-Mr)^2
        // order effects); away, quieter.
        let amb = Ambient::air_20c();
        let n = 8192;
        let sr = 16384.0_f64;
        let samples: Vec<f64> = (0..n)
            .map(|i| libm::sin(core::f64::consts::TAU * 40.0 * i as f64 / sr))
            .collect();
        let source = CompactSource {
            position: Vec3::new(-100.0, 0.0, 0.0),
            direction: Vec3::ZERO,
            strength: TimeSignal::new(samples, sr).unwrap(),
        };
        let obs = Vec3::new(100.0, 0.0, 0.0);
        let vel = Vec3::new(50.0, 0.0, 0.0); // Mach 0.15 toward +x
        let rms = |vel| {
            let mut acc = 0.0;
            let mut cnt = 0_usize;
            for k in 0..4000 {
                let t = 0.4 + k as f64 / sr;
                let p = monopole_pressure(&amb, &source, vel, obs, t);
                acc += p * p;
                cnt += 1;
            }
            libm::sqrt(acc / cnt as f64)
        };
        let approaching = rms(vel);
        let receding = rms(Vec3::new(-50.0, 0.0, 0.0));
        assert!(
            approaching > receding * 1.1,
            "approaching {approaching:.3e} vs receding {receding:.3e}"
        );
    }

    #[test]
    fn rotating_ring_has_bpf_tone() {
        let amb = Ambient::air_20c();
        let b = 4_u32;
        let omega = core::f64::consts::TAU * 50.0; // 50 rev/s
        let bpf = f64::from(b) * omega / core::f64::consts::TAU; // 200 Hz
        let radius = 0.25;
        let obs = Vec3::new(2.0, 0.0, 2.0);
        let sr = 8192.0_f64;
        let n = 8192;
        let mut series = Vec::with_capacity(n);
        for i in 0..n {
            let t = 0.5 + i as f64 / sr;
            series.push(rotating_ring_pressure(&amb, b, radius, omega, 1.0, obs, t));
        }
        // Goertzel amplitude at the BPF vs at half the BPF (should dominate).
        let goertzel = |freq: f64| -> f64 {
            let wk = core::f64::consts::TAU * freq / sr;
            let coeff = 2.0 * libm::cos(wk);
            let mut s1 = 0.0_f64;
            let mut s2 = 0.0_f64;
            for &x in &series {
                let s0 = x + coeff * s1 - s2;
                s2 = s1;
                s1 = s0;
            }
            libm::sqrt(s1 * s1 + s2 * s2 - coeff * s1 * s2)
        };
        let amp_bpf = goertzel(bpf);
        let amp_off = goertzel(bpf / 2.0);
        assert!(
            amp_bpf > 4.0 * amp_off.max(1e-12),
            "BPF amp {amp_bpf:.3e} vs off-tone {amp_off:.3e}"
        );
    }

    #[test]
    fn stationary_ring_radiates_nothing() {
        // omega = 0 used to divide by zero; a fixed constant-strength ring
        // has no time-varying geometry, so its far field is exactly zero.
        let amb = Ambient::air_20c();
        let p = rotating_ring_pressure(&amb, 4, 0.25, 0.0, 1.0, Vec3::new(2.0, 0.0, 2.0), 0.5);
        assert_eq!(p, 0.0);
    }

    #[test]
    fn undefined_far_field_is_nan_not_inf() {
        // Observer at the source: a = 0 has no far-field value.
        let amb = Ambient::air_20c();
        let source = CompactSource {
            position: Vec3::ZERO,
            direction: Vec3::ZERO,
            strength: TimeSignal::new(vec![1.0; 64], 1000.0).unwrap(),
        };
        let p = monopole_pressure(&amb, &source, Vec3::ZERO, Vec3::ZERO, 0.01);
        assert!(p.is_nan(), "got {p}");
    }

    #[test]
    fn ring_pressure_is_directional_and_finite() {
        // Weak check: with the convective factor the far field depends on the
        // observer side, and stays finite. An analytic benchmark for the
        // amplitude is still open (see todo.md 4.2).
        let amb = Ambient::air_20c();
        let omega = core::f64::consts::TAU * 200.0; // Tip Mach ~0.12 at R = 0.25
        let rms = |obs: Vec3| {
            let sr = 16384.0_f64;
            let mut acc = 0.0;
            for k in 0..4096 {
                let t = 0.5 + k as f64 / sr;
                let p = rotating_ring_pressure(&amb, 2, 0.25, omega, 1.0, obs, t);
                acc += p * p;
            }
            libm::sqrt(acc / 4096.0)
        };
        let toward = rms(Vec3::new(0.0, 5.0, 0.0));
        let away = rms(Vec3::new(0.0, -5.0, 0.0));
        assert!(toward.is_finite() && away.is_finite());
        assert!(toward != away, "no directivity: {toward:e} vs {away:e}");
    }

    #[test]
    fn dipole_directivity_of_rotating_loading() {
        // A dipole source: observer along the force axis hears more than an
        // observer at 90 degrees (stationary case).
        let amb = Ambient::air_20c();
        let n = 8192;
        let sr = 16384.0_f64;
        let samples: Vec<f64> = (0..n)
            .map(|i| libm::sin(core::f64::consts::TAU * 60.0 * i as f64 / sr))
            .collect();
        let source = CompactSource {
            position: Vec3::ZERO,
            direction: Vec3::new(0.0, 0.0, 1.0),
            strength: TimeSignal::new(samples, sr).unwrap(),
        };
        let on_axis = Vec3::new(0.0, 0.0, 200.0);
        let off_axis = Vec3::new(200.0, 0.0, 0.0);
        let rms = |obs: Vec3| {
            let mut acc = 0.0;
            for k in 0..4000 {
                let t = 0.5 + k as f64 / sr;
                let p = dipole_pressure(&amb, &source, Vec3::ZERO, obs, t);
                acc += p * p;
            }
            libm::sqrt(acc / 4000.0)
        };
        let on = rms(on_axis);
        let off = rms(off_axis);
        assert!(off < on * 0.15, "on {on:.3e} off {off:.3e}");
    }
}
