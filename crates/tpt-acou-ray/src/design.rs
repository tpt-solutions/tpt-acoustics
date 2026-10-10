//! Differentiable and inverse room design on the Sabine model.
//!
//! * [`sabine_with_gradient`] gives the Sabine reverberation time and its
//!   exact gradient with respect to the three room dimensions and the six
//!   per-wall absorption coefficients.
//! * [`absorption_for_rt60`] inverts the model for a uniform absorption.
//! * [`scale_absorption_to_rt60`] keeps a given wall-by-wall distribution and
//!   scales it to meet a target, which is the simplest inverse design that
//!   respects a chosen finish layout.
//! * [`rt60_uncertainty`] propagates per-wall absorption uncertainty by a
//!   seeded Monte Carlo.
//!
//! All of this is the Sabine model, valid for diffuse fields and moderate
//! absorption. Use the ray tracer for other cases.

use alloc::vec::Vec;

use tpt_acou_core::error::AcousticError;

use crate::rng::Xoshiro256;

/// Sabine constant in s/m.
const SABINE_CONSTANT: f64 = 0.161;

/// Validates a room size and per-wall absorption set.
fn check_room(dims: [f64; 3], alphas: &[f64; 6]) -> Result<(), AcousticError> {
    if dims.iter().any(|d| !(d.is_finite() && *d > 0.0)) {
        return Err(AcousticError::InvalidArgument(
            "room dimensions must be positive",
        ));
    }
    if alphas
        .iter()
        .any(|a| !(a.is_finite() && (0.0..=1.0).contains(a)))
    {
        return Err(AcousticError::InvalidArgument(
            "absorption must lie in [0, 1]",
        ));
    }
    Ok(())
}

/// Wall areas in the order used by the ray and ISM modules.
fn wall_areas(dims: [f64; 3]) -> [f64; 6] {
    let [lx, ly, lz] = dims;
    [ly * lz, ly * lz, lx * lz, lx * lz, lx * ly, lx * ly]
}

/// Sabine reverberation time `T = 0.161 V / A` and its gradient.
///
/// Returns `(T, ∂T/∂(Lx, Ly, Lz), ∂T/∂α₀…α₅)`. With `A = Σ Sᵢ αᵢ`:
/// `∂T/∂αᵢ = −T·Sᵢ / A` and `∂T/∂L = 0.161·∂V/∂L / A − T·(∂A/∂L) / A`.
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] for non-positive dimensions, an
/// absorption outside `[0, 1]`, or a room with no absorption (A = 0).
pub fn sabine_with_gradient(
    dims: [f64; 3],
    alphas: [f64; 6],
) -> Result<(f64, [f64; 3], [f64; 6]), AcousticError> {
    check_room(dims, &alphas)?;
    let [lx, ly, lz] = dims;
    let areas = wall_areas(dims);
    let a_total: f64 = areas.iter().zip(&alphas).map(|(s, a)| s * a).sum();
    if a_total <= 0.0 {
        return Err(AcousticError::InvalidArgument(
            "room has no absorption; Sabine time is infinite",
        ));
    }
    let volume = lx * ly * lz;
    let t = SABINE_CONSTANT * volume / a_total;

    // Per-axis absorption sums: a_x multiplies the x-face area Ly·Lz, etc.
    let a_x = alphas[0] + alphas[1];
    let a_y = alphas[2] + alphas[3];
    let a_z = alphas[4] + alphas[5];
    let d_a_dlx = a_y * lz + a_z * ly;
    let d_a_dly = a_x * lz + a_z * lx;
    let d_a_dlz = a_x * ly + a_y * lx;
    let d_v_dlx = ly * lz;
    let d_v_dly = lx * lz;
    let d_v_dlz = lx * ly;
    let grad_dims = [
        SABINE_CONSTANT * d_v_dlx / a_total - t * d_a_dlx / a_total,
        SABINE_CONSTANT * d_v_dly / a_total - t * d_a_dly / a_total,
        SABINE_CONSTANT * d_v_dlz / a_total - t * d_a_dlz / a_total,
    ];
    let mut grad_alpha = [0.0; 6];
    for i in 0..6 {
        grad_alpha[i] = -t * areas[i] / a_total;
    }
    Ok((t, grad_dims, grad_alpha))
}

/// Uniform absorption coefficient that gives Sabine time `target_t60_s` in a
/// room of the given size.
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] for invalid dimensions or a
/// target, and [`AcousticError::InvalidGeometry`] if the required uniform
/// absorption lies outside `(0, 1]` (the target cannot be met with a single
/// coefficient).
pub fn absorption_for_rt60(dims: [f64; 3], target_t60_s: f64) -> Result<f64, AcousticError> {
    if !(target_t60_s.is_finite() && target_t60_s > 0.0) {
        return Err(AcousticError::InvalidArgument(
            "target T60 must be positive",
        ));
    }
    check_room(dims, &[0.0; 6])?;
    let [lx, ly, lz] = dims;
    let volume = lx * ly * lz;
    let surface: f64 = wall_areas(dims).iter().sum();
    let alpha = SABINE_CONSTANT * volume / (target_t60_s * surface);
    if alpha > 1.0 {
        return Err(AcousticError::InvalidGeometry(
            "target T60 is too short for any uniform absorption ≤ 1",
        ));
    }
    Ok(alpha)
}

/// Scales a per-wall absorption distribution by one factor so the Sabine time
/// equals `target_t60_s`. The relative wall weights are kept, which makes this
/// a layout-preserving inverse design.
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] for invalid inputs, a
/// distribution with no absorption, or a target that needs a scale that puts
/// some coefficient above 1.
pub fn scale_absorption_to_rt60(
    dims: [f64; 3],
    alphas: [f64; 6],
    target_t60_s: f64,
) -> Result<[f64; 6], AcousticError> {
    if !(target_t60_s.is_finite() && target_t60_s > 0.0) {
        return Err(AcousticError::InvalidArgument(
            "target T60 must be positive",
        ));
    }
    let (t0, _, _) = sabine_with_gradient(dims, alphas)?;
    // T ∝ 1/s when all α are multiplied by s.
    let s = t0 / target_t60_s;
    let mut scaled = [0.0; 6];
    for (out, a) in scaled.iter_mut().zip(&alphas) {
        *out = a * s;
    }
    if scaled.iter().any(|a| *a > 1.0) {
        return Err(AcousticError::InvalidGeometry(
            "scaled absorption exceeds 1 on some wall; target is out of reach for this layout",
        ));
    }
    Ok(scaled)
}

/// Summary of a Monte Carlo reverberation-time distribution.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Uncertainty {
    /// Mean of the samples, seconds.
    pub mean_s: f64,
    /// Sample standard deviation, seconds.
    pub std_s: f64,
    /// 5th percentile, seconds.
    pub p05_s: f64,
    /// 95th percentile, seconds.
    pub p95_s: f64,
}

/// Monte Carlo spread of the Sabine time when each wall's absorption is
/// independently normal with mean `alpha_mean` and standard deviation
/// `alpha_sd`, clipped to `[0.001, 0.999]`. The generator is seeded, so the
/// result is reproducible.
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] for invalid dimensions, a mean
/// outside `[0, 1]`, a negative spread, or zero samples.
pub fn rt60_uncertainty(
    dims: [f64; 3],
    alpha_mean: f64,
    alpha_sd: f64,
    samples: usize,
    seed: u64,
) -> Result<Uncertainty, AcousticError> {
    check_room(dims, &[alpha_mean; 6])?;
    if !(alpha_sd.is_finite() && alpha_sd >= 0.0) {
        return Err(AcousticError::InvalidArgument(
            "spread must be finite and ≥ 0",
        ));
    }
    if samples == 0 {
        return Err(AcousticError::InvalidArgument(
            "at least one sample is needed",
        ));
    }
    let mut rng = Xoshiro256::new(seed);
    let mut times: Vec<f64> = Vec::with_capacity(samples);
    for _ in 0..samples {
        let mut alphas = [0.0; 6];
        for a in alphas.iter_mut() {
            *a = (alpha_mean + alpha_sd * rng.next_normal()).clamp(0.001, 0.999);
        }
        let (t, _, _) = sabine_with_gradient(dims, alphas)?;
        times.push(t);
    }
    let n = times.len() as f64;
    let mean = times.iter().sum::<f64>() / n;
    let var = times.iter().map(|t| (t - mean) * (t - mean)).sum::<f64>() / n;
    times.sort_by(|a, b| a.partial_cmp(b).unwrap_or(core::cmp::Ordering::Equal));
    let pick = |q: f64| times[libm::round((times.len() - 1) as f64 * q) as usize];
    Ok(Uncertainty {
        mean_s: mean,
        std_s: libm::sqrt(var),
        p05_s: pick(0.05),
        p95_s: pick(0.95),
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    const DIMS: [f64; 3] = [10.0, 8.0, 6.0];
    const ALPHAS: [f64; 6] = [0.10, 0.12, 0.20, 0.15, 0.08, 0.30];

    fn t_of(dims: [f64; 3], alphas: [f64; 6]) -> f64 {
        sabine_with_gradient(dims, alphas).unwrap().0
    }

    #[test]
    fn sabine_matches_the_library_room() {
        // The ray crate's ShoeboxRoom uses the same model; compare directly.
        let (t, _, _) = sabine_with_gradient(DIMS, [0.15; 6]).unwrap();
        let room = crate::ism::ShoeboxRoom::uniform(
            DIMS,
            crate::material::BoundaryMaterial::new(
                "uniform",
                alloc::vec![(125.0, 0.15), (4000.0, 0.15)],
                alloc::vec![(125.0, 0.0), (4000.0, 0.0)],
            )
            .unwrap(),
        )
        .unwrap();
        assert!((t - room.sabine_rt60(1000.0)).abs() < 1e-9);
    }

    #[test]
    fn gradients_match_finite_differences() {
        let (_, grad_dims, grad_alpha) = sabine_with_gradient(DIMS, ALPHAS).unwrap();
        let h = 1e-6;
        for axis in 0..3 {
            let mut plus = DIMS;
            let mut minus = DIMS;
            plus[axis] += h;
            minus[axis] -= h;
            let fd = (t_of(plus, ALPHAS) - t_of(minus, ALPHAS)) / (2.0 * h);
            assert!(
                (grad_dims[axis] - fd).abs() < 1e-6 * fd.abs().max(1.0),
                "dim {axis}: {} vs {fd}",
                grad_dims[axis]
            );
        }
        for i in 0..6 {
            let mut plus = ALPHAS;
            let mut minus = ALPHAS;
            plus[i] += h;
            minus[i] -= h;
            let fd = (t_of(DIMS, plus) - t_of(DIMS, minus)) / (2.0 * h);
            assert!(
                (grad_alpha[i] - fd).abs() < 1e-6 * fd.abs().max(1.0),
                "alpha {i}: {} vs {fd}",
                grad_alpha[i]
            );
        }
    }

    #[test]
    fn uniform_inverse_round_trips() {
        let alpha = absorption_for_rt60(DIMS, 1.0).unwrap();
        let t = t_of(DIMS, [alpha; 6]);
        assert!((t - 1.0).abs() < 1e-12, "{t}");
    }

    #[test]
    fn scaling_keeps_the_layout_and_hits_the_target() {
        let target = 0.8;
        let scaled = scale_absorption_to_rt60(DIMS, ALPHAS, target).unwrap();
        assert!((t_of(DIMS, scaled) - target).abs() < 1e-12);
        // Relative wall weights are unchanged.
        let ratio = scaled[5] / scaled[0];
        assert!((ratio - ALPHAS[5] / ALPHAS[0]).abs() < 1e-12);
    }

    #[test]
    fn unreachable_targets_are_reported() {
        // A very short target needs α > 1 everywhere.
        assert!(absorption_for_rt60(DIMS, 0.01).is_err());
        assert!(scale_absorption_to_rt60(DIMS, ALPHAS, 0.01).is_err());
        assert!(sabine_with_gradient(DIMS, [0.0; 6]).is_err());
    }

    #[test]
    fn uncertainty_is_reproducible_and_brackets_the_mean() {
        let a = rt60_uncertainty(DIMS, 0.2, 0.02, 2000, 9).unwrap();
        let b = rt60_uncertainty(DIMS, 0.2, 0.02, 2000, 9).unwrap();
        assert_eq!(a, b);
        assert!(a.std_s > 0.0);
        assert!(a.p05_s < a.mean_s && a.mean_s < a.p95_s, "{a:?}");
        // With no spread the distribution collapses to the nominal time.
        let nominal = rt60_uncertainty(DIMS, 0.2, 0.0, 10, 1).unwrap();
        assert!(nominal.std_s < 1e-12);
        assert!((nominal.mean_s - t_of(DIMS, [0.2; 6])).abs() < 1e-9);
    }
}
