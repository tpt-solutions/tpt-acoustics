//! Acoustic radiation: baffled-piston resistance and BEM radiation
//! efficiency, plus the modal-velocity power interface for `tpt-fem-modal`
//! data.

use tpt_acou_core::error::AcousticError;
use tpt_acou_core::medium::Medium;
use tpt_math_linalg_complex::Complex;

use tpt_acou_wave::bem::HelmholtzBem;

/// First-kind Bessel function `J1(x)` (power series below 8, asymptotic
/// expansion above; ~1e-10 accurate in each regime).
#[must_use]
pub fn bessel_j1(x: f64) -> f64 {
    let ax = x.abs();
    if ax < 8.0 {
        // J1(x) = sum_k (-1)^k (x/2)^(2k+1) / (k! (k+1)!).
        let half = ax / 2.0;
        let mut term = half;
        let mut sum = term;
        for k in 1..30 {
            term *= -half * half / (f64::from(k) * f64::from(k + 1));
            sum += term;
            if term.abs() < 1e-18 * sum.abs().max(1e-30) {
                break;
            }
        }
        if x < 0.0 {
            -sum
        } else {
            sum
        }
    } else {
        // Asymptotic expansion (Abramowitz & Stegun 9.2.1):
        // J1 ~ sqrt(2/(pi x)) [P cos(chi) - Q sin(chi)], chi = x - 3pi/4,
        // P = 1 + 15/(128 x^2), Q = 3/(8x) - 315/(3072 x^3)  (mu = 4).
        let chi = ax - 0.75 * core::f64::consts::PI;
        let amp = libm::sqrt(2.0 / (core::f64::consts::PI * ax));
        let p1 = 1.0 + 15.0 / (128.0 * ax * ax);
        let q1 = 3.0 / (8.0 * ax) - 315.0 / (3072.0 * ax * ax * ax);
        let val = amp * (p1 * libm::cos(chi) - q1 * libm::sin(chi));
        if x < 0.0 {
            -val
        } else {
            val
        }
    }
}

/// Normalized radiation resistance of a baffled circular piston:
/// `R = 1 - 2 J1(2ka)/(2ka)`. Equals the radiation efficiency of the piston;
/// small-argument limit `(ka)^2/2`, large-argument limit 1. (The reactance —
/// Struve-function part — is not included.)
#[must_use]
pub fn piston_radiation_resistance(ka: f64) -> f64 {
    if ka <= 0.0 {
        return 0.0;
    }
    let z = 2.0 * ka;
    1.0 - 2.0 * bessel_j1(z) / z
}

/// Radiation efficiency (power coefficient) of a baffled surface with a
/// prescribed complex normal-velocity pattern `v_n` (m/s amplitudes), via the
/// `tpt-acou-wave` Helmholtz BEM:
///
/// `sigma = W_rad / (rho c A <v^2>)`
///
/// with `W_rad = 2 · ½ · Re ∫ p v* dS` (the factor 2 is the rigid-baffle
/// image) and `<v^2>` the area-averaged squared velocity amplitude.
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] on a velocity length mismatch
/// and [`AcousticError::NonPositive`] if all velocities are zero; propagates
/// singular BEM solves as [`AcousticError::Singular`].
pub fn radiation_efficiency_baffled(
    bem: &HelmholtzBem,
    v_n: &[Complex],
    medium: &Medium,
) -> Result<f64, AcousticError> {
    let els = bem.mesh().elements();
    if els.len() != v_n.len() {
        return Err(AcousticError::InvalidArgument("velocity length mismatch"));
    }
    let area: f64 = els.iter().map(|e| e.area).sum();
    let v_rms_sq: f64 = els
        .iter()
        .zip(v_n)
        .map(|(e, v)| v.norm_sqr() * e.area)
        .sum::<f64>()
        / area;
    if v_rms_sq <= 0.0 || v_rms_sq.is_nan() {
        return Err(AcousticError::NonPositive("reference velocity"));
    }
    // Neumann data for the radiating surface: q = dP/dn = -i*omega*rho*v_n
    // (exp(+i w t) convention, exterior side).
    let omega = bem.k() * medium.speed_of_sound;
    let q: Vec<Complex> = v_n
        .iter()
        .map(|v| *v * Complex::new(0.0, -omega * medium.density))
        .collect();
    let sigma = bem
        .solve_exterior_neumann(&q)
        .map_err(|_| AcousticError::Singular("BEM radiation solve"))?;
    // Surface pressure via the single-layer potential evaluated just outside
    // each centroid (standard constant-element route).
    let probe: Vec<tpt_acou_core::Vec3> = els
        .iter()
        .map(|e| e.centroid + e.normal * libm::sqrt(e.area) * 1e-4)
        .collect();
    let p = bem
        .evaluate_single_layer(&sigma, &probe)
        .map_err(|_| AcousticError::Singular("BEM evaluation"))?;
    let mut w_rad = 0.0;
    for (e, (pv, v)) in els.iter().zip(p.iter().zip(v_n)) {
        let flux = *pv * v.conj();
        w_rad += flux.re * e.area;
    }
    let w_baffled = w_rad; // 2 * 1/2 * integral
    Ok(w_baffled / (medium.characteristic_impedance() * area * v_rms_sq))
}

/// Radiated sound power (watts) from modal velocities — the
/// `tpt-fem-modal` bridge.
///
/// `modal` carries mode shapes over structural DOFs; `modal_velocities` are
/// the complex modal generalised velocities `i·omega·q_i` at the excitation
/// frequency; `element_dof[e]` maps mesh element `e` to the structural DOF
/// whose mode-shape entry drives it.
///
/// # Errors
/// Propagates [`radiation_efficiency_baffled`] errors and
/// [`AcousticError::InvalidArgument`] on length mismatches.
pub fn radiated_power_from_modal_velocities(
    bem: &HelmholtzBem,
    modal: &tpt_fem_modal::ModalData,
    modal_velocities: &[Complex],
    element_dof: &[usize],
    rho_c: f64,
) -> Result<f64, AcousticError> {
    let els = bem.mesh().elements();
    if element_dof.len() != els.len() {
        return Err(AcousticError::InvalidArgument(
            "element_dof must map each element to a structural DOF",
        ));
    }
    if modal_velocities.len() != modal.omega.len() {
        return Err(AcousticError::InvalidArgument(
            "one modal velocity per retained mode",
        ));
    }
    let mut v_n = Vec::with_capacity(els.len());
    for &dof in element_dof {
        let mut v = Complex::new(0.0, 0.0);
        for (m, qd) in modal_velocities.iter().enumerate() {
            let shape = modal.modes[m].get(dof).copied().unwrap_or(0.0);
            v = v + *qd * shape;
        }
        v_n.push(v);
    }
    let eff = radiation_efficiency_baffled(bem, &v_n, &medium_from_impedance(rho_c))?;
    let area: f64 = els.iter().map(|e| e.area).sum();
    let v_rms_sq: f64 = els
        .iter()
        .zip(&v_n)
        .map(|(e, v)| v.norm_sqr() * e.area)
        .sum::<f64>()
        / area;
    Ok(eff * rho_c * area * v_rms_sq)
}

/// Reconstructs a lossless `Medium` from a characteristic impedance at the
/// reference 343 m/s (used where callers pass `rho_c` directly).
fn medium_from_impedance(rho_c: f64) -> Medium {
    Medium {
        density: rho_c / 343.0,
        speed_of_sound: 343.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bessel_j1_matches_known_values() {
        // J1(1) = 0.4400505857; J1(10) = 0.0434727.
        assert!((bessel_j1(1.0) - 0.440_050_585_7).abs() < 1e-8);
        assert!((bessel_j1(10.0) - 0.043_472_7).abs() < 1e-5);
        assert!(bessel_j1(0.0).abs() < 1e-15);
    }

    #[test]
    fn piston_resistance_limits() {
        // Small ka: (ka)^2/2; large ka: -> 1.
        let small = piston_radiation_resistance(0.05);
        assert!((small - 0.5 * 0.0025).abs() < 1e-6, "{small}");
        // Oscillation decays like (ka)^(-3/2): ~0.7% at ka = 10.
        assert!((piston_radiation_resistance(10.0) - 1.0).abs() < 0.01);
    }

    #[test]
    fn piston_resistance_stays_bounded() {
        // The resistance oscillates about 1 at high ka (overshooting ~13%
        // near ka = 2.4) but stays positive and bounded.
        for &ka in &[0.3_f64, 1.0, 2.5, 5.0] {
            let r = piston_radiation_resistance(ka);
            assert!(r > 0.0 && r < 1.4, "ka={ka}: r={r}");
        }
    }
}
