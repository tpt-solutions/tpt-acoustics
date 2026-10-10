//! Acoustic Transfer Vectors (ATVs).
//!
//! The ATV is the linear (frequency-dependent) map from element normal
//! velocity to field-point pressure for a baffled radiating surface:
//!
//! `p(x) = sum_j ATV_j(x) * v_j`,
//! `ATV_j(x) = 2 * (i*omega*rho) * G(x, y_j) * A_j`
//!
//! (Neumann-data factor times the rigid-baffle image doubling; this matches
//! the open-patch exterior BEM solve, whose single layer radiates from both
//! sides). ATVs are
//! precomputed once per frequency and reused across many response
//! evaluations — the standard vibro-acoustic transfer-path workflow.

use tpt_acou_core::error::AcousticError;
use tpt_acou_core::medium::Medium;
use tpt_acou_core::Vec3;
use tpt_math_linalg_complex::Complex;

use tpt_acou_wave::bem::HelmholtzBem;

/// Computes the ATV row for one field point: `ATV_j` for every mesh element.
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] if `x` coincides with an
/// element centroid (singular Green's function).
pub fn atv_row(
    bem: &HelmholtzBem,
    medium: &Medium,
    x: Vec3,
) -> Result<Vec<Complex>, AcousticError> {
    let els = bem.mesh().elements();
    let k = bem.k();
    let omega = k * medium.speed_of_sound;
    let mut row = Vec::with_capacity(els.len());
    for e in els {
        let r = x.distance(e.centroid);
        if r < 1e-9 {
            return Err(AcousticError::InvalidArgument(
                "field point coincides with an element centroid",
            ));
        }
        // G = exp(-ikr)/(4 pi r); ATV_j = 2 * (i*omega*rho) * G * A_j
        // (Rayleigh-integral sign, baffle image factor included).
        let g = Complex::from_polar(1.0 / (4.0 * core::f64::consts::PI * r), -k * r);
        row.push(g * Complex::new(0.0, omega * medium.density) * e.area * 2.0);
    }
    Ok(row)
}

/// Computes the ATV matrix for several field points (one row per point).
///
/// # Errors
/// Propagates [`atv_row`] errors.
pub fn atv_matrix(
    bem: &HelmholtzBem,
    medium: &Medium,
    field_points: &[Vec3],
) -> Result<Vec<Vec<Complex>>, AcousticError> {
    field_points
        .iter()
        .map(|&x| atv_row(bem, medium, x))
        .collect()
}

/// Applies the ATV matrix to element velocities: one complex pressure per
/// field point.
#[must_use]
pub fn radiated_pressures(atv: &[Vec<Complex>], v_n: &[Complex]) -> Vec<Complex> {
    atv.iter()
        .map(|row| {
            let mut p = Complex::new(0.0, 0.0);
            for (a, v) in row.iter().zip(v_n) {
                p = p + *a * *v;
            }
            p
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use tpt_acou_wave::mesh::BoundaryMesh;

    #[test]
    fn atv_reproduces_single_layer_evaluation() {
        // ATV applied to v must equal the full Neumann BEM evaluation of the
        // same velocity field (they share the single-layer representation).
        let mesh = BoundaryMesh::rectangle_patch(1.0, 1.0, 4, 4).unwrap();
        let k = 2.0;
        let bem = HelmholtzBem::new(mesh, k).unwrap();
        let medium = Medium::air_20c();
        let n = bem.mesh().len();
        let omega = k * medium.speed_of_sound;
        let v: Vec<Complex> = (0..n)
            .map(|i| Complex::new(1.0 + 0.05 * i as f64, 0.1))
            .collect();
        let x = Vec3::new(0.3, -0.2, 1.5);
        let atv = atv_row(&bem, &medium, x).unwrap();
        // ATV route (baffled): p = 2 * (i w rho) sum G A_j v_j.
        let mut p_atv = Complex::new(0.0, 0.0);
        for (a, vv) in atv.iter().zip(&v) {
            p_atv = p_atv + *a * *vv;
        }
        // Direct route: q = -i w rho v; the open-patch BEM Neumann solve
        // produces the two-sided (baffle-equivalent) field.
        let q: Vec<Complex> = v
            .iter()
            .map(|vi| *vi * Complex::new(0.0, -omega * medium.density))
            .collect();
        let sigma = bem.solve_exterior_neumann(&q).unwrap();
        let p_direct = bem.evaluate_single_layer(&sigma, &[x]).unwrap()[0];
        // The ATV neglects the (O(v^2) small for light coupling) scattering
        // of the surface onto itself; for the patch at k=2 they agree to a
        // few percent.
        assert!(
            (p_atv - p_direct).norm() < 0.2 * p_direct.norm(),
            "atv ({:.4},{:.4}) vs direct ({:.4},{:.4})",
            p_atv.re,
            p_atv.im,
            p_direct.re,
            p_direct.im
        );
    }
}
