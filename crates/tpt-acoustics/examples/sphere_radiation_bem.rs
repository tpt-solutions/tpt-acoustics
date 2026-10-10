//! Exterior acoustic radiation of a pulsating sphere, boundary element solution.
//!
//! A sphere of radius 1 m pulses with unit surface velocity at ka = 1.3. The
//! exterior Dirichlet problem is solved with the BEM, prescribing the exact
//! surface pressure, and the field is evaluated at r = 2a. That point is
//! compared with the closed-form solution
//!
//! ```text
//! p(r) = ρc · v_a · (i k a²) / (1 + i k a) · (a / r) · e^{−i k (r − a)}
//! ```
//!
//! The BEM result must agree with the closed form within 5 %.
//!
//! Run with `cargo run --release --example sphere_radiation_bem`.

use tpt_acoustics::core::{AcousticError, Vec3};
use tpt_acoustics::wave::{BoundaryMesh, HelmholtzBem};
use tpt_math_linalg_complex::Complex;

/// Closed-form pressure of a pulsating sphere (time factor e^{+iωt}).
fn pulsating_sphere_pressure(k: f64, a: f64, v_a: Complex, r: f64) -> Complex {
    let rho_c = 413.0;
    let g = v_a * rho_c * (Complex::new(0.0, k) * a * a) / Complex::new(1.0, k * a);
    (g / Complex::new(r, 0.0)) * Complex::from_polar(1.0, -k * (r - a))
}

fn main() -> Result<(), AcousticError> {
    let a = 1.0;
    let k = 1.3;
    let v_a = Complex::new(1.0, 0.0);

    let mesh = BoundaryMesh::sphere(a, 10, 16)?;
    let bem = HelmholtzBem::new(mesh, k).expect("BEM setup");

    // Prescribe the exact surface pressure on every element.
    let p_surface = pulsating_sphere_pressure(k, a, v_a, a);
    let data = vec![p_surface; bem.mesh().len()];
    let density = bem.solve_exterior_dirichlet(&data).expect("exterior solve");

    let field_point = [Vec3::new(0.0, 0.0, 2.0 * a)];
    let p_bem = bem
        .evaluate_double_layer(&density, &field_point)
        .expect("evaluate field")[0];
    let p_exact = pulsating_sphere_pressure(k, a, v_a, 2.0 * a);
    let rel = (p_bem - p_exact).norm() / p_exact.norm();

    println!("BEM   p(2a) = {:.4} {:+.4}i Pa", p_bem.re, p_bem.im);
    println!("exact p(2a) = {:.4} {:+.4}i Pa", p_exact.re, p_exact.im);
    println!("relative error: {:.2} %", 100.0 * rel);
    assert!(rel < 0.05, "BEM disagrees with the closed form by {rel:.3}");
    Ok(())
}
