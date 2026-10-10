//! Collocation Boundary Element Method for the Helmholtz equation.
//!
//! Constant flat triangular elements with centroid collocation. The key
//! structural property of this discretization: the *self* terms of the double
//! layer `∂G/∂n_y` and of the adjoint double layer `∂G/∂n_x` vanish exactly,
//! because the kernel is antisymmetric in the element plane (`(x − y)·n = 0`
//! there). The four canonical problems therefore need **no singular
//! quadrature at all**:
//!
//! | problem | representation | boundary equation |
//! | --- | --- | --- |
//! | exterior Dirichlet | `p = Dμ` | `(−½I + D)μ = g` |
//! | interior Dirichlet | `p = Dμ` | `(+½I + D)μ = g` |
//! | exterior Neumann | `p = Sσ` | `(−½I + K′)σ = q` |
//! | interior Neumann | `p = Sσ` | `(+½I + K′)σ = q` |
//!
//! (Jump signs verified against the closed-form pulsating-sphere and
//! constant-value sphere solutions; with outward normals and
//! `G = e^{−ikr}/4πr`, `γ^+ Dμ = (−½ + K)μ` and
//! `γ^- Dμ = (+½ + K)μ`.)
//!
//! where `D_ij = ∂G(x_i, y_j)/∂n_y · A_j`, `K′_ij = ∂G(x_i, y_j)/∂n_x · A_j`
//! (`K′_ii = D_ii = 0`), and `G` is evaluated off-surface only.
//!
//! Excluding the (well-documented) interior-resonance non-uniqueness of the
//! exterior formulations, all four systems are diagonally dominant and solved
//! with the complex partial-pivot LU from `tpt-math-linalg-complex`.

use alloc::vec::Vec;

use tpt_acou_core::error::AcousticError;
use tpt_math_linalg_complex::{Complex, ComplexDMatrix, ComplexDVector};

use crate::mesh::BoundaryMesh;
use tpt_acou_core::Vec3;

/// Failure modes of a BEM solve.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BemError {
    /// The LU factorization found a numerically singular matrix (typically an
    /// interior eigenfrequency of the complementary domain).
    Singular,
    /// The boundary data length does not match the mesh element count.
    DataLengthMismatch,
}

impl core::fmt::Display for BemError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            BemError::Singular => f.write_str("singular BEM system (interior resonance?)"),
            BemError::DataLengthMismatch => {
                f.write_str("boundary data length must equal the element count")
            }
        }
    }
}

impl core::error::Error for BemError {}

impl From<AcousticError> for BemError {
    fn from(_value: AcousticError) -> Self {
        // Degenerate-mesh and non-positive-wavenumber errors both surface as
        // construction failures; the mapping keeps `?` ergonomic at call sites.
        BemError::Singular
    }
}

/// A `k`-bound BEM operator bundle over a fixed mesh.
pub struct HelmholtzBem {
    mesh: BoundaryMesh,
    k: f64,
}

impl HelmholtzBem {
    /// Binds the mesh and wavenumber (rad/m).
    ///
    /// # Errors
    /// Returns [`AcousticError::NonPositive`] for `k <= 0`.
    pub fn new(mesh: BoundaryMesh, k: f64) -> Result<Self, AcousticError> {
        if k <= 0.0 {
            return Err(AcousticError::NonPositive("wavenumber"));
        }
        Ok(Self { mesh, k })
    }

    /// The mesh.
    #[must_use]
    pub const fn mesh(&self) -> &BoundaryMesh {
        &self.mesh
    }

    /// The wavenumber.
    #[must_use]
    pub const fn k(&self) -> f64 {
        self.k
    }

    /// Free-space Green's function `G(r) = e^{−ikr}/(4πr)`.
    fn green(&self, r: f64) -> Complex {
        let mag = 1.0 / (4.0 * core::f64::consts::PI * r);
        Complex::from_polar(mag, -self.k * r)
    }

    /// Double-layer matrix `D_ij = ∂G(x_i, y_j)/∂n_j · A_j`, self terms zero.
    ///
    /// `∂G/∂n_y = −G·(ik + 1/r)·((x−y)·n_j)/r`.
    #[must_use]
    pub fn double_layer_matrix(&self) -> ComplexDMatrix {
        let n = self.mesh.len();
        let els = self.mesh.elements();
        let mut data = Vec::with_capacity(n * n);
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    data.push(Complex::new(0.0, 0.0));
                    continue;
                }
                let x = els[i].centroid;
                let (y, nj, aj) = (els[j].centroid, els[j].normal, els[j].area);
                let d = x - y;
                let r = d.norm();
                let g = self.green(r);
                let dn = d.dot(nj) / r;
                data.push(
                    g * (Complex::new(0.0, self.k) + Complex::new(1.0 / r, 0.0)) * (-dn * aj),
                );
            }
        }
        ComplexDMatrix::from_row_slice(n, n, &data)
    }

    /// Adjoint double-layer matrix `K′_ij = ∂G(x_i, y_j)/∂n_x · A_j`, self
    /// terms zero. `∂G/∂n_x = −G·(ik + 1/r)·((x−y)·n_i)/r`.
    #[must_use]
    pub fn adjoint_double_layer_matrix(&self) -> ComplexDMatrix {
        let n = self.mesh.len();
        let els = self.mesh.elements();
        let mut data = Vec::with_capacity(n * n);
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    data.push(Complex::new(0.0, 0.0));
                    continue;
                }
                let x = els[i].centroid;
                let (y, ni, aj) = (els[j].centroid, els[i].normal, els[j].area);
                let d = x - y;
                let r = d.norm();
                let g = self.green(r);
                let dn = d.dot(ni) / r;
                data.push(
                    g * (Complex::new(0.0, self.k) + Complex::new(1.0 / r, 0.0)) * (-dn * aj),
                );
            }
        }
        ComplexDMatrix::from_row_slice(n, n, &data)
    }

    /// Single-layer potential `p(x) = Σ_j G(x, y_j)·σ_j·A_j` evaluated at
    /// field points *off* the boundary.
    pub fn evaluate_single_layer(
        &self,
        sigma: &[Complex],
        points: &[Vec3],
    ) -> Result<Vec<Complex>, BemError> {
        self.check_data(sigma)?;
        let els = self.mesh.elements();
        Ok(points
            .iter()
            .map(|&x| {
                let mut acc = Complex::new(0.0, 0.0);
                for (j, e) in els.iter().enumerate() {
                    let r = x.distance(e.centroid);
                    acc = acc + self.green(r) * sigma[j] * e.area;
                }
                acc
            })
            .collect())
    }

    /// Double-layer potential `p(x) = Σ_j ∂G(x, y_j)/∂n_j·μ_j·A_j` at field
    /// points off the boundary.
    pub fn evaluate_double_layer(
        &self,
        mu: &[Complex],
        points: &[Vec3],
    ) -> Result<Vec<Complex>, BemError> {
        self.check_data(mu)?;
        let els = self.mesh.elements();
        Ok(points
            .iter()
            .map(|&x| {
                let mut acc = Complex::new(0.0, 0.0);
                for (j, e) in els.iter().enumerate() {
                    let d = x - e.centroid;
                    let r = d.norm();
                    let g = self.green(r);
                    let dn = d.dot(e.normal) / r;
                    acc = acc
                        + g * (Complex::new(0.0, self.k) + Complex::new(1.0 / r, 0.0))
                            * (-dn)
                            * e.area
                            * mu[j];
                }
                acc
            })
            .collect())
    }

    /// Exterior Dirichlet solve: `p = Dμ` with `(−½I + D)μ = g`.
    ///
    /// # Errors
    /// [`BemError::DataLengthMismatch`] or [`BemError::Singular`].
    pub fn solve_exterior_dirichlet(&self, g: &[Complex]) -> Result<Vec<Complex>, BemError> {
        self.check_data(g)?;
        let n = self.mesh.len();
        let d = self.double_layer_matrix();
        let mut data = Vec::with_capacity(n * n);
        for i in 0..n {
            for j in 0..n {
                let mut v = d[(i, j)];
                if i == j {
                    v = v + Complex::new(-0.5, 0.0);
                }
                data.push(v);
            }
        }
        let a = ComplexDMatrix::from_row_slice(n, n, &data);
        let mu = a
            .solve(&ComplexDVector::from_vec(g.to_vec()))
            .map_err(|_| BemError::Singular)?;
        Ok((0..n).map(|i| mu[i]).collect())
    }

    /// Interior Dirichlet solve: `p = Dμ` with `(+½I + D)μ = g`.
    ///
    /// # Errors
    /// See [`HelmholtzBem::solve_exterior_dirichlet`].
    pub fn solve_interior_dirichlet(&self, g: &[Complex]) -> Result<Vec<Complex>, BemError> {
        self.check_data(g)?;
        let n = self.mesh.len();
        let d = self.double_layer_matrix();
        let mut data = Vec::with_capacity(n * n);
        for i in 0..n {
            for j in 0..n {
                let mut v = d[(i, j)];
                if i == j {
                    v = v + Complex::new(0.5, 0.0);
                }
                data.push(v);
            }
        }
        let a = ComplexDMatrix::from_row_slice(n, n, &data);
        let mu = a
            .solve(&ComplexDVector::from_vec(g.to_vec()))
            .map_err(|_| BemError::Singular)?;
        Ok((0..n).map(|i| mu[i]).collect())
    }

    /// Exterior Neumann solve: `p = Sσ` with `(−½I + K′)σ = q`, where `q` is
    /// the outward normal derivative ∂p/∂n (m·s⁻¹ scale, i.e. `iωρ·v_n` for
    /// the acoustic velocity convention of `tpt-acou-core`).
    ///
    /// # Errors
    /// See [`HelmholtzBem::solve_exterior_dirichlet`].
    pub fn solve_exterior_neumann(&self, q: &[Complex]) -> Result<Vec<Complex>, BemError> {
        self.check_data(q)?;
        let n = self.mesh.len();
        let kp = self.adjoint_double_layer_matrix();
        let mut data = Vec::with_capacity(n * n);
        for i in 0..n {
            for j in 0..n {
                let mut v = kp[(i, j)];
                if i == j {
                    v = v + Complex::new(-0.5, 0.0);
                }
                data.push(v);
            }
        }
        let a = ComplexDMatrix::from_row_slice(n, n, &data);
        let sigma = a
            .solve(&ComplexDVector::from_vec(q.to_vec()))
            .map_err(|_| BemError::Singular)?;
        Ok((0..n).map(|i| sigma[i]).collect())
    }

    /// Interior Neumann solve: `p = Sσ` with `(+½I + K′)σ = q`.
    ///
    /// # Errors
    /// See [`HelmholtzBem::solve_exterior_dirichlet`].
    pub fn solve_interior_neumann(&self, q: &[Complex]) -> Result<Vec<Complex>, BemError> {
        self.check_data(q)?;
        let n = self.mesh.len();
        let kp = self.adjoint_double_layer_matrix();
        let mut data = Vec::with_capacity(n * n);
        for i in 0..n {
            for j in 0..n {
                let mut v = kp[(i, j)];
                if i == j {
                    v = v + Complex::new(0.5, 0.0);
                }
                data.push(v);
            }
        }
        let a = ComplexDMatrix::from_row_slice(n, n, &data);
        let sigma = a
            .solve(&ComplexDVector::from_vec(q.to_vec()))
            .map_err(|_| BemError::Singular)?;
        Ok((0..n).map(|i| sigma[i]).collect())
    }

    /// Radiated sound power `W = ½·Re(Σ p_i·v_n,i*·A_i)` in watts, from
    /// boundary pressure and normal velocity samples (RMS-compatible complex
    /// amplitudes).
    pub fn radiated_power(&self, p: &[Complex], v_n: &[Complex]) -> Result<f64, BemError> {
        self.check_data(p)?;
        self.check_data(v_n)?;
        let els = self.mesh.elements();
        let mut w = 0.0;
        for (i, e) in els.iter().enumerate() {
            let pv = p[i] * v_n[i].conj();
            w += pv.re * e.area;
        }
        Ok(0.5 * w)
    }

    fn check_data(&self, data: &[Complex]) -> Result<(), BemError> {
        if data.len() != self.mesh.len() {
            return Err(BemError::DataLengthMismatch);
        }
        Ok(())
    }
}

/// Full-sphere interior validation helper: `j0(x) = sin(x)/x` with the
/// removable singularity handled.
#[must_use]
pub fn spherical_bessel_j0(x: f64) -> f64 {
    if x.abs() < 1e-8 {
        1.0 - x * x / 6.0
    } else {
        libm::sin(x) / x
    }
}

/// `j1(x)` spherical Bessel function (series + recursion).
#[must_use]
pub fn spherical_bessel_j1(x: f64) -> f64 {
    if x.abs() < 1e-8 {
        x / 3.0
    } else {
        libm::sin(x) / (x * x) - libm::cos(x) / x
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::BoundaryMesh;
    use alloc::vec;

    const TOL: f64 = 1e-9;

    /// Analytic pulsating-sphere radiation (time factor e^{+iωt}):
    /// `p(r) = ρc·v_a·(ika²)/(1+ika) · (a/r)·e^{−ik(r−a)}`.
    fn pulsating_sphere_exact(k: f64, a: f64, v_a: Complex, r: f64) -> Complex {
        let rho_c = 413.0;
        let g = v_a * rho_c * (Complex::new(0.0, k) * a * a) / Complex::new(1.0, k * a);
        (g / Complex::new(r, 0.0)) * Complex::from_polar(1.0, -k * (r - a))
    }

    #[test]
    fn double_layer_self_terms_are_zero() {
        let mesh = BoundaryMesh::sphere(1.0, 6, 8).unwrap();
        let bem = HelmholtzBem::new(mesh, 2.0).unwrap();
        let d = bem.double_layer_matrix();
        for i in 0..bem.mesh().len() {
            assert!(d[(i, i)].norm() < 1e-15);
        }
        let kp = bem.adjoint_double_layer_matrix();
        for i in 0..bem.mesh().len() {
            assert!(kp[(i, i)].norm() < 1e-15);
        }
    }

    #[test]
    fn exterior_pulsating_sphere_neumann() {
        let a = 1.0;
        let k = 1.7;
        let v_a = Complex::new(1.0, 0.0);
        let mesh = BoundaryMesh::sphere(a, 10, 16).unwrap();
        let bem = HelmholtzBem::new(mesh, k).unwrap();
        // q = ∂p/∂r at r=a from the exact solution:
        // p(r) = C·e^{−ik(r−a)}/r, ∂p/∂r = C·e^{−ik(r−a)}·(−ik − 1/r).
        let c = pulsating_sphere_exact(k, a, v_a, a) * a; // C = p(a)·a
        let q_at_surface = |_: usize| {
            c * Complex::from_polar(1.0, 0.0)
                * (Complex::new(0.0, -k) + Complex::new(-1.0 / a, 0.0))
        };
        let n = bem.mesh().len();
        let q: Vec<Complex> = (0..n).map(q_at_surface).collect();
        let sigma = bem.solve_exterior_neumann(&q).unwrap();
        // Evaluate at r = 2a on the z axis and compare.
        let points = [Vec3::new(0.0, 0.0, 2.0 * a)];
        let p_bem = bem.evaluate_single_layer(&sigma, &points).unwrap();
        let p_exact = pulsating_sphere_exact(k, a, v_a, 2.0 * a);
        let rel = (p_bem[0] - p_exact).norm() / p_exact.norm();
        assert!(rel < 0.05, "BEM {p_bem:?} vs exact {p_exact:?} (rel {rel})");
    }

    #[test]
    fn interior_sphere_dirichlet() {
        // Interior of a sphere with constant boundary value 1:
        // p(r) = j0(kr)/j0(ka). Choose ka = 1.0 (no interior eigenfrequency:
        // those sit at ka = π, 2π, ...).
        let a = 1.0;
        let k = 1.0;
        let mesh = BoundaryMesh::sphere(a, 10, 16).unwrap();
        let bem = HelmholtzBem::new(mesh, k).unwrap();
        let n = bem.mesh().len();
        let g = vec![Complex::new(1.0, 0.0); n];
        let mu = bem.solve_interior_dirichlet(&g).unwrap();
        let points = [Vec3::new(0.0, 0.0, 0.0), Vec3::new(0.4, 0.0, 0.0)];
        // p = −Dμ for the interior representation.
        let p_bem = bem.evaluate_double_layer(&mu, &points).unwrap();
        let expected: Vec<f64> = points
            .iter()
            .map(|p| {
                let r = p.norm();
                spherical_bessel_j0(k * r) / spherical_bessel_j0(k * a)
            })
            .collect();
        for (i, e) in expected.iter().enumerate() {
            let got = p_bem[i].re;
            assert!(
                (got - e).abs() < 0.05,
                "point {i}: bem {got:.4} vs exact {e:.4}"
            );
        }
    }

    #[test]
    fn exterior_sphere_dirichlet() {
        // Prescribe the pulsating-sphere surface pressure, solve exterior
        // Dirichlet, evaluate at r = 2a.
        let a = 1.0;
        let k = 1.3;
        let v_a = Complex::new(1.0, 0.0);
        let mesh = BoundaryMesh::sphere(a, 10, 16).unwrap();
        let bem = HelmholtzBem::new(mesh, k).unwrap();
        let n = bem.mesh().len();
        let p_a = pulsating_sphere_exact(k, a, v_a, a);
        let g = vec![p_a; n];
        let mu = bem.solve_exterior_dirichlet(&g).unwrap();
        let points = [Vec3::new(0.0, 0.0, 2.0 * a)];
        let p_bem = bem.evaluate_double_layer(&mu, &points).unwrap();
        let p_exact = pulsating_sphere_exact(k, a, v_a, 2.0 * a);
        let rel = (p_bem[0] - p_exact).norm() / p_exact.norm();
        assert!(rel < 0.05, "BEM {p_bem:?} vs exact {p_exact:?} (rel {rel})");
    }

    #[test]
    fn radiated_power_of_pulsating_sphere() {
        // W = ½·Re(p·v*)·S evaluated with the EXACT surface quantities must
        // equal the analytic radiated power π ρc v² a² (ka)²/(1+(ka)²)·...
        // Actually W = ½·ρc·v²·(ka)²/(1+(ka)²)·4πa². Check the helper on
        // exact data.
        let a = 1.0;
        let k = 2.0;
        let rho_c = 413.0;
        let v = 1.0;
        let p_a = pulsating_sphere_exact(k, a, Complex::new(v, 0.0), a);
        let mesh = BoundaryMesh::sphere(a, 12, 20).unwrap();
        let bem = HelmholtzBem::new(mesh, k).unwrap();
        let n = bem.mesh().len();
        let p = vec![p_a; n];
        let q_val = p_a * (Complex::new(0.0, -k) + Complex::new(-1.0 / a, 0.0));
        // v_n = q/(−iωρ)... for power purposes use v = q/(iωρ) with ρc folded:
        // q = −iωρv → v = q/(−iωρ). With ρc = 413 and c = 343: ω = kc.
        let c = 343.0;
        let omega = k * c;
        let v_n = q_val / Complex::new(0.0, -omega) / Complex::new(rho_c / c, 0.0);
        let v_vec = alloc::vec![v_n; n];
        let w = bem.radiated_power(&p, &v_vec).unwrap();
        let w_exact = 0.5 * rho_c * v * v * (k * a * k * a) / (1.0 + (k * a) * (k * a))
            * 4.0
            * core::f64::consts::PI
            * a
            * a;
        let rel = (w - w_exact).abs() / w_exact;
        assert!(rel < 0.08, "w {w:.6} vs exact {w_exact:.6} (rel {rel})");
    }

    #[test]
    fn data_length_mismatch_is_rejected() {
        let mesh = BoundaryMesh::sphere(1.0, 4, 6).unwrap();
        let bem = HelmholtzBem::new(mesh, 1.0).unwrap();
        assert_eq!(
            bem.solve_exterior_dirichlet(&[Complex::new(1.0, 0.0); 3])
                .unwrap_err(),
            BemError::DataLengthMismatch
        );
    }

    #[test]
    fn bessel_helpers_match_series() {
        assert!(spherical_bessel_j0(0.0) - 1.0 < TOL);
        assert!((spherical_bessel_j0(1.0) - libm::sin(1.0)).abs() < TOL);
        assert!((spherical_bessel_j1(1e-10) - 1e-10 / 3.0).abs() < 1e-15);
    }
}
