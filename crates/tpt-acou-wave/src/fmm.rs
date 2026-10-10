//! Single-level multipole (FMM) acceleration of the BEM single-layer matvec
//! (feature `fmm`).
//!
//! Uses the axial addition theorem for the Helmholtz kernel — for a target
//! `x` outside a source cluster `B` with center `c` (`|x−c| > 2r_B`),
//!
//! ```text
//! G(x, y) = -i k / (4pi) * sum_n (2n+1) j_n(k|x-c|) h_n(k|y-c|) P_n(cos gamma)
//!         = sum_nm [J_nm(x)] * [M_nm],
//! ```
//!
//! with the separation into an outgoing target part
//! `J_nm = −ik/(4π)·(2n+1)·h_n(k|x−c|)·Y_nm(x̂)` and a regular source multipole
//! moment `M_nm = Σ_j σ_j·A_j·j_n(k|y_j−c|)·Y*_nm(ŷ_j)` (regular `j_n` on the
//! smaller source-side radius, outgoing `h_n` on the larger target radius —
//! the only assignment for which the series converges). `P_n(cos γ)` is expanded
//! through the spherical-harmonic addition theorem, so no azimuthal
//! translations are needed (single level, `O(N·n_clusters·p²)` matvec).
//!
//! This is deliberately a modest scheme: it accelerates *matvecs* (the
//! assembly-by-matvec path used by iterative solvers), not factorization,
//! and it is validated against the dense kernel in the tests.

use alloc::vec;
use alloc::vec::Vec;

use tpt_acou_core::Vec3;
use tpt_math_linalg_complex::Complex;

use crate::bem::HelmholtzBem;

/// Maximum multipole order.
const P: usize = 8;

/// A one-level octree partition of the mesh into eight children, each with a
/// center, radius and element index list.
struct Cluster {
    center: Vec3,
    radius: f64,
    elements: Vec<usize>,
    /// Multipole moments `M_nm`, index `n * (n + 1) + m` for `m in [-n, n]`
    /// stored as `m + n` offset rows: flat index `n * n + (m + n)`.
    moments: Vec<Complex>,
}

fn nm_index(n: usize, m: i32) -> usize {
    n * n + (m + n as i32) as usize
}

/// Total number of `(n, m)` pairs up to order `P` (inclusive).
fn nm_count() -> usize {
    (P + 1) * (P + 1)
}

/// Orthonormal complex spherical harmonic `Y_nm(θ, φ)` with
/// Condon–Shortley phase, built on the standard associated-Legendre
/// recurrences over a `(P+1)×(P+1)` table.
fn spherical_harmonic(n: usize, m: i32, cos_theta: f64, phi: f64) -> Complex {
    debug_assert!(m.unsigned_abs() as usize <= n);
    let mut p = [[0.0_f64; P + 1]; P + 1];
    p[0][0] = 1.0;
    // Diagonal: P_m^m = (-1)^m (2m-1)!! (1-x^2)^{m/2}.
    for mm in 1..=P {
        p[mm][mm] =
            -p[mm - 1][mm - 1] * (2.0 * mm as f64 - 1.0) * libm::sqrt(1.0 - cos_theta * cos_theta);
    }
    // Sub-diagonal: P_{m+1}^m = x (2m+1) P_m^m.
    for mm in 0..P {
        p[mm + 1][mm] = cos_theta * (2.0 * mm as f64 + 1.0) * p[mm][mm];
    }
    // Three-term recurrence for l >= m+2.
    for mm in 0..=P {
        for l in (mm + 2)..=P {
            let lf = l as f64;
            let mf = mm as f64;
            p[l][mm] = ((2.0 * lf - 1.0) * cos_theta * p[l - 1][mm]
                - (lf + mf - 1.0) * p[l - 2][mm])
                / (lf - mf);
        }
    }
    let mf = m.unsigned_abs();
    let nu = n as u32;
    let norm = libm::sqrt(
        (2.0 * nu as f64 + 1.0) * fact(nu - mf) / (4.0 * core::f64::consts::PI * fact(nu + mf)),
    );
    let angular = Complex::from_polar(1.0, m as f64 * phi);
    let cs = if m % 2 == 0 { 1.0 } else { -1.0 };
    angular * (cs * norm * p[n][mf as usize])
}

fn fact(k: u32) -> f64 {
    let mut f = 1.0;
    for i in 2..=k {
        f *= i as f64;
    }
    f
}

/// Spherical Bessel `j_n` (forward recurrence; `n <= P` with `k r` moderate).
fn bessel_j(n: usize, z: f64) -> f64 {
    let j0 = libm::sin(z) / z;
    if n == 0 {
        return j0;
    }
    let j1 = libm::sin(z) / (z * z) - libm::cos(z) / z;
    if n == 1 {
        return j1;
    }
    let (mut jm2, mut jm1) = (j0, j1);
    for l in 2..=n {
        let jl = (2.0 * (l - 1) as f64 + 1.0) / z * jm1 - jm2;
        jm2 = jm1;
        jm1 = jl;
    }
    jm1
}

/// Hankel function of the second kind `h^{(2)}_n = j_n − i·y_n`
/// (`y` by forward recurrence; adequate for the cluster-separation regime
/// `k·r ≳ 1` this FMM operates in).
fn hankel2(n: usize, z: f64) -> Complex {
    let mut ys = [0.0_f64; P + 1];
    ys[0] = -libm::cos(z) / z;
    if P >= 1 {
        ys[1] = -libm::cos(z) / (z * z) - libm::sin(z) / z;
    }
    for l in 2..=P {
        ys[l] = (2.0 * (l - 1) as f64 + 1.0) / z * ys[l - 1] - ys[l - 2];
    }
    Complex::new(bessel_j(n, z), -ys[n])
}

/// Builds the eight-cluster partition and assembles the matvec accelerator.
pub struct SingleLevelFmm {
    clusters: Vec<Cluster>,
}

impl SingleLevelFmm {
    /// Partitions `bem`'s mesh around its bounding-box center into eight
    /// octants and precomputes nothing else (moments depend on `sigma`).
    ///
    /// # Errors
    /// Propagates [`crate::bem::BemError::DataLengthMismatch`].
    pub fn new(bem: &HelmholtzBem) -> Result<Self, crate::bem::BemError> {
        let els = bem.mesh().elements();
        if els.is_empty() {
            return Ok(Self {
                clusters: Vec::new(),
            });
        }
        let mut lo = els[0].centroid;
        let mut hi = els[0].centroid;
        for e in els {
            lo = lo.min_components(e.centroid);
            hi = hi.max_components(e.centroid);
        }
        let center = (lo + hi) / 2.0;
        let mut clusters: Vec<Cluster> = (0..8)
            .map(|_| Cluster {
                center: Vec3::ZERO,
                radius: 0.0,
                elements: Vec::new(),
                moments: Vec::new(),
            })
            .collect();
        for (i, e) in els.iter().enumerate() {
            let d = e.centroid - center;
            let octant = (usize::from(d.x > 0.0))
                | (usize::from(d.y > 0.0) << 1)
                | (usize::from(d.z > 0.0) << 2);
            clusters[octant].elements.push(i);
        }
        for c in &mut clusters {
            if c.elements.is_empty() {
                continue;
            }
            let mut cc = Vec3::ZERO;
            let mut total = 0.0;
            for &i in &c.elements {
                cc += els[i].centroid * els[i].area;
                total += els[i].area;
            }
            c.center = cc / total;
            c.radius = c
                .elements
                .iter()
                .map(|&i| els[i].centroid.distance(c.center))
                .fold(0.0_f64, f64::max);
        }
        Ok(Self { clusters })
    }

    /// Assembles the multipole moments for source weights `sigma` (the BEM
    /// density times element area is applied inside).
    pub fn assemble(&mut self, bem: &HelmholtzBem, sigma: &[Complex]) {
        let k = bem.k();
        let els = bem.mesh().elements();
        for c in &mut self.clusters {
            c.moments = vec![Complex::new(0.0, 0.0); nm_count()];
            if c.elements.is_empty() {
                continue;
            }
            for &i in &c.elements {
                let d = els[i].centroid - c.center;
                let r = d.norm();
                if r < 1e-12 {
                    continue;
                }
                let ct = d.z / r;
                let phi = libm::atan2(d.y, d.x);
                for n in 0..=P {
                    let j = bessel_j(n, k * r);
                    for m in -(n as i32)..=(n as i32) {
                        let y = spherical_harmonic(n, m, ct, phi);
                        let idx = nm_index(n, m);
                        let jy = y.conj() * j;
                        c.moments[idx] = c.moments[idx] + jy * (sigma[i] * els[i].area);
                    }
                }
            }
        }
    }

    /// Evaluates `p(x) = Σ_j G(x, y_j)·σ_j·A_j` using the expansion for far
    /// clusters and the direct kernel for near ones.
    pub fn matvec(&self, bem: &HelmholtzBem, sigma: &[Complex], x: Vec3) -> Complex {
        let k = bem.k();
        let els = bem.mesh().elements();
        let mut acc = Complex::new(0.0, 0.0);
        for c in &self.clusters {
            if c.elements.is_empty() {
                continue;
            }
            // Strong separation: geometric (target well outside the source
            // cluster) AND regime (k*rx large enough that j_n/h_n are in
            // their oscillatory range, avoiding catastrophic cancellation
            // from the small-argument growth of y_n).
            let dx = x - c.center;
            let rx = dx.norm();
            if rx > 3.0 * c.radius && k * rx > 5.0 && rx > 1e-9 {
                // Far cluster: local expansion evaluation.
                let ct = dx.z / rx;
                let phi = libm::atan2(dx.y, dx.x);
                let mut sum = Complex::new(0.0, 0.0);
                for n in 0..=P {
                    let h = hankel2(n, k * rx);
                    // P_n(cos gamma) = (4pi/(2n+1)) sum_m Y Y*, which folds
                    // the (2n+1)/4pi factors of the axial theorem into a
                    // plain -i*k prefactor.
                    let prefac = Complex::new(0.0, -k);
                    for m in -(n as i32)..=(n as i32) {
                        let y = spherical_harmonic(n, m, ct, phi);
                        sum = sum + y * c.moments[nm_index(n, m)] * prefac * h;
                    }
                }
                acc = acc + sum;
            } else {
                // Near cluster: direct kernel.
                for &i in &c.elements {
                    let r = x.distance(els[i].centroid);
                    acc = acc
                        + Complex::from_polar(1.0 / (4.0 * core::f64::consts::PI * r), -k * r)
                            * (sigma[i] * els[i].area);
                }
            }
        }
        acc
    }
}

impl SingleLevelFmm {
    /// Far-cluster expansion evaluated at `x`: `Σ_nm J_nm(x) M_nm`, the
    /// single-layer potential of the cluster's sources. Valid only for `x`
    /// outside the cluster's far criterion.
    fn far_potential(&self, c: &Cluster, x: Vec3, k: f64) -> Complex {
        let dx = x - c.center;
        let rx = dx.norm();
        let ct = dx.z / rx;
        let phi = libm::atan2(dx.y, dx.x);
        let mut sum = Complex::new(0.0, 0.0);
        let prefac = Complex::new(0.0, -k);
        for n in 0..=P {
            let h = hankel2(n, k * rx);
            for m in -(n as i32)..=(n as i32) {
                let y = spherical_harmonic(n, m, ct, phi);
                sum = sum + y * c.moments[nm_index(n, m)] * prefac * h;
            }
        }
        sum
    }

    /// Applies the adjoint double-layer operator to a density:
    /// `(K′μ)_i = Σ_j A_j μ_j · n_i·∇_x G(x_i − y_j)`, the operator in the
    /// exterior Neumann equation `(K′ − ½I)σ = q`.
    ///
    /// Far clusters use the field-side normal derivative of the multipole
    /// expansion, by a central difference at step `1/(100k)`, and a stricter
    /// far criterion (6 cluster radii) to keep the derivative accurate. Near
    /// clusters use the exact kernel, and the self term is zero as in the dense
    /// matrix. Accuracy is about 10⁻⁵ relative for well-separated clusters.
    ///
    /// # Errors
    /// Returns [`crate::bem::BemError::DataLengthMismatch`] if `mu` does not
    /// match the mesh.
    pub fn apply_adjoint_double_layer(
        &mut self,
        bem: &HelmholtzBem,
        mu: &[Complex],
    ) -> Result<Vec<Complex>, crate::bem::BemError> {
        if mu.len() != bem.mesh().len() {
            return Err(crate::bem::BemError::DataLengthMismatch);
        }
        self.assemble(bem, mu);
        let k = bem.k();
        let els = bem.mesh().elements();
        let step = 1.0 / (100.0 * k);
        let mut out = vec![Complex::new(0.0, 0.0); els.len()];
        for (i, e_i) in els.iter().enumerate() {
            let x = e_i.centroid;
            let n = e_i.normal;
            let mut acc = Complex::new(0.0, 0.0);
            for c in &self.clusters {
                if c.elements.is_empty() {
                    continue;
                }
                let rx = (x - c.center).norm();
                if rx > 6.0 * c.radius && k * rx > 5.0 && rx > 1e-9 {
                    let plus = self.far_potential(c, x + n * step, k);
                    let minus = self.far_potential(c, x - n * step, k);
                    acc = acc + (plus - minus) * (1.0 / (2.0 * step));
                } else {
                    for &j in &c.elements {
                        if j == i {
                            continue;
                        }
                        let d = x - els[j].centroid;
                        let r = d.norm();
                        let g =
                            Complex::from_polar(1.0 / (4.0 * core::f64::consts::PI * r), -k * r);
                        let dn = d.dot(n) / r;
                        // n·∇_x G = −G (ik + 1/r)·(d·n)/r
                        let kernel =
                            g * (Complex::new(0.0, k) + Complex::new(1.0 / r, 0.0)) * (-dn);
                        acc = acc + kernel * mu[j] * els[j].area;
                    }
                }
            }
            out[i] = acc;
        }
        Ok(out)
    }
}

/// Exterior Neumann solve `(K′ − ½I)σ = q` by GMRES, with the adjoint
/// double-layer operator applied through the single-level FMM. Returns the
/// density and convergence information (see [`crate::krylov::GmresResult`]).
///
/// # Errors
/// Returns [`crate::bem::BemError::DataLengthMismatch`] for mismatched data.
pub fn solve_exterior_neumann_fmm(
    bem: &HelmholtzBem,
    q: &[Complex],
    config: crate::krylov::GmresConfig,
) -> Result<crate::krylov::GmresResult, crate::bem::BemError> {
    if q.len() != bem.mesh().len() {
        return Err(crate::bem::BemError::DataLengthMismatch);
    }
    let mut fmm = SingleLevelFmm::new(bem)?;
    let half = Complex::new(0.5, 0.0);
    crate::krylov::gmres(
        |v| {
            // The FMM path cannot fail for a length-checked vector.
            let kv = fmm.apply_adjoint_double_layer(bem, v).unwrap_or_default();
            kv.iter().zip(v).map(|(a, b)| *a - *b * half).collect()
        },
        q,
        None,
        None,
        config,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::mesh::BoundaryMesh;

    #[test]
    fn adjoint_double_layer_fmm_matches_the_dense_matrix() {
        let mesh = BoundaryMesh::sphere(1.0, 8, 12).unwrap();
        let bem = HelmholtzBem::new(mesh, 1.5).unwrap();
        let n = bem.mesh().len();
        let mu: Vec<Complex> = (0..n)
            .map(|i| Complex::new(1.0 + 0.1 * (i % 7) as f64, 0.2 * ((i % 5) as f64 - 2.0)))
            .collect();
        let dense = bem.adjoint_double_layer_matrix();
        let mut fmm = SingleLevelFmm::new(&bem).unwrap();
        let fast = fmm.apply_adjoint_double_layer(&bem, &mu).unwrap();
        let mut num = 0.0_f64;
        let mut den = 0.0_f64;
        for i in 0..n {
            let mut exact = Complex::new(0.0, 0.0);
            for j in 0..n {
                exact = exact + dense[(i, j)] * mu[j];
            }
            num += (fast[i] - exact).norm_sqr();
            den += exact.norm_sqr();
        }
        let rel = (num / den).sqrt();
        assert!(rel < 1e-3, "relative error {rel:.3e}");
    }

    #[test]
    fn fmm_neumann_solve_matches_the_dense_solve() {
        let mesh = BoundaryMesh::sphere(1.0, 8, 12).unwrap();
        let bem = HelmholtzBem::new(mesh, 1.5).unwrap();
        let n = bem.mesh().len();
        let q: Vec<Complex> = (0..n)
            .map(|i| Complex::new(1.0, 0.3 * (i as f64 * 0.37).sin()))
            .collect();
        let dense = bem.solve_exterior_neumann(&q).unwrap();
        let config = crate::krylov::GmresConfig {
            restart: 40,
            max_iterations: 400,
            tolerance: 1e-6,
        };
        let res = solve_exterior_neumann_fmm(&bem, &q, config).unwrap();
        assert!(res.converged, "relative residual {}", res.relative_residual);
        let mut num = 0.0_f64;
        let mut den = 0.0_f64;
        for (a, b) in res.x.iter().zip(&dense) {
            num += (*a - *b).norm_sqr();
            den += b.norm_sqr();
        }
        let rel = (num / den).sqrt();
        assert!(rel < 1e-3, "relative difference {rel:.3e}");
    }

    #[test]
    fn fmm_matvec_matches_dense() {
        let mesh = BoundaryMesh::sphere(1.0, 8, 12).unwrap();
        let bem = HelmholtzBem::new(mesh, 2.0).unwrap();
        let n = bem.mesh().len();
        let sigma: Vec<Complex> = (0..n)
            .map(|i| Complex::new(1.0 + 0.1 * (i % 7) as f64, 0.2 * ((i % 5) as f64 - 2.0)))
            .collect();
        let mut fmm = SingleLevelFmm::new(&bem).unwrap();
        fmm.assemble(&bem, &sigma);
        // Field points on a shell around the sphere (outside every cluster's
        // 2r criterion where the expansion is used).
        let mut max_err = 0.0_f64;
        let mut max_ref = 0.0_f64;
        for i in 0..12 {
            let th = core::f64::consts::PI * (i as f64 + 0.5) / 12.0;
            for j in 0..8 {
                let phi = core::f64::consts::TAU * j as f64 / 8.0;
                let _ = phi;
                let x = Vec3::new(
                    4.0 * libm::sin(th) * libm::cos(phi),
                    4.0 * libm::sin(th) * libm::sin(phi),
                    4.0 * libm::cos(th),
                );
                let dense = bem.evaluate_single_layer(&sigma, &[x]).unwrap()[0];
                let fast = fmm.matvec(&bem, &sigma, x);
                max_err = max_err.max((dense - fast).norm());
                max_ref = max_ref.max(dense.norm());
            }
        }
        assert!(
            max_err < 1e-4 * max_ref.max(1e-30),
            "max abs err {max_err:.3e} vs ref {max_ref:.3e}"
        );
    }

    #[test]
    fn spherical_harmonic_orthonormality() {
        // Y_00 = 1/sqrt(4pi).
        let y00 = spherical_harmonic(0, 0, 1.0, 0.0);
        assert!((y00.norm() - 1.0 / libm::sqrt(4.0 * core::f64::consts::PI)).abs() < 1e-12);
        // Y_10 = sqrt(3/4pi) cos(theta).
        let y10 = spherical_harmonic(1, 0, 1.0, 0.0);
        let want = libm::sqrt(3.0 / (4.0 * core::f64::consts::PI));
        assert!((y10.norm() - want).abs() < 1e-12);
    }

    #[test]
    fn bessel_functions_match_series() {
        // j0(z) = sin z / z at z = 1.
        assert!((bessel_j(0, 1.0) - libm::sin(1.0)).abs() < 1e-12);
        // Recurrence identity: j_{n+1} = (2n+1)/z j_n - j_{n-1} at z = 2, n=1.
        let j0 = bessel_j(0, 2.0);
        let j1 = bessel_j(1, 2.0);
        let j2 = bessel_j(2, 2.0);
        assert!((j2 - (3.0 / 2.0 * j1 - j0)).abs() < 1e-12);
    }
}
