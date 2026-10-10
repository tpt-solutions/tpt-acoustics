//! Restarted, right-preconditioned GMRES for complex linear systems.
//!
//! The operator is supplied as a closure `apply(x) = A·x`, so the same solver
//! runs on a dense BEM matrix or on a fast-multipole matvec without ever
//! forming the matrix. A Jacobi (diagonal) preconditioner can be supplied as
//! the inverse diagonal; the iteration is right-preconditioned, so the
//! residual it monitors is the true residual of `A x = b`.
//!
//! Givens rotations are complex-valued: with `c` real and `s` complex each
//! rotation is unitary, which keeps the least-squares residual exact.

use alloc::vec;
use alloc::vec::Vec;

use tpt_math_linalg_complex::Complex;

use crate::bem::BemError;

/// Outcome of a GMRES solve.
#[derive(Debug, Clone, PartialEq)]
pub struct GmresResult {
    /// Approximate solution.
    pub x: Vec<Complex>,
    /// Number of Arnoldi steps taken across all restart cycles.
    pub iterations: usize,
    /// Final relative residual `‖b − A x‖ / ‖b‖`.
    pub relative_residual: f64,
    /// Whether the relative residual reached the requested tolerance.
    pub converged: bool,
}

/// Tuning for [`gmres`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GmresConfig {
    /// Krylov subspace size before a restart.
    pub restart: usize,
    /// Maximum total Arnoldi steps.
    pub max_iterations: usize,
    /// Target relative residual.
    pub tolerance: f64,
}

impl Default for GmresConfig {
    fn default() -> Self {
        Self {
            restart: 30,
            max_iterations: 300,
            tolerance: 1e-10,
        }
    }
}

fn dot(a: &[Complex], b: &[Complex]) -> Complex {
    // Conjugate-linear in the first argument: ⟨a, b⟩ = Σ conj(aᵢ) bᵢ.
    a.iter()
        .zip(b)
        .fold(Complex::new(0.0, 0.0), |acc, (x, y)| acc + x.conj() * *y)
}

fn norm(a: &[Complex]) -> f64 {
    libm::sqrt(a.iter().map(|x| x.norm_sqr()).sum::<f64>())
}

/// Complex Givens rotation `(c, s)` that zeroes `b` in `[a; b]`:
/// `c = |a| / r`, `s = (a/|a|)·conj(b) / r`, with `r = √(|a|² + |b|²)`.
fn givens(a: Complex, b: Complex) -> (f64, Complex) {
    let (na, nb) = (a.norm(), b.norm());
    let r = libm::sqrt(na * na + nb * nb);
    if r == 0.0 {
        return (1.0, Complex::new(0.0, 0.0));
    }
    if na == 0.0 {
        return (0.0, b.conj() * (1.0 / nb));
    }
    let phase = a * (1.0 / na);
    (na / r, phase * b.conj() * (1.0 / r))
}

/// Applies the rotation `(c, s)` to the pair `(x, y)`.
fn rotate(c: f64, s: Complex, x: Complex, y: Complex) -> (Complex, Complex) {
    (x * c + s * y, s.conj() * (-1.0) * x + y * c)
}

/// Solves `A x = b` by restarted GMRES.
///
/// `apply` computes `A·v` for a vector `v`. `inverse_diagonal`, if given,
/// is the element-wise inverse of the diagonal of `A` (a Jacobi right
/// preconditioner). `x0` is the initial guess; zero is used if it is `None`.
///
/// Non-convergence is not an error: the result reports `converged = false`
/// together with the residual reached.
///
/// # Errors
/// [`BemError::DataLengthMismatch`] if `x0` or `inverse_diagonal` does not
/// have the length of `b`.
pub fn gmres<F>(
    mut apply: F,
    b: &[Complex],
    x0: Option<&[Complex]>,
    inverse_diagonal: Option<&[Complex]>,
    config: GmresConfig,
) -> Result<GmresResult, BemError>
where
    F: FnMut(&[Complex]) -> Vec<Complex>,
{
    let n = b.len();
    if x0.is_some_and(|v| v.len() != n) || inverse_diagonal.is_some_and(|v| v.len() != n) {
        return Err(BemError::DataLengthMismatch);
    }
    let precondition = |v: &[Complex]| -> Vec<Complex> {
        match inverse_diagonal {
            Some(d) => v.iter().zip(d).map(|(a, m)| *a * *m).collect(),
            None => v.to_vec(),
        }
    };

    let mut x: Vec<Complex> = match x0 {
        Some(v) => v.to_vec(),
        None => vec![Complex::new(0.0, 0.0); n],
    };
    let b_norm = norm(b);
    if b_norm == 0.0 {
        return Ok(GmresResult {
            x: vec![Complex::new(0.0, 0.0); n],
            iterations: 0,
            relative_residual: 0.0,
            converged: true,
        });
    }

    let restart = config.restart.max(1);
    let mut iterations = 0;
    let mut relative_residual;

    loop {
        // Residual at the start of this cycle.
        let ax = apply(&x);
        let r: Vec<Complex> = b.iter().zip(&ax).map(|(bi, ai)| *bi - *ai).collect();
        let beta = norm(&r);
        relative_residual = beta / b_norm;
        if relative_residual <= config.tolerance || iterations >= config.max_iterations {
            return Ok(GmresResult {
                x,
                iterations,
                relative_residual,
                converged: relative_residual <= config.tolerance,
            });
        }

        let mut basis: Vec<Vec<Complex>> = Vec::with_capacity(restart + 1);
        basis.push(r.iter().map(|v| *v * (1.0 / beta)).collect());
        // Hessenberg matrix H[i][j], i ≤ j + 1, and the rotated RHS g.
        let mut h = vec![vec![Complex::new(0.0, 0.0); restart]; restart + 1];
        let mut cs: Vec<(f64, Complex)> = Vec::with_capacity(restart);
        let mut g = vec![Complex::new(0.0, 0.0); restart + 1];
        g[0] = Complex::new(beta, 0.0);

        let mut steps = 0;
        for j in 0..restart {
            // Right preconditioning: w = A M⁻¹ v_j.
            let z = precondition(&basis[j]);
            let mut w = apply(&z);
            for i in 0..=j {
                let hij = dot(&basis[i], &w);
                h[i][j] = hij;
                for (wk, vk) in w.iter_mut().zip(&basis[i]) {
                    *wk = *wk - hij * *vk;
                }
            }
            let h_next = norm(&w);
            h[j + 1][j] = Complex::new(h_next, 0.0);
            if h_next > 0.0 {
                basis.push(w.iter().map(|v| *v * (1.0 / h_next)).collect());
            } else {
                // Happy breakdown: the Krylov space is invariant.
                basis.push(vec![Complex::new(0.0, 0.0); n]);
            }

            // Apply the previous rotations to the new column.
            for i in 0..j {
                let (c, s) = cs[i];
                let (u, v) = rotate(c, s, h[i][j], h[i + 1][j]);
                h[i][j] = u;
                h[i + 1][j] = v;
            }
            // New rotation annihilating h[j+1][j].
            let (c, s) = givens(h[j][j], h[j + 1][j]);
            cs.push((c, s));
            let (u, v) = rotate(c, s, h[j][j], h[j + 1][j]);
            h[j][j] = u;
            h[j + 1][j] = v;
            let (gu, gv) = rotate(c, s, g[j], g[j + 1]);
            g[j] = gu;
            g[j + 1] = gv;

            steps = j + 1;
            iterations += 1;
            relative_residual = g[j + 1].norm() / b_norm;
            if relative_residual <= config.tolerance
                || iterations >= config.max_iterations
                || h_next == 0.0
            {
                break;
            }
        }

        // Back-substitute the triangular system H y = g.
        let mut y = vec![Complex::new(0.0, 0.0); steps];
        for i in (0..steps).rev() {
            let mut acc = g[i];
            for k in (i + 1)..steps {
                acc = acc - h[i][k] * y[k];
            }
            y[i] = acc * h[i][i].inv().unwrap_or(Complex::new(0.0, 0.0));
        }

        // x ← x + M⁻¹ (V y).
        let mut update = vec![Complex::new(0.0, 0.0); n];
        for (i, yi) in y.iter().enumerate() {
            for (uk, vk) in update.iter_mut().zip(&basis[i]) {
                *uk = *uk + *yi * *vk;
            }
        }
        let correction = precondition(&update);
        for (xk, ck) in x.iter_mut().zip(&correction) {
            *xk = *xk + *ck;
        }

        if relative_residual <= config.tolerance || iterations >= config.max_iterations {
            // Re-check the true residual on the way out.
            let ax = apply(&x);
            let res: Vec<Complex> = b.iter().zip(&ax).map(|(bi, ai)| *bi - *ai).collect();
            let true_rel = norm(&res) / b_norm;
            return Ok(GmresResult {
                x,
                iterations,
                relative_residual: true_rel,
                converged: true_rel <= config.tolerance,
            });
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(re: f64, im: f64) -> Complex {
        Complex::new(re, im)
    }

    /// Diagonally dominant complex test matrix, n × n, row-major.
    fn test_matrix(n: usize) -> Vec<Vec<Complex>> {
        (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        if i == j {
                            c(4.0 + i as f64 * 0.1, 0.5)
                        } else {
                            c(
                                0.3 / (1.0 + (i as f64 - j as f64).abs()),
                                -0.2 * ((i + j) as f64).sin(),
                            )
                        }
                    })
                    .collect()
            })
            .collect()
    }

    fn matvec(a: &[Vec<Complex>], v: &[Complex]) -> Vec<Complex> {
        a.iter()
            .map(|row| {
                row.iter()
                    .zip(v)
                    .fold(c(0.0, 0.0), |acc, (x, y)| acc + *x * *y)
            })
            .collect()
    }

    #[test]
    fn solves_a_diagonally_dominant_system() {
        let n = 12;
        let a = test_matrix(n);
        let x_true: Vec<Complex> = (0..n).map(|i| c(1.0 + i as f64, -0.5 * i as f64)).collect();
        let b = matvec(&a, &x_true);
        let res = gmres(|v| matvec(&a, v), &b, None, None, GmresConfig::default()).unwrap();
        assert!(res.converged, "residual {}", res.relative_residual);
        for (xi, ti) in res.x.iter().zip(&x_true) {
            assert!((*xi - *ti).norm() < 1e-8, "{xi:?} vs {ti:?}");
        }
    }

    #[test]
    fn restarts_do_not_break_convergence() {
        // A restart length smaller than the problem forces several cycles.
        let n = 20;
        let a = test_matrix(n);
        let b: Vec<Complex> = (0..n).map(|i| c((i as f64).cos(), 0.1)).collect();
        let config = GmresConfig {
            restart: 4,
            max_iterations: 400,
            tolerance: 1e-9,
        };
        let res = gmres(|v| matvec(&a, v), &b, None, None, config).unwrap();
        assert!(res.converged);
        assert!(res.iterations > 4, "expected more than one cycle");
        let back = matvec(&a, &res.x);
        let err = norm(
            &back
                .iter()
                .zip(&b)
                .map(|(p, q)| *p - *q)
                .collect::<Vec<_>>(),
        );
        assert!(err / norm(&b) < 1e-8);
    }

    #[test]
    fn jacobi_preconditioner_gives_the_same_solution() {
        let n = 10;
        let a = test_matrix(n);
        let inv_diag: Vec<Complex> = (0..n).map(|i| a[i][i].inv().unwrap()).collect();
        let b: Vec<Complex> = (0..n).map(|i| c(1.0, i as f64 * 0.2)).collect();
        let plain = gmres(|v| matvec(&a, v), &b, None, None, GmresConfig::default()).unwrap();
        let precond = gmres(
            |v| matvec(&a, v),
            &b,
            None,
            Some(&inv_diag),
            GmresConfig::default(),
        )
        .unwrap();
        assert!(plain.converged && precond.converged);
        for (p, q) in plain.x.iter().zip(&precond.x) {
            assert!((*p - *q).norm() < 1e-8);
        }
    }

    #[test]
    fn zero_rhs_returns_zero() {
        let a = test_matrix(3);
        let b = vec![c(0.0, 0.0); 3];
        let res = gmres(|v| matvec(&a, v), &b, None, None, GmresConfig::default()).unwrap();
        assert!(res.converged);
        assert!(res.x.iter().all(|v| v.norm() == 0.0));
    }

    #[test]
    fn reports_non_convergence_without_error() {
        let a = test_matrix(8);
        let b: Vec<Complex> = (0..8).map(|i| c(i as f64, 1.0)).collect();
        let config = GmresConfig {
            restart: 2,
            max_iterations: 2,
            tolerance: 1e-14,
        };
        let res = gmres(|v| matvec(&a, v), &b, None, None, config).unwrap();
        assert!(!res.converged);
        assert!(res.relative_residual > 1e-14);
    }

    #[test]
    fn rejects_mismatched_lengths() {
        let a = test_matrix(3);
        let b = vec![c(1.0, 0.0); 3];
        let short = vec![c(0.0, 0.0); 2];
        assert_eq!(
            gmres(
                |v| matvec(&a, v),
                &b,
                Some(&short),
                None,
                GmresConfig::default()
            ),
            Err(BemError::DataLengthMismatch)
        );
    }

    #[test]
    fn matches_the_dense_bem_exterior_solve() {
        // Same system as `HelmholtzBem::solve_exterior_dirichlet`, solved
        // iteratively with the operator applied as a matvec.
        let mesh = crate::mesh::BoundaryMesh::sphere(1.0, 6, 8).unwrap();
        let bem = crate::bem::HelmholtzBem::new(mesh, 1.3).unwrap();
        let n = bem.mesh().len();
        let d = bem.double_layer_matrix();
        let a: Vec<Vec<Complex>> = (0..n)
            .map(|i| {
                (0..n)
                    .map(|j| {
                        if i == j {
                            d[(i, j)] + c(-0.5, 0.0)
                        } else {
                            d[(i, j)]
                        }
                    })
                    .collect()
            })
            .collect();
        let g: Vec<Complex> = (0..n).map(|i| c(1.0, 0.2 * (i as f64).sin())).collect();

        let dense = bem.solve_exterior_dirichlet(&g).unwrap();
        let config = GmresConfig {
            restart: n,
            max_iterations: 4 * n,
            tolerance: 1e-11,
        };
        let res = gmres(|v| matvec(&a, v), &g, None, None, config).unwrap();
        assert!(res.converged, "residual {}", res.relative_residual);
        let num: f64 = res
            .x
            .iter()
            .zip(&dense)
            .map(|(p, q)| (*p - *q).norm_sqr())
            .sum::<f64>()
            .sqrt();
        let den: f64 = dense.iter().map(|q| q.norm_sqr()).sum::<f64>().sqrt();
        assert!(num / den < 1e-7, "relative difference {}", num / den);
    }

    #[test]
    fn givens_rotation_zeroes_the_second_entry() {
        let (a, b) = (c(3.0, -1.0), c(-2.0, 4.0));
        let (cc, s) = givens(a, b);
        let (u, v) = rotate(cc, s, a, b);
        assert!(v.norm() < 1e-12, "second entry {v:?}");
        assert!((u.norm() - (a.norm_sqr() + b.norm_sqr()).sqrt()).abs() < 1e-12);
    }
}
