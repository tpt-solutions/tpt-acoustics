//! Interior acoustic cavity FEM (feature `fem`).
//!
//! Structured hexahedral meshes of rectangular cavities discretized with
//! 8-node trilinear pressure elements and 2×2×2 Gauss quadrature; assembly
//! goes through `tpt-fem-sparse`'s duplicate-summing [`tpt_fem_sparse::Coo`]
//! accumulator (the tpt-fem integration named in spec.txt §4). Cavity modes
//! are extracted with shift-inverse iteration on `(K − σM)φ = Mφ`, i.e. the
//! discrete Helmholtz eigenproblem `Kφ = ω²Mφ` for the rigid-walled cavity.
//!
//! Rigid walls give a zero-eigenvalue constant-pressure mode, so callers pick
//! a shift `σ = s²` (s in rad/s) below the mode of interest; the iteration
//! converges to the eigenpair nearest the shift.

use alloc::vec::Vec;

use tpt_acou_core::error::AcousticError;
use tpt_fem_sparse::{solve, Coo};

/// One assembled matrix entry (duplicate-summed on use).
type Triple = (usize, usize, f64);

/// A rectangular acoustic cavity discretized with trilinear hexahedra.
#[derive(Debug, Clone)]
pub struct CavityFem {
    dims: [f64; 3],
    cells: [usize; 3],
    k: Vec<Triple>,
    m: Vec<Triple>,
    n_dofs: usize,
}

/// Number of Gauss points per axis (2×2×2 exact for trilinear bricks).
const GAUSS: [f64; 2] = [-0.577_350_269_189_625_8, 0.577_350_269_189_625_8];
/// Gauss weight per axis.
const WEIGHT: f64 = 1.0;

impl CavityFem {
    /// Builds a box cavity `lx × ly × lz` with `cells` hexahedra per axis.
    ///
    /// # Errors
    /// Returns [`AcousticError::NonPositive`] for non-positive extents or
    /// zero cell counts.
    pub fn box_cavity(dims: [f64; 3], cells: [usize; 3]) -> Result<Self, AcousticError> {
        for d in dims {
            if d <= 0.0 {
                return Err(AcousticError::NonPositive("cavity extent"));
            }
        }
        for c in cells {
            if c == 0 {
                return Err(AcousticError::NonPositive("cavity cell count"));
            }
        }
        let [lx, ly, lz] = dims;
        let [nx, ny, nz] = cells;
        let n_dofs = (nx + 1) * (ny + 1) * (nz + 1);
        let node = |i: usize, j: usize, k: usize| i + (nx + 1) * (j + (ny + 1) * k);

        let hx = lx / nx as f64;
        let hy = ly / ny as f64;
        let hz = lz / nz as f64;
        // Jacobian determinant of the axis-aligned brick mapping.
        let det_j = (hx / 2.0) * (hy / 2.0) * (hz / 2.0);

        // Corner sign table for the trilinear hex.
        let signs: [(f64, f64, f64); 8] = [
            (-1.0, -1.0, -1.0),
            (1.0, -1.0, -1.0),
            (1.0, 1.0, -1.0),
            (-1.0, 1.0, -1.0),
            (-1.0, -1.0, 1.0),
            (1.0, -1.0, 1.0),
            (1.0, 1.0, 1.0),
            (-1.0, 1.0, 1.0),
        ];

        let mut k_coo: Vec<Triple> = Vec::new();
        let mut m_coo: Vec<Triple> = Vec::new();

        for kz in 0..nz {
            for jy in 0..ny {
                for ix in 0..nx {
                    let conn = [
                        node(ix, jy, kz),
                        node(ix + 1, jy, kz),
                        node(ix + 1, jy + 1, kz),
                        node(ix, jy + 1, kz),
                        node(ix, jy, kz + 1),
                        node(ix + 1, jy, kz + 1),
                        node(ix + 1, jy + 1, kz + 1),
                        node(ix, jy + 1, kz + 1),
                    ];
                    let mut ke = [[0.0_f64; 8]; 8];
                    let mut me = [[0.0_f64; 8]; 8];
                    for &xi in &GAUSS {
                        for &eta in &GAUSS {
                            for &zeta in &GAUSS {
                                for a in 0..8 {
                                    let (sa, ta, ua) = signs[a];
                                    let na = 0.125
                                        * (1.0 + sa * xi)
                                        * (1.0 + ta * eta)
                                        * (1.0 + ua * zeta);
                                    // dN/d(x,y,z) via the chain rule with the
                                    // diagonal Jacobian.
                                    let dnx = 0.125
                                        * sa
                                        * (1.0 / (hx / 2.0))
                                        * (1.0 + ta * eta)
                                        * (1.0 + ua * zeta);
                                    let dny = 0.125
                                        * (1.0 + sa * xi)
                                        * (ta * (1.0 / (hy / 2.0)))
                                        * (1.0 + ua * zeta);
                                    let dnz = 0.125
                                        * (1.0 + sa * xi)
                                        * (1.0 + ta * eta)
                                        * (ua * (1.0 / (hz / 2.0)));
                                    for b in 0..8 {
                                        let (tb, ub, wb) = signs[b];
                                        let nb = 0.125
                                            * (1.0 + tb * xi)
                                            * (1.0 + ub * eta)
                                            * (1.0 + wb * zeta);
                                        let dnbx = 0.125
                                            * tb
                                            * (1.0 / (hx / 2.0))
                                            * (1.0 + ub * eta)
                                            * (1.0 + wb * zeta);
                                        let dnby = 0.125
                                            * (1.0 + tb * xi)
                                            * (ub * (1.0 / (hy / 2.0)))
                                            * (1.0 + wb * zeta);
                                        let dnbz = 0.125
                                            * (1.0 + tb * xi)
                                            * (1.0 + ub * eta)
                                            * (wb * (1.0 / (hz / 2.0)));
                                        ke[a][b] +=
                                            (dnx * dnbx + dny * dnby + dnz * dnbz) * det_j * WEIGHT;
                                        me[a][b] += na * nb * det_j * WEIGHT;
                                    }
                                }
                            }
                        }
                    }
                    for a in 0..8 {
                        for b in 0..8 {
                            k_coo.push((conn[a], conn[b], ke[a][b]));
                            m_coo.push((conn[a], conn[b], me[a][b]));
                        }
                    }
                }
            }
        }

        Ok(Self {
            dims,
            cells,
            k: k_coo,
            m: m_coo,
            n_dofs,
        })
    }

    /// Cavity dimensions.
    #[must_use]
    pub const fn dims(&self) -> [f64; 3] {
        self.dims
    }

    /// Cells per axis.
    #[must_use]
    pub const fn cells(&self) -> [usize; 3] {
        self.cells
    }

    /// Number of pressure DOFs.
    #[must_use]
    pub const fn n_dofs(&self) -> usize {
        self.n_dofs
    }

    /// Stiffness matrix, assembled into a duplicate-summing COO.
    #[must_use]
    pub fn stiffness(&self) -> Coo {
        let mut coo = Coo::with_capacity(self.k.len());
        for &(r, c, v) in &self.k {
            coo.push(r, c, v);
        }
        coo
    }

    /// Mass matrix.
    #[must_use]
    pub fn mass(&self) -> Coo {
        let mut coo = Coo::with_capacity(self.m.len());
        for &(r, c, v) in &self.m {
            coo.push(r, c, v);
        }
        coo
    }

    /// Eigenpair `(k², φ)` of the discrete Helmholtz operator
    /// `Kφ = k²Mφ` nearest the shift `sigma = s²` (s in rad/m), by
    /// shift-inverse iteration `(K − σM)φ_{n+1} = Mφ_n`.
    /// Returns `(k_rad_per_m, shape)`; the physical frequency is
    /// `f = c·k / 2π`. (The assembly carries no `c²`, so the generalized
    /// eigenvalue is a squared *wavenumber*.)
    ///
    /// # Errors
    /// Propagates solver errors from `tpt-fem-sparse` (singular shifted
    /// system, i.e. `sigma` exactly at an eigenvalue).
    pub fn mode_nearest(
        &self,
        sigma: f64,
        iterations: usize,
    ) -> Result<(f64, Vec<f64>), tpt_fem_sparse::SparseError> {
        let n = self.n_dofs;
        // Shifted operator (K − σM) as a fresh Coo.
        let mut shifted = Coo::with_capacity(self.k.len() + self.m.len());
        for &(r, c, v) in &self.k {
            shifted.push(r, c, v);
        }
        for &(r, c, v) in &self.m {
            shifted.push(r, c, -sigma * v);
        }

        // Start from a smooth ramp (non-orthogonal to the fundamental modes).
        let mut phi: Vec<f64> = (0..n)
            .map(|i| {
                let k = i as f64;
                libm::sin(0.7 * k + 0.3) + 0.5 * libm::cos(0.11 * k)
            })
            .collect();
        let mut rhs = coo_matvec(&self.m, &phi);

        for _ in 0..iterations {
            phi = solve(&shifted, &rhs)?;
            // Normalize (max-abs) and recompute the mass-weighted vector.
            let norm = phi.iter().fold(0.0_f64, |a, &v| a.max(v.abs()));
            for v in &mut phi {
                *v /= norm;
            }
            rhs = coo_matvec(&self.m, &phi);
        }
        // Rayleigh quotient k² = φᵀKφ / φᵀMφ.
        let kphi = coo_matvec(&self.k, &phi);
        let mphi = &rhs;
        let num: f64 = phi.iter().zip(&kphi).map(|(a, b)| a * b).sum();
        let den: f64 = phi.iter().zip(mphi).map(|(a, b)| a * b).sum();
        Ok((num / den, phi))
    }
}

/// Dense matvec over assembled triples (duplicate-summed).
fn coo_matvec(triples: &[Triple], x: &[f64]) -> Vec<f64> {
    let mut y = alloc::vec![0.0; x.len()];
    for &(r, c, v) in triples {
        y[r] += v * x[c];
    }
    y
}

/// Analytic rigid-cavity modal frequency `f = c/2·sqrt((nx/Lx)² + (ny/Ly)² +
/// (nz/Lz)²)` for mode indices `(nx, ny, nz) ≥ 0` (not all zero).
/// The wavenumber is `k = 2πf/c = π·sqrt((nx/Lx)² + (ny/Ly)² + (nz/Lz)²)`.
#[must_use]
pub fn analytic_mode_frequency(dims: [f64; 3], c: f64, mode: [u32; 3]) -> f64 {
    let [lx, ly, lz] = dims;
    let [mx, my, mz] = mode;
    let a = mx as f64 / lx;
    let b = my as f64 / ly;
    let d = mz as f64 / lz;
    0.5 * c * libm::sqrt(a * a + b * b + d * d)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fundamental_mode_of_box_cavity() {
        // 1 × 0.8 × 0.6 m cavity, (1,0,0) fundamental at c/(2·Lx).
        let dims = [1.0, 0.8, 0.6];
        let c = 343.0;
        let fem = CavityFem::box_cavity(dims, [8, 7, 5]).unwrap();
        // Exact wavenumber of the (1,0,0) mode: k = π/Lx.
        let k_exact = core::f64::consts::PI / dims[0];
        // Shift just below the target mode so the zero (constant-pressure)
        // mode is on the far side of the shift.
        let sigma = 0.9 * k_exact;
        let (k2, phi) = fem.mode_nearest(sigma * sigma, 40).unwrap();
        let f_fem = c * libm::sqrt(k2) / core::f64::consts::TAU;
        let f_exact = analytic_mode_frequency(dims, c, [1, 0, 0]);
        let rel = (f_fem - f_exact).abs() / f_exact;
        assert!(rel < 0.02, "fem {f_fem:.2} Hz vs exact {f_exact:.2} Hz");
        // Shape: the rigid-cavity (1,0,0) mode is cos(πx/Lx) along x
        // (zero normal gradient at both walls), sampled at the mid-y/z nodes.
        let [nx, ny, _nz] = fem.cells();
        let k_mid = _nz / 2;
        let mut xs = Vec::new();
        let mut vals = Vec::new();
        for i in 0..=nx {
            let node = i + (nx + 1) * (ny / 2 + (ny + 1) * k_mid);
            xs.push(i as f64 / nx as f64);
            vals.push(phi[node]);
        }
        let sign = if vals[0] < 0.0 { -1.0 } else { 1.0 };
        for (i, &v) in vals.iter().enumerate() {
            let x = xs[i];
            let expect = sign * libm::cos(core::f64::consts::PI * x);
            assert!(
                (v - expect).abs() < 0.05,
                "x={x:.3}: fem {v:.4} vs exact {expect:.4}"
            );
        }
    }

    #[test]
    fn analytic_frequencies_match_textbook() {
        // 25 °C air, 2.5 × 4 × 5 m room: (1,0,0) at 68.6 Hz.
        let f = analytic_mode_frequency([2.5, 4.0, 5.0], 343.0, [1, 0, 0]);
        assert!((f - 68.6).abs() < 0.1);
    }

    #[test]
    fn rejects_invalid_cavity() {
        assert!(CavityFem::box_cavity([0.0, 1.0, 1.0], [2, 2, 2]).is_err());
        assert!(CavityFem::box_cavity([1.0, 1.0, 1.0], [2, 0, 2]).is_err());
    }
}
