//! # tpt-acou-wave
//!
//! Low-to-mid frequency wave acoustics: from-scratch frequency-domain solvers
//! for the Helmholtz equation.
//!
//! * [`mesh`] — boundary meshes (triangulated surfaces) and element views.
//! * [`bem`] — collocation Boundary Element Method for interior/exterior
//!   Helmholtz problems with *constant flat elements*: the double-layer and
//!   adjoint-double-layer self terms vanish identically (the kernels are odd
//!   in the element plane), so no singular quadrature is required anywhere.
//! * [`fem_cavity`] *(feature `fem`, default)* — interior acoustic cavity FEM
//!   (trilinear hexahedra) assembled with `tpt-fem-sparse`, with shift-inverse
//!   iteration for cavity modes.
//! * `fmm` *(feature `fmm`)* — single-level multipole acceleration of the
//!   BEM single-layer matvec for larger meshes.
//!
//! # Conventions
//!
//! Time factor `exp(+iωt)`; the outgoing free-field Green's function is
//! `G(r) = e^{−ikr}/(4πr)` (see `tpt-acou-core::wave`). Surface normals point
//! *out of the interior domain*. Complex amplitudes are RMS-compatible:
//! time-averaged intensity is `½·Re(p·v*)`.
//!
//! # Known limitations (documented, not hidden)
//!
//! * Exterior Dirichlet/Neumann collocation equations are non-unique at the
//!   interior eigenfrequencies of the complementary domain; CHIEF points /
//!   Burton–Miller regularization are not implemented — choose `k` away from
//!   those frequencies (the validation tests do).
//! * The FMM is a single-level multipole scheme: it accelerates the
//!   single-layer *matvec* (assembly-by-matvec), not the full factorization.
#![no_std]
#![forbid(unsafe_code)]
extern crate alloc;

pub mod bem;
pub mod krylov;
pub mod mesh;

#[cfg(feature = "fem")]
pub mod fem_cavity;

#[cfg(feature = "fmm")]
pub mod fmm;

pub use bem::{BemError, HelmholtzBem};
pub use mesh::{BoundaryMesh, ElementView};
