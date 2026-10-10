//! # tpt-acou-aero
//!
//! Prediction of flow-generated noise via the acoustic analogies, implemented
//! from scratch (no CFD dependency in the core — source data is consumed
//! through plain numeric types; the `tpt-physics` adapter stays deferred
//! under "Integration & Synergies" in todo.md, see Open Decisions).
//!
//! * [`lighthill`] — Lighthill's acoustic analogy: the quadrupole source
//!   tensor `T_ij = rho u_i u_j + (p - c0^2 rho') delta_ij - tau_ij` and the
//!   far-field density from a volume distribution of sources.
//! * [`curle`] — Curle's extension: the compact dipole radiated by unsteady
//!   surface forces (vortex shedding, bluff bodies).
//! * [`fwh`] — Ffowcs Williams–Hawkings far-field forms for moving monopole
//!   (thickness) and dipole (loading) sources, including rotating machinery:
//!   blade-passing-frequency tones from rotating point sources with the
//!   exact retarded-time geometry.
//!
//! # Conventions
//!
//! Time factor `exp(+i omega t)` for the frequency-domain helpers; the
//! time-domain kernels use plain retarded times `t_r = t - a/c0`. Far-field
//! compact-source forms are valid for observers at `a >> source extent`.

#![no_std]
#![forbid(unsafe_code)]
extern crate alloc;

pub mod curle;
pub mod fwh;
pub mod lighthill;
pub mod source;
