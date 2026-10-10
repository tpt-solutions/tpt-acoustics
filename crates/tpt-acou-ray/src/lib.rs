//! # tpt-acou-ray
//!
//! High-frequency geometric acoustics for room modeling and auralization.
//!
//! * [`ism`] — Image Source Method: exact early reflections for convex
//!   polyhedral (shoebox) rooms, with explicit reflection-order bookkeeping
//!   (the order bound is [`no_std`](crate)-safe and Kani-provable).
//! * [`ray`] — Monte Carlo energy ray tracing for late reverberation, with
//!   deterministic seeding (a small in-house xoshiro RNG — no `rand`
//!   dependency), Schroeder backward-integration decay curves and RT60
//!   estimates.
//! * [`bvh`] / [`scene`] — from-scratch BVH scene intersection over triangle
//!   meshes (Möller–Trumbore). This is the native replacement for the
//!   (not-yet-existing) `tpt-percept-cloud` spatial-index dependency.
//! * [`material`] — frequency-dependent absorption and scattering coefficients.
//! * [`differentiable`] — forward-mode dual-number hooks: gradients of the ISM
//!   transfer function with respect to room geometry, for inverse room design.
//!
//! All of this is `no_std` + `alloc` compatible.
#![no_std]
#![forbid(unsafe_code)]
extern crate alloc;

pub mod bvh;
pub mod differentiable;
pub mod ism;
pub mod material;
pub mod ray;
pub mod rng;
pub mod scene;
pub use tpt_acou_core::vec3;

pub use material::BoundaryMaterial;
pub use rng::Xoshiro256;
pub use scene::{Scene, SceneTriangle};
pub use vec3::Vec3;
