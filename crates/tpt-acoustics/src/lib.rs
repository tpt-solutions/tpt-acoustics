//! # tpt-acoustics (umbrella)
//!
//! Feature-gated re-exports of the `tpt-acoustics` workspace. Every domain
//! ships behind its own feature so downstreams pull exactly the crates they
//! need (all enabled by default):
//!
//! | feature | crate | contents |
//! | --- | --- | --- |
//! | `core` | [`tpt_acou_core`] | wave equations, media, unit-safe quantities |
//! | `ray` | [`tpt_acou_ray`] | ISM, Monte Carlo ray tracing, BVH scenes |
//! | `wave` | [`tpt_acou_wave`] | Helmholtz BEM/FEM, single-level FMM |
//! | `vibro` | [`tpt_acou_vibro`] | radiation efficiency, TL, STC, ATVs |
//! | `aero` | [`tpt_acou_aero`] | Lighthill, Curle, FW-H |
//! | `psycho` | [`tpt_acou_psycho`] | weightings, HRTF, Ambisonics, VBAP, loudness, masking |
//! | `verify` | [`tpt_acou_verify`] | proptest strategies and Kani harnesses |
//!
//! # Example
//!
//! Early reflections of a shoebox room and the Sabine reverberation time:
//!
//! ```
//! use tpt_acoustics::{ray, core};
//!
//! let room = ray::ism::ShoeboxRoom::uniform(
//!     [8.0, 5.0, 3.0],
//!     ray::material::BoundaryMaterial::concrete_painted(),
//! )
//! .unwrap();
//! let src = ray::vec3::Vec3::new(1.0, 1.0, 1.5);
//! let rcv = ray::vec3::Vec3::new(6.0, 4.0, 2.0);
//! let images = room.image_sources(src, 2).unwrap();
//! assert!(images.len() > 1);
//! assert!(room.sabine_rt60(1000.0) > 0.5);
//! ```
#![forbid(unsafe_code)]

#[cfg(feature = "aero")]
pub use tpt_acou_aero as aero;
#[cfg(feature = "core")]
pub use tpt_acou_core as core;
#[cfg(feature = "psycho")]
pub use tpt_acou_psycho as psycho;
#[cfg(feature = "ray")]
pub use tpt_acou_ray as ray;
#[cfg(feature = "verify")]
pub use tpt_acou_verify as verify;
#[cfg(feature = "vibro")]
pub use tpt_acou_vibro as vibro;
#[cfg(feature = "wave")]
pub use tpt_acou_wave as wave;
