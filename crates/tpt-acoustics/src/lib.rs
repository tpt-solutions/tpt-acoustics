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
//! | `verify` | `tpt_acou_verify` | proptest strategies and Kani harnesses |
//!
//! # Example
//!
//! Early reflections of a shoebox room and the Sabine reverberation time:
//!
//! ```
//! use tpt_acoustics::ray;
//!
//! let room = ray::ism::ShoeboxRoom::uniform(
//!     [8.0, 5.0, 3.0],
//!     ray::material::BoundaryMaterial::concrete_painted(),
//! )
//! .unwrap();
//! let src = ray::vec3::Vec3::new(1.0, 1.0, 1.5);
//! let images = room.image_sources(src, 2).unwrap();
//! assert!(images.len() > 1);
//! assert!(room.sabine_rt60(1000.0) > 0.5);
//! ```
//!
//! # Quick start
//!
//! `use tpt_acoustics::prelude::*` brings in the common types. This program
//! builds an energy echogram, turns it into a 48 kHz impulse response, and
//! reports clarity:
//!
//! ```
//! # fn main() -> Result<(), tpt_acoustics::prelude::AcousticError> {
//! use tpt_acoustics::prelude::*;
//!
//! let room = ShoeboxRoom::uniform([8.0, 5.0, 3.0], BoundaryMaterial::concrete_painted())?;
//! let air = Medium::air_20c();
//! let src = Vec3::new(1.5, 1.2, 1.4);
//! let rcv = Vec3::new(6.2, 3.7, 1.6);
//!
//! let echogram = room.echogram(src, rcv, 1000.0, 12, &air)?;
//! let fs = 48_000.0;
//! let rir = rir_from_echogram(&echogram, fs, 1)?;
//! let c80 = clarity_db(&rir, fs, 80.0).expect("late energy present");
//! assert!(c80.is_finite());
//! # Ok(())
//! # }
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

/// Builder for a validated shoebox source–receiver setup. See [`builder`].
#[cfg(all(feature = "core", feature = "ray"))]
pub mod builder;

/// Common types for a first program: `use tpt_acoustics::prelude::*;`.
///
/// Available when both `core` and `ray` are enabled (the default).
#[cfg(all(feature = "core", feature = "ray"))]
pub mod prelude {
    pub use crate::builder::{RoomBuilder, Simulation};
    pub use tpt_acou_core::error::AcousticError;
    pub use tpt_acou_core::medium::Medium;
    pub use tpt_acou_core::quantity::{Pressure, SoundPressureLevel};
    pub use tpt_acou_core::vec3::Vec3;
    pub use tpt_acou_ray::ism::ShoeboxRoom;
    pub use tpt_acou_ray::iso3382::{clarity_db, definition, reverberation_times};
    pub use tpt_acou_ray::material::BoundaryMaterial;
    pub use tpt_acou_ray::ray::{trace_shoebox, RayTracerConfig};
    pub use tpt_acou_ray::rir::{rir_from_echogram, write_wav_f32};
}
