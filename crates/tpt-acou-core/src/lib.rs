//! # tpt-acou-core
//!
//! The foundational layer of `tpt-acoustics`: governing equations of acoustics,
//! acoustic media, and unit-safe quantity types.
//!
//! This crate is `no_std` + `alloc` compatible: it never links `std`, so it
//! can be used on embedded (`thumbv6m-none-eabi`), WASM, and hosted targets
//! alike. Transcendental float math goes through the `libm` crate.
//!
//! # Contents
//!
//! * [`quantity`] — phantom-typed acoustic quantities (pressure, particle
//!   velocity, power, intensity) and sound pressure levels. Unit mismatches are
//!   compile errors.
//! * [`medium`] — fluids and solids: density, speed of sound, bulk modulus,
//!   characteristic impedance; standard media (air, water, …).
//! * [`impedance`] — complex-valued specific acoustic impedance and admittance,
//!   reflection factors and normal-incidence absorption.
//! * [`wave`] — frequency/wavenumber types, the acoustic wave equation and its
//!   standard reductions (Helmholtz, linearized Euler), plane and spherical
//!   waves, the free-field Green's function.
//! * [`bands`] — octave and third-octave band tables (ANSI/IEC nominal values).
//!
//! # Time-harmonic convention
//!
//! Throughout the `tpt-acoustics` family the time factor is
//! `exp(+i ω t)`: a wave travelling in the +x direction is
//! `p(x, t) = Re{ p̂ e^{i(ωt − k x)} }`, and the outgoing free-field Green's
//! function of the Helmholtz operator is `G(r) = e^{−ikr} / (4π r)`.
#![no_std]
#![forbid(unsafe_code)]
// The optional `materials` feature bridges to `tpt-mat-database`, which is
// std-linked (serde_json); everything else here stays pure core + alloc.
#[cfg(feature = "materials")]
extern crate alloc;

pub mod bands;
pub mod error;
pub mod impedance;
pub mod medium;
pub mod quantity;
pub mod vec3;
pub mod wave;

pub use error::AcousticError;
pub use quantity::{
    Intensity, ParticleVelocity, Power, Pressure, Quantity, SoundPressureLevel, Speed,
};
pub use vec3::Vec3;
pub use wave::{Frequency, Wavenumber};
