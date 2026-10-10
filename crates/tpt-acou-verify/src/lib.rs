//! # tpt-acou-verify
//!
//! Mathematical verification of acoustic algorithms across the workspace:
//!
//! * [`strategies`] — reusable `proptest` strategies generating valid room
//!   geometries, boundary materials, impedance values and observation points.
//! * Kani bounded-model-checking harnesses in the `proofs` module (compiled
//!   only under `cfg(kani)`; run via `cargo kani` in CI).
//!
//! Physical-law property tests (acoustic reciprocity, energy conservation in
//! lossless ray tracing, Green's-function scaling) live in the private
//! `properties` module, compiled with `cargo test`.
//!
//! This crate is deliberately std-linked (proptest needs it); the crates it
//! verifies stay `no_std`-clean.
#![forbid(unsafe_code)]

pub mod strategies;

#[cfg(test)]
mod properties;

/// Kani proofs (compiled only by the Kani verifier under `cfg(kani)`).
#[cfg(kani)]
pub mod proofs;
