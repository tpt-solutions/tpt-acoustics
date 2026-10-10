//! # tpt-acou-vibro
//!
//! Coupling between structural vibration and acoustic radiation.
//!
//! * [`plate`] — analytic simply-supported thin-plate modal model
//!   (Kirchhoff–Love theory): natural frequencies, mode shapes, modal masses,
//!   and point harmonic response.
//! * [`radiation`] — baffled-piston radiation resistance (closed form via the
//!   Bessel-function series) and radiation efficiency computed through the
//!   `tpt-acou-wave` BEM for arbitrary surface velocity patterns.
//! * [`atv`] — Acoustic Transfer Vectors: the linear map from surface normal
//!   velocity to field-point pressure, built on the Helmholtz Green's
//!   function.
//! * [`transmission`] — transmission loss (normal- and diffuse-field mass
//!   laws, exact diffuse integral) and the ASTM E413 Sound Transmission Class
//!   (STC) contour rating.
//!
//! Modal data produced by `tpt-fem`'s `tpt-fem-modal` (sparse generalized
//! eigenanalysis of a structural FE model) plugs directly into
//! [`radiation::radiated_power_from_modal_velocities`], completing the
//! spec.txt §5 `tpt-fem → tpt-acou-vibro` synergy (the Phase-1 integration
//! test exercises the full chain).
#![forbid(unsafe_code)]

pub mod atv;
pub mod plate;
pub mod radiation;
pub mod transmission;
