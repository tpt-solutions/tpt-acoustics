//! # tpt-acou-psycho
//!
//! Modeling human sound perception for audio engineering and VR/AR.
//!
//! * [`weighting`] — standard frequency weightings: A/B/C (IEC 61672 analytic
//!   magnitude curves plus a stable digital biquad cascade via the bilinear
//!   transform of the real-pole ANSI prototype) and ITU-R BS.468-4 noise
//!   weighting (nominal-table interpolation).
//! * [`spectrum`] — FFT-based band levels and weighted levels built on
//!   `tpt-dsp-core` (github.com/tpt-solutions/tpt-dsp), which also supplies
//!   the Bark mapping, threshold-in-quiet and spreading function backing
//!   [`loudness`] and [`masking`].
//! * [`hrtf`] — spherical-head interaural time differences (Woodworth
//!   formula) and HRTF/HRIR grids with directional interpolation.
//! * [`ambisonics`] — SN3D B-format encoding to 2nd order and velocity
//!   decoding to arbitrary loudspeaker sets (least-squares decoder).
//! * [`vbap`] — Vector Base Amplitude Panning gains for 2-D and 3-D
//!   loudspeaker layouts.
//! * [`loudness`] — phon/sone conversion (ISO 532-1 mapping) and a simplified
//!   Bark-domain stationary loudness model (Zwicker-style power law with the
//!   standard 0.08 / 0.23 constants; the threshold-in-quiet is anchored with
//!   A-weighting — an engineering approximation, *not* a certified ISO 532-1
//!   implementation).
//! * [`masking`] — simultaneous-masking thresholds via the Johnston/Zwicker
//!   spreading function with tone and noise masker offsets.
#![no_std]
#![forbid(unsafe_code)]
extern crate alloc;

pub mod ambisonics;
pub mod hrtf;
pub mod loudness;
pub mod masking;
pub mod spectrum;
pub mod vbap;
pub mod weighting;
