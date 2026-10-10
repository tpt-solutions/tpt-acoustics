# Changelog

All notable changes to this project are recorded here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions follow
[Semantic Versioning](https://semver.org/). Releases are cut by release-plz from
the commit history; see `release-plz.toml`.

## [Unreleased]

### Added

- **Rooms and ray tracing:** ISO 3382 metrics (T20, T30, EDT, C50, C80, D50, G),
  Schroeder frequency, Speech Transmission Index (IEC 60268-16 structure, male
  and female weights), sample-based room impulse responses from the echogram,
  a hybrid impulse response with the ray tail calibrated at the crossover, a
  broadband hybrid (modal below the Schroeder frequency, image sources above),
  auralization by FFT convolution, and float WAV export.
- **Rigid-box modal transfer function**, with reciprocity and resonance checks.
- **Impedance walls** in the image-source transfer function (oblique complex
  reflection from the unfolded image path).
- **Design:** analytic Sabine time and gradients with respect to room dimensions
  and per-wall absorption, uniform and layout-preserving absorption inversion,
  and a seeded Monte Carlo uncertainty estimate.
- **Octave-band absorption database** (three finishes, each from a named chart).
- **Geometry import:** Wavefront OBJ and STL (ASCII and binary) into ray-tracing
  scenes.
- **Parallel ray tracing** behind the `parallel` feature (rayon). Every ray has
  its own seeded stream, so results do not depend on the thread count.
- **Boundary elements:** restarted, right-preconditioned complex GMRES; the
  adjoint double-layer operator applied through the single-level FMM; and the
  exterior Neumann solve through that operator.
- **Psychoacoustics:** Zwicker spreading; Moore–Glasberg tonal loudness
  (calibrated to the 1 kHz / 40 dB anchor); ITU-R BS.468-4 table completed from
  the official response column.
- **Vibroacoustics:** ASTM E413 STC contour corrected.
- **Umbrella crate:** `prelude`, `RoomBuilder` and `Simulation` (echogram,
  impulse response, T30, STI, G), and examples for each domain.
- **Command line** (`tpt-acou`): `rt60`, `stc`, `weight`, `study`, `design`,
  `uq` and `auralize`. Studies cover room, partition, rotor and compliance.
- **Optional serde** derives on core value types (`serde` feature).
- **Tooling:** criterion benchmarks; a `todo.md` check that ticked items name
  code that exists; Dependabot; issue and PR templates; release-plz and
  semver-check configuration; mdBook layout over `docs/`.

### Changed

- The ray tracer returns `Result`, and every ray has its own seeded stream
  (results differ from 1.0 for the same seed but keep the same statistics).
- `rotating_monopole_ring` no longer takes an unused `omega` parameter
  (**breaking**).
- `ShoeboxRoom::echogram` takes a `&Medium` argument (**breaking**).
- `SoundPressureLevel` gained fallible `try_` constructors.
- The `verify` feature is no longer a default feature of the umbrella crate.
- Sibling `tpt-*` dependencies are pinned git revisions, not relative paths.

### Fixed

- Ray particles no longer stop at their first receiver capture, which had biased
  late decay and T30.
- The FW-H rotating-ring pressure now includes the `1/(1 − M_r)` factor, and it
  guards `omega = 0` and undefined far-field points.
- STC ratings use the ASTM E413 contour shape.
- The ITU-R BS.468-4 table had the wrong 3.15 kHz value and was missing points.
- The Zwicker documentation overstated the 1 kHz / 40 dB loudness; it is now
  reported as measured.

## [1.0.0]

Initial release of the workspace: geometric acoustics, wave acoustics
(BEM/FEM), vibroacoustics, aeroacoustics and psychoacoustics, with unit,
property and integration tests.
