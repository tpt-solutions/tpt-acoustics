# tpt-acoustics — Phased Build Tracker

> Pure-Rust, AI-native computational acoustics library: geometric room acoustics, wave acoustics
> (BEM/FEM Helmholtz), aeroacoustics, vibroacoustics, and psychoacoustic modeling — no C/C++ FFI,
> no external solver dependencies. Dual-licensed **MIT OR Apache-2.0**. © TPT Solutions.
> See [spec.txt](spec.txt) for the full design.

## Open Decisions

- **Checkpoint:** spec.txt section 1 currently reads "License: MIT ... no Apache-2.0-only crates".
  The project is dual-licensed **MIT OR Apache-2.0** (matching every other `tpt-*` project); the
  "no Apache-2.0-only" *dependency* rule stays as-is (still enforced via `deny.toml`), but the
  spec's own license line needs updating to avoid reading as MIT-only. See Phase 0.
- **Checkpoint:** `tpt-materials`, `tpt-fem` (+ `tpt-fem-sparse`, `tpt-fem-modal`), and
  `tpt-physics` are external dependencies (crates.io / GitHub, not yet vendored into this
  workspace). Confirm exact crate names, versions, and repo URLs before wiring them into any
  `Cargo.toml`, and re-check their own license (must be MIT or MIT/Apache-2.0 dual — no
  Apache-2.0-only) before adding.
- **Checkpoint:** spec.txt refers to `tpt-percept` / `tpt-percept-cloud`; the actual sibling
  project is `tpt-perception`, whose sub-crate is `tpt-percept-cloud`. Depend on
  `tpt-perception`'s `tpt-percept-cloud` crate (BVH/k-d tree spatial indexing) for
  `tpt-acou-ray`'s scene intersection.
- **Checkpoint:** `tpt-eve` / `tpt-anima` (AI-native agents mentioned in spec section 5) do not
  exist yet in the workspace. `tpt-acoustics` exposes differentiable hooks for them to consume
  later; no dependency in the other direction, so this doesn't block any phase here.

## Per-Crate Checklist Template

Applied to every crate below (referenced by name, not repeated per crate):

1. Scaffold crate under `crates/`, add to workspace `Cargo.toml`.
2. Add only required dependencies (internal `tpt-*` crates or approved MIT/Apache-2.0 deps —
   check `deny.toml`).
3. Implement the feature(s) listed for the crate.
4. Unit tests for the feature(s).
5. `rustdoc` on all public items.
6. `cargo fmt` + `cargo clippy --all-targets -- -D warnings` clean.
7. `cargo deny check` clean (no Apache-2.0-only deps).
8. `no_std` + `alloc` compatibility verified (where applicable per spec).

---

## Phase 0 — Foundation & Geometric Acoustics (Weeks 1-3)

### Workspace scaffolding
- [ ] `git init`
- [ ] Root workspace `Cargo.toml` (members list under `crates/`, shared `[workspace.package]`
      with `license = "MIT OR Apache-2.0"`)
- [ ] `deny.toml` — allow-list MIT / Apache-2.0 / dual-licensed deps, deny Apache-2.0-only crates
- [ ] `LICENSE-MIT` file (copyright TPT Solutions)
- [ ] `LICENSE-APACHE` file (copyright TPT Solutions)
- [ ] Reconcile spec.txt section 1 wording ("License: MIT") with the dual-license decision
- [ ] `rustfmt.toml`, shared clippy lint config
- [ ] `rust-toolchain.toml`
- [ ] `README.md` (purpose, non-goals, workspace layout from spec section 3)
- [ ] `CONTRIBUTING.md`
- [ ] CI pipeline: fmt, clippy, test
- [ ] CI: `no_std` / `alloc` compatibility checks where applicable

### `tpt-acou-core`
- [ ] Acoustic wave equation, linearized Euler equations, Helmholtz equation formulations
- [ ] Unit-safe wrapper: Pressure (Pa)
- [ ] Unit-safe wrapper: Intensity (W/m²)
- [ ] Unit-safe wrapper: SPL (dB), proper reference pressures
- [ ] Phantom types distinguishing pressure / particle velocity / power (compile-time unit safety)
- [ ] Complex-valued acoustic impedance and admittance representations
- [ ] Depends on `tpt-math` (internal) and `tpt-materials` (external — see Open Decisions) for
      density/bulk modulus
- [ ] `no_std` + `alloc` compatible
- [ ] Per-Crate Checklist Template

### `tpt-acou-ray`
- [ ] Image Source Method (ISM): exact early reflections for convex polyhedral (shoebox) rooms
- [ ] Monte Carlo ray tracing for late reverberation
- [ ] Spatial scene intersection via `tpt-perception`'s `tpt-percept-cloud` (BVH/k-d tree)
- [ ] Frequency-dependent absorption and scattering coefficients at boundaries
- [ ] Differentiable ray tracing hooks (for later inverse room geometry estimation)
- [ ] Per-Crate Checklist Template

### `tpt-acou-verify` (initial)
- [ ] proptest strategies for valid room geometries and impedance boundary conditions
- [ ] Property test: reciprocity holds in acoustic transfer functions
- [ ] Property test: energy is conserved in lossless ray tracing
- [ ] Per-Crate Checklist Template

---

## Phase 1 — Wave Acoustics & Vibroacoustics (Weeks 4-6)

### `tpt-acou-wave`
- [ ] From-scratch Boundary Element Method (BEM) solver, exterior Helmholtz equation
- [ ] From-scratch BEM solver, interior Helmholtz equation
- [ ] Dense matrix assembly for BEM
- [ ] Hierarchical (Fast Multipole Method / FMM) matrix assembly for large-scale BEM
- [ ] Integration with `tpt-fem` for interior acoustic cavity modes (coupled
      acoustic-structural eigenvalue problems) — see Open Decisions
- [ ] Per-Crate Checklist Template

### `tpt-acou-vibro`
- [ ] Radiation efficiency calculation
- [ ] Acoustic transfer vectors
- [ ] Forced response of plates/shells to acoustic excitation
- [ ] Transmission Loss (TL) calculation
- [ ] Sound Transmission Class (STC) calculation
- [ ] Consume mode shapes and FRFs from `tpt-fem-modal` — see Open Decisions
- [ ] Per-Crate Checklist Template

### Integration test
- [ ] Compute transmission loss of a simple partition using `tpt-fem` + `tpt-acou-vibro`

---

## Phase 2 — Aeroacoustics & Spatial Audio (Weeks 7-9)

### `tpt-acou-aero`
- [ ] Lighthill's acoustic analogy
- [ ] Curle's extension for solid boundaries
- [ ] Ffowcs Williams-Hawkings (FW-H) equation solver, monopole sources
- [ ] FW-H equation solver, dipole sources (rotating machinery: fans, rotors, propellers)
- [ ] Integration with `tpt-physics` (CFD) to extract Lighthill stress tensor source terms —
      see Open Decisions
- [ ] Per-Crate Checklist Template

### `tpt-acou-psycho`
- [ ] A-weighting filter
- [ ] B-weighting filter
- [ ] C-weighting filter
- [ ] ITU-R 468 weighting filter
- [ ] Weighting filters built on `tpt-dsp` FFTs/filters
- [ ] HRTF (Head-Related Transfer Function) support
- [ ] Ambisonics (B-format) encoding/decoding
- [ ] VBAP (Vector Base Amplitude Panning)
- [ ] Loudness model: Zwicker
- [ ] Loudness model: Moore
- [ ] Masking curves
- [ ] Per-Crate Checklist Template

### `tpt-acou-verify` (hardening)
- [ ] Kani harness: Image Source Method doesn't exceed maximum reflection orders
- [ ] Kani harness: ISM edge cases (grazing angles, non-convex room validation)
- [ ] Kani harness: CFL condition stability in time-domain wave solvers
- [ ] Invariant check: acoustic reciprocity (H_AB = H_BA)
- [ ] Invariant check: energy conservation in lossless BEM formulations

---

## Phase 3 — Umbrella, Docs & Release v1.0

### `tpt-acoustics` (umbrella crate)
- [ ] Feature-gated re-exports of all sub-crates (core, ray, wave, vibro, aero, psycho, verify)
- [ ] Crate-level documentation and usage examples
- [ ] Per-Crate Checklist Template

### Release readiness
- [ ] Full Kani proof suite wired into CI across all crates
- [ ] proptest coverage pass across all physical invariants (reciprocity, energy conservation,
      CFL stability)
- [ ] `cargo deny check` clean across full workspace (dual MIT/Apache-2.0, no Apache-only deps)
- [ ] Docs build clean (`cargo doc --no-deps`)
- [ ] Version tag `v1.0`

---

## Integration & Synergies (tracked, not phase-gated)

- [ ] `tpt-fem`: structural vibration modes / harmonic response feed `tpt-acou-vibro`;
      `tpt-fem-sparse` feeds `tpt-acou-wave`'s interior acoustic FEM
- [ ] `tpt-perception` (`tpt-percept-cloud`): BVH/octree spatial indexing for `tpt-acou-ray`
- [ ] `tpt-dsp`: FFTs and digital filters for `tpt-acou-psycho` weighting filters and for
      converting time-domain aeroacoustic pressure signals to frequency spectra
- [ ] `tpt-physics`: unsteady CFD flow fields as source terms for `tpt-acou-aero`'s FW-H solver
- [ ] AI-native agent hooks: differentiable acoustic transfer functions exposed for `tpt-eve` /
      `tpt-anima` neural acoustic rendering, geometry optimization, and spatial audio training
