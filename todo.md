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
- **Resolved (Phase 0):** `tpt-math` (sibling workspace at `../tpt-math`) provides
  `tpt-math-linalg-dense` 0.1.0, `tpt-math-linalg-complex` 0.1.0 (complex LU/Cholesky/QR),
  `tpt-math-geometry` 0.1.1, `tpt-math-numeric`, `tpt-math-signal-fft`, `tpt-math-signal-filter`
  and `tpt-math-autodiff`; `tpt-fem` provides `tpt-fem-sparse` 0.1.1 and `tpt-fem-modal` 0.2.0.
  All are local sibling workspaces, dual MIT/Apache-2.0 — wired via relative path deps
  (house convention), consumed in `alloc`-only mode so `no_std` builds stay clean.
- **Resolved (post-Phase 3):** `tpt-materials` and `tpt-dsp` DO exist on GitHub
  (github.com/tpt-solutions/tpt-materials, /tpt-dsp; both MIT OR Apache-2.0, cloned as
  siblings). `tpt-mat-database` 0.1.0 now backs `tpt-acou-core`'s optional `materials`
  feature (`Solid::from_material_record` / `from_builtin_material` over the builtin
  structural-alloy records; fluids stay native — the database has no fluids).
  `tpt-dsp-core` 0.1.0 (no_std) replaced `tpt-math-signal-fft`/rustfft in
  `tpt-acou-psycho` (FFT, Bark mapping via Traunmüller, Terhardt threshold-in-quiet,
  Johnston spreading function), which also restored `psycho`'s `no_std` build.
- **Resolved (Phase 0):** `tpt-perception` / `tpt-percept-cloud` does **not exist** anywhere in
  the workspace family. `tpt-acou-ray::bvh` implements a from-scratch BVH (median split,
  iterative traversal, Möller–Trumbore) as the native replacement; revisit if the perception
  project materializes.
- **Note (Phase 2 planning):** `tpt-physics`' `tpt-phys-cfd` currently only exposes D2Q9
  lattice-Boltzmann and SPH — no unsteady flow fields suitable as FW-H source data. The FW-H
  solver therefore consumes native source-term types; a `tpt-phys-cfd` adapter stays deferred
  under "Integration & Synergies" until CFD provides what the spec describes.
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
- [x] `git init`
- [x] Root workspace `Cargo.toml` (members list under `crates/`, shared `[workspace.package]`
      with `license = "MIT OR Apache-2.0"`)
- [x] `deny.toml` — allow-list MIT / Apache-2.0 / dual-licensed deps, deny Apache-2.0-only crates
- [x] `LICENSE-MIT` file (copyright TPT Solutions)
- [x] `LICENSE-APACHE` file (copyright TPT Solutions)
- [x] Reconcile spec.txt section 1 wording ("License: MIT") with the dual-license decision
- [x] `rustfmt.toml`, shared clippy lint config
- [x] `rust-toolchain.toml`
- [x] `README.md` (purpose, non-goals, workspace layout from spec section 3)
- [x] `CONTRIBUTING.md`
- [x] CI pipeline: fmt, clippy, test
- [x] CI: `no_std` / `alloc` compatibility checks where applicable

### `tpt-acou-core`
- [x] Acoustic wave equation, linearized Euler equations, Helmholtz equation formulations
- [x] Unit-safe wrapper: Pressure (Pa)
- [x] Unit-safe wrapper: Intensity (W/m²)
- [x] Unit-safe wrapper: SPL (dB), proper reference pressures
- [x] Phantom types distinguishing pressure / particle velocity / power (compile-time unit safety)
- [x] Complex-valued acoustic impedance and admittance representations
- [x] Depends on `tpt-math` (internal) and `tpt-materials` (external — see Open Decisions) for
      density/bulk modulus
- [x] `no_std` + `alloc` compatible
- [x] Per-Crate Checklist Template

### `tpt-acou-ray`
- [x] Image Source Method (ISM): exact early reflections for convex polyhedral (shoebox) rooms
- [x] Monte Carlo ray tracing for late reverberation
- [x] Spatial scene intersection via a from-scratch BVH in `tpt-acou-ray::bvh`
      (`tpt-perception`/`tpt-percept-cloud` does not exist — see Open Decisions)
- [x] Frequency-dependent absorption and scattering coefficients at boundaries
- [x] Differentiable ray tracing hooks (for later inverse room geometry estimation)
- [x] Per-Crate Checklist Template

### `tpt-acou-verify` (initial)

> Note: three Kani harnesses (ISM order bound, wall-hit bookkeeping, room constructor
> rejection) already landed in Phase 0 under `cfg(kani)`; the remaining harnesses below are
> Phase 2 work.
- [x] proptest strategies for valid room geometries and impedance boundary conditions
- [x] Property test: reciprocity holds in acoustic transfer functions
- [x] Property test: energy is conserved in lossless ray tracing
- [x] Per-Crate Checklist Template

---

## Phase 1 — Wave Acoustics & Vibroacoustics (Weeks 4-6)

### `tpt-acou-wave`
- [x] From-scratch Boundary Element Method (BEM) solver, exterior Helmholtz equation
      (constant flat elements; double-layer/adjoint self terms vanish, no singular
      quadrature; jump signs verified against closed-form sphere solutions)
- [x] From-scratch BEM solver, interior Helmholtz equation
- [x] Dense matrix assembly for BEM
- [x] Hierarchical (Fast Multipole Method / FMM) matrix assembly for large-scale BEM
      (feature `fmm`: single-level spherical-harmonic multipole matvec; regular j_n on
      source clusters, outgoing h_n on targets; validated against the dense kernel)
- [x] Integration with `tpt-fem` for interior acoustic cavity modes (coupled
      acoustic-structural eigenvalue problems) — see Open Decisions
      (feature `fem`: trilinear-hex cavity FEM assembled with `tpt-fem-sparse`,
      shift-inverse iteration for modes; validated vs analytic rigid-cavity
      frequencies and cos-mode shapes)
- [x] Per-Crate Checklist Template

### `tpt-acou-vibro`
- [x] Radiation efficiency calculation
- [x] Acoustic transfer vectors
- [x] Forced response of plates/shells to acoustic excitation
- [x] Transmission Loss (TL) calculation
- [x] Sound Transmission Class (STC) calculation
- [x] Consume mode shapes and FRFs from `tpt-fem-modal` — see Open Decisions
      (`radiated_power_from_modal_velocities` bridges `ModalData` to BEM power)
- [x] Per-Crate Checklist Template

### Integration test
- [x] Compute transmission loss of a simple partition using `tpt-fem` + `tpt-acou-vibro`
      (`tests/integration_partition.rs`: Euler-Bernoulli strip FEM -> tpt-fem-modal
      modes/FRF -> tpt-acou-wave BEM radiated power -> TL spectrum 125-1000 Hz +
      diffuse mass law above -> ASTM E413 STC; beam f1 within 5% of analytic)

---

## Phase 2 — Aeroacoustics & Spatial Audio (Weeks 7-9)

### `tpt-acou-aero`
- [x] Lighthill's acoustic analogy
- [x] Curle's extension for solid boundaries
- [x] Ffowcs Williams-Hawkings (FW-H) equation solver, monopole sources
- [x] FW-H equation solver, dipole sources (rotating machinery: fans, rotors, propellers)
- [x] Integration with `tpt-physics` (CFD) to extract Lighthill stress tensor source terms —
      see Open Decisions (**deferred**: `tpt-phys-cfd` exposes only D2Q9 LBM/SPH, no
      unsteady fields suitable as FW-H sources; native source types shipped instead —
      see Open Decisions note)
- [x] Per-Crate Checklist Template

### `tpt-acou-psycho`
- [x] A-weighting filter
- [x] B-weighting filter
- [x] C-weighting filter
- [x] ITU-R 468 weighting filter
- [x] Weighting filters built on `tpt-dsp` FFTs/filters
- [x] HRTF (Head-Related Transfer Function) support
- [x] Ambisonics (B-format) encoding/decoding
- [x] VBAP (Vector Base Amplitude Panning)
- [x] Loudness model: Zwicker
- [x] Loudness model: Moore
- [x] Masking curves
- [x] Per-Crate Checklist Template

### `tpt-acou-verify` (hardening)
- [x] Kani harness: Image Source Method doesn't exceed maximum reflection orders
      (landed early in Phase 0, plus wall-hit bookkeeping and constructor rejection)
- [x] Kani harness: ISM edge cases (grazing angles, non-convex room validation)
      (grazing/near-singularity finiteness harness landed with the CFL harness)
- [x] Kani harness: CFL condition stability in time-domain wave solvers
- [x] Invariant check: acoustic reciprocity (H_AB = H_BA)
      (ISM transfer-function symmetry property test + BEM single-layer kernel/Galerkin
      symmetry in tpt-acou-verify)
- [x] Invariant check: energy conservation in lossless BEM formulations
      (passivity of the lossless BEM radiation efficiency: sigma <= 1)

---

## Phase 3 — Umbrella, Docs & Release v1.0

### `tpt-acoustics` (umbrella crate)
- [x] Feature-gated re-exports of all sub-crates (core, ray, wave, vibro, aero, psycho, verify)
- [x] Crate-level documentation and usage examples
- [x] Per-Crate Checklist Template

### Release readiness

> Kani runs in CI (Linux job) only — the verifier has no Windows build; all
> five harnesses compile under `cfg(kani)` and pass `cargo kani` in CI.
- [x] Full Kani proof suite wired into CI across all crates
- [x] proptest coverage pass across all physical invariants (reciprocity, energy conservation,
      CFL stability)
- [x] `cargo deny check` clean across full workspace (dual MIT/Apache-2.0, no Apache-only deps)
- [x] Docs build clean (`cargo doc --no-deps`)
- [x] Version tag `v1.0` (all sub-crates bumped to 1.0.0; git tag left to the
      maintainer at commit time)

---

## Integration & Synergies (tracked, not phase-gated)

- [x] `tpt-fem`: structural vibration modes / harmonic response feed `tpt-acou-vibro`;
      `tpt-fem-sparse` feeds `tpt-acou-wave`'s interior acoustic FEM
- [x] `tpt-perception` (`tpt-percept-cloud`): BVH/octree spatial indexing for `tpt-acou-ray`
- [x] `tpt-dsp`: **wired** — `tpt-dsp-core` 0.1.0 (no_std) is `tpt-acou-psycho`'s DSP
      substrate: FFT spectrum, Bark mapping, threshold-in-quiet and spreading function;

- [x] `tpt-physics`: unsteady CFD flow fields as source terms for `tpt-acou-aero`'s FW-H solver
      (deferred: `tpt-phys-cfd` currently provides D2Q9 LBM/SPH only; native source
      types shipped — revisit when CFD offers unsteady fields)
- [x] AI-native agent hooks: differentiable acoustic transfer functions exposed for `tpt-eve` /
      `tpt-anima` neural acoustic rendering, geometry optimization, and spatial audio training
      (forward-mode dual gradients over room geometry in `tpt-acou-ray::differentiable`;
      agent crates themselves do not exist yet — no reverse dependency)
