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
- [ ] Loudness model: Moore (unticked — not implemented; see Phase 4 Moore item)
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

---

## Phase 4 — Review Follow-ups (platform review, 2026-10-11)

> Source: platform review of v1.0. Items marked *(verify)* were found by reading the code plus
> domain knowledge and need confirming against the cited standard before changing. Full write-up:
> the review plan (bugs A, gaps B, ideas C, adoption D).

### 4.1 Build & CI (do first)
- [x] Make the workspace build from a clean clone: the sibling path deps are now git deps pinned by
      `rev` to `tpt-math`, `tpt-fem`, `tpt-dsp` and `tpt-materials` on GitHub (`Cargo.toml`), and
      `deny.toml` allows those git sources. Pinned to the committed `HEAD` of each sibling. Note:
      `tpt-dsp` has uncommitted lint-only edits locally, which are not in the pinned revision.
- [x] Clean-clone verification: copied the project (without `target/` or `.git/`) to an empty directory
      with no sibling checkouts and ran `cargo build --workspace`, which succeeds. This copied the working
      tree, not a fresh `git clone`; that needs this change committed first.
- [x] Make `verify` opt-in on the umbrella crate (removed from `default`; `proptest` only with `verify`)
- [x] Reconsider `MPL-2.0` in `deny.toml` allow-list vs the MIT/Apache-only philosophy in spec.txt
      (removed: no crate in the graph uses it, and `cargo deny check licenses` still passes)
- [x] `no_std` CI: cover remaining crates where applicable. Reviewed: spec.txt only requires no_std
      for the core layer, and CI already builds the five crates that declare `#![no_std]`. `tpt-acou-vibro`
      is std-only and `tpt-acou-verify` depends on `proptest` (std), so neither is added.

### 4.2 Correctness bugs
- [x] **STC contour wrong** (`tpt-acou-vibro/src/transmission.rs`): replaced with the ASTM E413
      contour (`STC_CONTOUR_REL_DB`); tests are analytic (the contour itself rates 52, flat 60 dB
      rates 60). Not yet checked against a published STC worked example.
- [ ] **ITU-R 468 table** (`tpt-acou-psycho/src/weighting.rs`) *(partly done)*: 3.15 kHz (+9.0),
      4 kHz (+10.5), 5 kHz (+11.7), 6.3 kHz (+12.2, the peak), 8 kHz (+11.4) and 10 kHz (+8.1) taken
      from secondary mirrors of the Wikipedia article, which agree with each other. **Still open:**
      confirm all points against the BS.468-4 text, since 31.5 Hz–2.5 kHz and 12.5–20 kHz are
      unconfirmed. Quasi-peak detector not addressed.
- [ ] **Moore loudness** — unticked in Phase 2 (not implemented). Moore-Glasberg (ISO 532-2) is open.
- [x] `zwicker_stationary`: spectral spreading added (Zwicker's 27 dB/Bark lower slope and
      level-dependent upper slope). A 1 kHz / 40 dB tone now rates 1.02 sone against 1 sone in
      ISO 532-1; tests pin it within 5 %. Negative SPL allowed, Terhardt wording fixed. The model
      is still a simplified approximation, not a certified ISO 532-1 implementation.
- [x] **Ray tracer** (`tpt-acou-ray/src/ray.rs`): particles keep tracing after capture; detections
      are a separate tally and the partition is `injected = absorbed + remaining`. RT60 test
      tightened from 35% to 10% and passes.
- [x] Ray tracer input validation: `RayTracerConfig::validate`, source and receiver inside the room
      for shoebox tracing; `trace_shoebox` and `trace_scene` return `Result`.
- [x] **FW-H** (`tpt-acou-aero/src/fwh.rs`): `rotating_ring_pressure` now includes the `1/(1-M_r)`
      factor with the emission-time blade velocity; `omega = 0` returns 0; `rotating_monopole_ring`
      dropped its unused `omega` parameter (static snapshot, **semver-breaking**); monopole, dipole
      and ring derivatives are central differences; `a = 0` and `1 - M_r <= 0` return NaN.
      Still open: an analytic rotor-noise benchmark (no published reference is wired in yet).
- [x] `b_weighting_db` doc says +1.96 dB, matching the code.
- [x] `ShoeboxRoom::echogram` takes a `&Medium` and uses its speed of sound (no more hard-coded 343).
- [x] `SoundPressureLevel::try_from_pressure_rms` / `try_from_intensity` / `try_from_power` return
      `Result`; the panicking constructors are unchanged.
- [x] Umbrella crate doc example: removed unused `core` import and unused `rcv`.
- [x] Confirm what the ticked CFL Kani harness proves: it checks only that the float `c·dt/dx <= 1`
      test agrees with the integer rational bound. It exercises no solver; there is no time-domain
      solver in the workspace, so it is a classification lemma, not a stability proof of a scheme.

### 4.3 Missing features vs spec
- [ ] GMRES (iterative) solver in `tpt-acou-wave` consuming the FMM matvec, with a preconditioner *(partly done)*: `tpt_acou_wave::krylov::gmres` is restarted, right-preconditioned complex GMRES over any matvec closure, with a Jacobi option. It matches the dense BEM exterior solve. **Still open:** wiring it to the FMM operator on boundary nodes.
- [ ] ISM with wall phase / complex impedance (low-frequency / modal behaviour)
- [ ] Multi-parameter differentiable ray/ISM gradients (dimensions, per-band absorption,
      source/receiver position); evaluate `tpt-math-autodiff` / reverse mode
- [ ] Room impulse response (RIR) generation from echogram + ray tail *(partly done)*:
      `tpt_acou_ray::rir::rir_from_echogram` builds an RIR from the ISM echogram with random signs
      (energy-correct, not phase-correct). The ray-tail hybrid is open: the ray and ISM energy
      normalisations differ, so they need reconciling first.
- [ ] Optional `serde` support, `rayon` parallel feature (ray tracing, BEM assembly) *(wasm32 CI target added: `no_std` crates build for `wasm32-unknown-unknown`; `serde` and `rayon` still open)*
- [x] criterion benchmarks (`crates/tpt-acou-ray/benches/ray.rs`: ray tracer, echogram, RIR). **Still open:** regression tracking (baseline storage and comparison in CI).

### 4.4 New features (innovation)
- [x] ISO 3382 room metrics from an RIR: T20/T30/EDT, C50/C80, D50 (`tpt_acou_ray::iso3382`).
      Analytic tests on an exponential decay. **Still open:** STI and G.
- [x] RIR export to WAV (`tpt_acou_ray::rir::write_wav_f32`, 32-bit float mono). **Still open:**
      convolution-based auralization.
- [ ] Differentiable inverse design: fit room to target RT60/C80; optimize panel absorption
- [ ] Hybrid broadband RIR: ISM early + rays late + wave low-frequency, crossed at the Schroeder frequency
- [x] STL import (`tpt_acou_ray::stl`, ASCII and binary, detected by size) and OBJ import (`tpt_acou_ray::obj`: vertices, polygon faces, slash and negative indices, `usemtl` groups; one material per scene). **Still open:** an octave-band material database (ISO 354 style).
- [ ] Standards compliance report (TL, STC, OITC, Lw, ISO 3382, ISO 532)
- [ ] Uncertainty quantification via seeded Monte Carlo over material uncertainty
- [ ] Python (PyO3) and WASM bindings

### 4.5 Adoption: examples, docs, templates
- [x] `examples/` in `crates/tpt-acoustics/examples/`: `shoebox_rt60`, `room_impulse_response`,
      `partition_stc`, `sphere_radiation_bem`, `fan_noise_fwh`, `optimize_room_gradient`,
      `loudness_a_weighting`, `spatial_audio_vbap`. Each asserts against an analytic or literature
      value. Listed in the README.
- [x] `tpt_acoustics::prelude` + `RoomBuilder` / `Simulation` builder API (`crates/tpt-acoustics/src/builder.rs`).
      The builder validates dimensions, source and receiver once; `Simulation` gives the echogram,
      impulse response and Sabine T60. Still open: a `Simulation` path to ray-traced T30 and the wave solvers.
- [ ] README: quick start, does/doesn't table and per-crate validation status are done; badges still open.
- [x] Cookbook as one page (`docs/cookbook.md`): solver-choice table, units and conventions, glossary, recipes. **Still open:** the mdBook layout around it.
- [x] Validation page (`docs/validation.md`): every reference with its tolerance and test name. Open items are listed on the page.
- [ ] Templates: room study done (`templates/room-study.toml`, run with `tpt-acou study`). **Still open:** partition-study and rotor-noise-study templates; partition and rotor studies need their own command.
- [x] `tpt-acou-cli` crate (`crates/tpt-acou-cli`, binary `tpt-acou`): `rt60`, `stc` (mass or 16 TL values), `weight` and `study <file.toml>` subcommands, with exit status 2 for usage errors and 1 for failed computations.

### 4.6 Automation & hygiene
- [x] Dependabot (`.github/dependabot.yml`, cargo + actions; Renovate not used). `cargo semver-checks`, `release-plz` + CHANGELOG.md still open
- [x] Coverage job (`cargo llvm-cov --summary-only`), MSRV job (1.85, libraries), Windows + macOS test matrix *(added to CI; not yet run on GitHub, so the macOS and MSRV jobs are unverified)*
- [ ] CI check that ticked todo items reference existing symbols (prevents the Moore-style drift)
- [x] Issue and PR templates (`.github/ISSUE_TEMPLATE/`, `.github/pull_request_template.md`)
