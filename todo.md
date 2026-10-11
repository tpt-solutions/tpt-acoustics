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
- [x] Loudness model: Moore (`tpt-acou-psycho/src/loudness_moore.rs`; see Phase 4 Moore item for its validation status)
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

> Source: platform review of v1.0. Status reflects the code as of the last commit plus the uncommitted
> work in this round. "Validated" means a test checks the result against a closed form, an independent
> numerical method, or a cited value. Items marked *(blocked)* need a source the repository does not have;
> they are listed under "Blockers" at the end of this phase.

### 4.1 Build & CI (do first)
- [x] Clean-clone build: sibling path deps are pinned git revisions (`Cargo.toml`), and `deny.toml` allows
      those git sources. Verified by building a copy with no sibling checkouts. Not yet verified by a real
      `git clone`, which needs this work committed.
- [x] `verify` is opt-in on the umbrella crate (removed from `default`).
- [x] `MPL-2.0` removed from the `deny.toml` allow-list (unused).
- [x] `no_std` CI covers every crate that declares `#![no_std]`, plus a wasm32 build. `vibro` (std-only) and
      `verify` (proptest) are out of scope per spec.txt.
- [ ] CI jobs not yet run on GitHub: macOS, MSRV, coverage, semver, bench artifact, wasm, release-plz.

### 4.2 Correctness bugs
- [x] **STC contour** (`tpt-acou-vibro/src/transmission.rs`): ASTM E413 shape, analytic tests.
- [ ] **STC worked example** from the published standard. *(blocked: paywalled)* Partial: the rating rule takes any
      user contour (`tpt_acou_vibro::transmission::RatingContour`, `rating_against_contour`; partition study
      `[partition.rating_contour]` with a required `source`). The built-in E413-shaped values are not checked against
      the standard, and the reports say so.
- [x] **ITU-R 468 table** (`tpt_acou_psycho::weighting::ITU468_TABLE_DB`): all 21 Table 1 response values match the
      official BS.468-4 PDF (free from ITU), checked on the rendered page. The 31.5 kHz row was missing and is added. The
      tolerance column is reproduced from the mask points by log-frequency interpolation
      (`tpt_acou_psycho::weighting::itu468_tolerance_db`, within 0.05 dB of all 20 ranged rows). The 31.5 kHz asymmetric
      tolerance is not represented.
- [x] **Quasi-peak detector, tone-burst limits only** (`tpt_acou_psycho::quasi_peak::QuasiPeakDetector`): a two-channel
      peak-hold model whose parameters were found by search so that the readings meet BS.468-4 Tables 2 and 3 (single and
      repeated 5 kHz bursts). The tests drive the sampled detector through the limits, and the smallest margin is 1.2
      points (50 ms burst). The standard gives no circuit or time constants, so the internal structure is ours. The other
      §2 requirements (overload, reversibility, overswing, calibration, input impedance) are not tested or claimed.
- [x] **Moore–Glasberg loudness** (`tpt-acou-psycho/src/loudness_moore.rs`): implemented from the published
      model structure (ERB filters, absolute threshold, compressive specific loudness), anchored so that 1 kHz at
      40 dB is 1 sone. Tested for the anchor, monotonicity, masking and the ERB scale. **Not validated** against
      ISO 532-2 reference values (licensed). Outer/middle-ear and temporal effects are not modelled.
- [x] `zwicker_stationary`: spectral spreading added. 1 kHz at 40 dB gives 1.02 sone (test pins it to 5%).
- [x] **Ray tracer**: particles keep tracing after capture; partition `injected = absorbed + remaining`; RT60 test
      at 10%. Per-ray seeding; the `parallel` feature gives identical output.
- [x] Ray tracer input validation (`RayTracerConfig::validate`, positions inside the room).
- [x] **FW-H** (`tpt-acou-aero/src/fwh.rs`): `1/(1−M_r)` factor; `omega = 0` returns 0; NaN where the far field is
      undefined; central differences with a step scaled to the blade-passing frequency. The previous step left the
      BPF harmonic about 10% low. **Validated** against an independent retarded-time quadrature of the thickness
      integral (`ring_bpf_matches_the_retarded_time_quadrature`, 2%). The Garrick–Watkins reference (NACA TN 3018)
      covers rotor *loading*, not thickness, so it is not the right reference for this model.
- [x] `b_weighting_db` doc matches the code (+1.96 dB).
- [x] `ShoeboxRoom::echogram` takes a `&Medium`.
- [x] `SoundPressureLevel::try_from_*` constructors.
- [x] Umbrella crate doc example cleaned up.
- [x] CFL Kani harness scope documented: a classification lemma, not a stability proof (no time-domain solver).

### 4.3 Missing features vs spec
- [x] **GMRES** (`tpt_acou_wave::krylov`): restarted, right-preconditioned, complex. Matches the dense BEM solve.
- [x] **FMM-backed Neumann solve**: the adjoint double-layer operator `K′` applied through the single-level FMM
      (field-side normal derivative of the expansion, 6-radius far criterion), and `solve_exterior_neumann_fmm`.
      Validated against the dense operator and the dense solve, both to a 1e-3 relative bound (`fmm::tests`,
      feature `fmm`). Interior Neumann and Dirichlet solves are not wired to the FMM.
- [x] **ISM with complex impedance** (`ShoeboxRoom::transfer_function_impedance`): locally reacting oblique
      reflection from the unfolded image path. Validated: stiff walls reproduce the rigid room. A matched wall is not
      perfectly absorbing at oblique incidence (stated in the docs).
- [ ] **Multi-parameter gradients of ray and ISM energy** (per-band absorption, source and receiver position). Done:
      ISM, all 15 partials (`tpt_acou_ray::differentiable::ism_energy_gradient`, checked against finite differences). Done:
      ray absorption, exact for fixed seeds (`tpt_acou_ray::ray::windowed_energy_absorption_gradient`, checked against
      finite differences of the full tracer to 1e-4). Still open: ray gradients with respect to source and receiver
      position, which change the capture events, so they are not differentiable. Scattering gradients are not done.
- [x] **Room impulse response** from the echogram (`rir_from_echogram`); **hybrid ISM–ray RIR** (`hybrid_rir`, tail
      calibrated to the ISM energy at the crossover, gain returned); **broadband hybrid** (modal below the Schroeder
      frequency, image sources above: `broadband::hybrid_broadband_ir`, checked against the two transfer functions in
      their bands).
- [x] **Three-way hybrid** (ISM early, ray late, modal low frequency in one response): `tpt_acou_ray::broadband::three_way_ir`.
      The signed broadband response (modal below f_S, ISM above) is kept up to the crossover, and a ray tail calibrated to
      its energy follows. The tail is broadband, so modal low-frequency decay past the crossover is lost. Not validated
      against a measured or wave-solved reference.
- [x] **Rigid-box modal transfer** (`modal::RigidBoxModes`): reciprocity exact; resonance location and damping
      checked. The image sum is not used as a reference, because it is only conditionally convergent for lossless
      rigid rooms.
- [x] **serde** (`tpt-acou-core`, `serde` feature; JSON round-trip test). `Medium` and `Solid` serialise only, so
      validation cannot be bypassed on read.
- [x] **rayon** (`tpt-acou-ray`, `parallel` feature). Per-ray seeds make the output independent of threads.
- [x] **criterion benchmarks** (`crates/tpt-acou-ray/benches/ray.rs`), and a CI job that keeps the report.
- [x] **Benchmark regression gate**: `bench-gate` job in `.github/workflows/ci.yml` benchmarks the PR base and head on one
      runner (criterion `--save-baseline` / `--baseline`), so no stored baseline is needed. `scripts/bench_gate.py` fails a
      benchmark only when its mean slowdown exceeds 15 % and the 95 % lower bound is above zero. The gate logic is tested
      (`--self-test`), and the criterion flow was run locally. Not run on GitHub yet.

### 4.4 New features (innovation)
- [x] **ISO 3382 metrics** from an RIR: T20, T30, EDT, C50, C80, D50 (closed-form tests), and **sound strength G**
      (free-field reference at 10 m; tested at 0 dB and +6.02 dB).
- [x] **Speech Transmission Index** (`tpt_acou_ray::sti`): octave bands, 14 modulation frequencies, IEC 60268-16
      structure. Validated against the closed-form index for an exponential noise tail (within 0.03). The male and
      female weights come from a secondary source; both sets sum to exactly 1 with their redundancy terms (checked).
      Not confirmed against the IEC edition. Octave filtering is rectangular in frequency, not the standard's
      band-pass filters. *(blocked: IEC table check)*
- [x] **RIR export** to float WAV, and **auralization** by FFT convolution (`rir::auralize`, CLI `auralize`).
- [x] **Differentiable inverse design** (`design`): closed-form uniform absorption for a target T60, layout-preserving
      scaling, analytic gradients. Not done: fitting a target C80, or per-panel optimisation against ray or ISM energy.
- [x] **Uncertainty quantification** (`design::rt60_uncertainty`, seeded Monte Carlo over per-wall absorption; CLI
      `uq`). Covers absorption uncertainty only, not geometry or source uncertainty.
- [x] **STL and OBJ import** (`tpt_acou_ray::stl`, `tpt_acou_ray::obj`). One material per scene; per-`usemtl`
      material assignment is not done.
- [x] **Octave-band absorption database** (`tpt_acou_ray::database`): three finishes, each taken from one named
      vendor chart, with the chart named in the data. Not a complete ISO 354 database. *(blocked: a cited, complete
      data source)*
- [x] **Standards compliance report**: the `compliance` study kind reports ISO 3382-1, IEC 60268-16 and ASTM E413
      results, and lists what it does not report (`run_compliance` in `crates/tpt-acou-cli/src/study.rs`, tested). The
      two ratings it does not compute, OITC (ASTM E1332 reference table, *blocked*) and Lw (needs measured source data,
      ISO 3741), stay open under Blockers.
- [x] **Python bindings** (`bindings/python`, PyO3): `cargo check` passes. Not built into a wheel or imported here
      (needs a Python environment and `maturin`).
- [x] **WASM bindings** (`crates/tpt-acou-wasm`): native tests pass, and the wasm32 build succeeds.

### 4.5 Adoption: examples, docs, templates
- [x] **Examples**: `shoebox_rt60`, `room_impulse_response`, `partition_stc`, `sphere_radiation_bem`, `fan_noise_fwh`,
      `optimize_room_gradient`, `loudness_a_weighting`, `spatial_audio_vbap`. Each asserts against a reference.
- [x] **Prelude, `RoomBuilder`, `Simulation`** (echogram, impulse response, ray-traced T30, STI, G).
      [x] A `Simulation` path to a wave solve: `Simulation::broadband_ir` (modal solution below the Schroeder frequency,
      image sources above) and `Simulation::three_way_impulse_response`. The modal solve is the rigid-box modal sum, not a BEM or FEM solve.
- [x] **README**: quick start, examples table, does/doesn't table, validation status, CI and licence badges.
- [x] **Cookbook** (`docs/cookbook.md`) and **validation page** (`docs/validation.md`).
- [x] **mdBook layout** (`book.toml`, `docs/SUMMARY.md`). Not built here (mdbook not installed).
- [x] **Templates**: room, partition, rotor and compliance studies in `templates/`.
- [x] **CLI** `tpt-acou`: `rt60`, `stc`, `weight`, `study`, `design`, `uq`, `auralize`.

### 4.6 Automation & hygiene
- [x] **Dependabot** (cargo and GitHub Actions).
- [x] **Semver check** on pull requests (`cargo-semver-checks-action`). Needs a baseline; not run on GitHub.
- [x] **release-plz** (`release-plz.toml`, workflow, `CHANGELOG.md`). `publish = false`, because the git dependencies
      cannot be published to crates.io. Not run on GitHub.
- [x] **Coverage** (`cargo llvm-cov --summary-only`), **MSRV** (1.85, libraries), **Windows and macOS** test matrix.
      macOS and MSRV not yet run on GitHub.
- [x] **Symbol check**: `crates/tpt-acoustics/tests/todo_symbols.rs` fails if a ticked item names code that does not
      exist, or a file that is missing.
- [x] **Issue and PR templates**.

## Phase 5 — Open, project-defined ratings (no standard text or tables used)

> Our own definitions, documented in `docs/open-rating.md`. They are not STC, OITC or any other standard, and must
> not be reported as such. Validation uses physics and closed-form cases only.

- [x] **TPT Isolation Index v1** (`tpt_acou_vibro::transmission::tii_rating`): a log-frequency reference line at 5 dB per
      octave, with the deficit rule. Tested against a line spectrum (rating = level + 2), against the diffuse mass law
      (rating step 5 dB per doubling of mass, observed), against a 3 dB shift, and against missing bands. Shown in the
      partition study, labelled as a project definition.
- [x] **User-supplied contours** for the deficit rule (`RatingContour`, `rating_against_contour`), so that a licensed
      table can be used by the person who holds it. The source is required and quoted in every report.
- [x] **Open rating v2** (`tpt_acou_vibro::transmission::tii_v2_rating`): the same line contour on 20 bands (63 Hz to
      5 kHz), with the total limit at 2 dB per band. Decision: v1 stays the headline value, and v2 is a comparison that
      shows low-frequency sensitivity. Tested against the line spectrum and the mass law, and shown in the partition
      report.

### Blockers (need a source the repository does not have)
- ~~**BS.468-4 printed table**~~ resolved: the free ITU PDF was checked (rendered). The **quasi-peak** detector still needs
  quasi-peak burst limits are met by a fitted model (see the quasi-peak item in Phase 4); the standard gives no circuit.
- **ASTM E1332 reference contour** (OITC): only 3 of 16 levels were found.
- **ASTM E413 worked example** for the STC check. The user-contour input is in place, so a copied table can be checked here.
- **ISO 532-2 reference values** (licensed) for Moore–Glasberg.
- **IEC 60268-16 weights**: confirm the edition and the male and female tables. Web search (round 3) found that
  IEC 60268-16:2020 removed the female set and changed the male test spectrum, so the female weights in `sti.rs` are
  from the 2003 edition and unverified. The male table in `sti.rs` is still a secondary source. Needs the standard itself.
- **Complete octave-band absorption data** (ISO 354 style), with a citable source for each value.
- **A published rotor *thickness* benchmark**: the Garrick–Watkins results found cover loading, not thickness.
