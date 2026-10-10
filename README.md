# tpt-acoustics

A pure-Rust, AI-native computational acoustics library: geometric room acoustics,
wave acoustics (BEM/FEM for the Helmholtz equation), aeroacoustics, vibroacoustics,
and psychoacoustic modeling — **no C/C++ FFI, no external solver dependencies**.

Dual-licensed **MIT OR Apache-2.0**. © TPT Solutions.
See [spec.txt](spec.txt) for the full design and [todo.md](todo.md) for the phased
build tracker.

## Quick start

Add the umbrella crate (from GitHub, since the `tpt-*` sibling dependencies are
pinned git revisions):

```toml
[dependencies]
tpt-acoustics = { git = "https://github.com/tpt-solutions/tpt-acoustics" }
```

Then:

```rust
use tpt_acoustics::prelude::*;

fn main() -> Result<(), AcousticError> {
    let room = ShoeboxRoom::uniform([8.0, 5.0, 3.0], BoundaryMaterial::concrete_painted())?;
    let air = Medium::air_20c();
    let echogram = room.echogram(Vec3::new(1.5, 1.2, 1.4), Vec3::new(6.2, 3.7, 1.6), 1000.0, 12, &air)?;
    let rir = rir_from_echogram(&echogram, 48_000.0, 1)?;
    println!("C80 = {:.1} dB", clarity_db(&rir, 48_000.0, 80.0).expect("late energy"));
    Ok(())
}
```

Runnable examples live in `crates/tpt-acoustics/examples/`. Each one asserts
against an analytic or literature value:

| Example | Shows | Reference |
| --- | --- | --- |
| `shoebox_rt60` | Ray-traced T30 of a shoebox | Sabine, within 10 % |
| `room_impulse_response` | Echogram → RIR → C50/C80/D50, WAV export | Ordering and bounds |
| `partition_stc` | Diffuse mass law and STC of a leaf | 6 dB/octave mass law |
| `sphere_radiation_bem` | Exterior BEM, pulsating sphere | Closed-form solution, within 5 % |
| `fan_noise_fwh` | Blade-passing tone of a rotor | BPF dominates its half-harmonic |
| `optimize_room_gradient` | Analytic gradient of room energy, Newton refinement | Finite differences; recovered length |
| `loudness_a_weighting` | A-weighting and phon/sone | IEC 61672 nominal points |
| `spatial_audio_vbap` | VBAP gains on a quad ring | Pulkki panning laws |

## What it does and does not do

| Does | Does not (yet) |
| --- | --- |
| Image-source early reflections for shoebox rooms | Wall phase or complex impedance in the ISM |
| Monte Carlo ray tracing with BVH scenes | Hybrid ISM/ray/wave RIR (ray tail not merged) |
| ISO 3382 T20/T30/EDT, C50/C80, D50 | STI or G |
| Helmholtz BEM (exterior/interior), single-level FMM, restarted GMRES with Jacobi preconditioning | GMRES driven by the FMM operator (the solver takes any matvec; the BEM wiring is open) |
| Cavity FEM modes | Time-domain solvers |
| FW-H monopole/dipole and rotating-ring thickness noise | Lighthill stress-tensor sources from CFD |
| Mass-law TL, ASTM E413 STC | Measured-spectrum STC corrections or OITC |
| A/B/C/ITU-R 468 weightings, VBAP, HRTF, Ambisonics | Moore loudness (not implemented) |

## Validation status

Per crate, what is checked against a reference and what is not. "Analytic"
means a closed-form result; "property" means a `proptest` invariant.

| Crate | Checked against | Open |
| --- | --- | --- |
| `tpt-acou-core` | Unit conversions, SPL round trips, impedance | — |
| `tpt-acou-ray` | Sabine T30 (10 %), ISO 3382 closed forms, energy partition (property) | Ray tail not calibrated against ISM levels |
| `tpt-acou-wave` | Pulsating sphere (BEM, 5 %), cavity FEM modes, FMM vs dense, GMRES vs dense BEM solve | GMRES on the FMM operator |
| `tpt-acou-vibro` | Mass law, ASTM E413 contour on analytic spectra, beam first mode vs analytic (5 %), FEM→modal→BEM partition TL | Published STC worked example |
| `tpt-acou-aero` | Monopole scaling, BPF tone, dipole directivity | Published rotor-noise benchmark |
| `tpt-acou-psycho` | IEC 61672 A-weighting, phon/sone anchor | ITU-R 468 table only partly confirmed; loudness is a simplified Zwicker model (~0.57 sone at 1 kHz/40 dB, not ISO 532-1) |
| `tpt-acou-verify` | Proptest invariants (reciprocity, energy, CFL), Kani harnesses | Kani runs in CI on Linux only |

The `todo.md` tracker is the authoritative list of open work. Its Phase 4 is
still in progress.

## Design principles

- **Zero external solver dependencies.** Every algorithm (Image Source Method,
  Helmholtz BEM, Ffowcs Williams–Hawkings) is implemented from scratch in pure
  Rust, on top of the in-house `tpt-math` numeric substrate.
- **Strict licensing.** Workspace crates are `MIT OR Apache-2.0`; no
  Apache-2.0-only dependency may enter the graph (enforced via `deny.toml` in CI).
- **Type-state & dimensional safety.** Phantom types distinguish acoustic
  pressure, particle velocity and power at compile time; sound pressure levels
  are handled strictly in logarithmic decibel scales with the proper 20 µPa
  reference.
- **Native formal verification.** Kani bounded model checking for critical
  numerical loops (ISM reflection orders, CFL stability) and property-based
  testing (`proptest`) for physical laws (reciprocity, energy conservation).
- **AI-native.** Differentiable hooks (forward-mode duals over geometry) let
  agents perform inverse room design and neural acoustic rendering.

## Workspace layout

| Crate | Purpose |
| --- | --- |
| `tpt-acou-core` | Foundation: wave equations, media, unit-safe types (`no_std + alloc`) |
| `tpt-acou-ray` | Geometric acoustics: ISM, Monte Carlo ray tracing, BVH scene intersection |
| `tpt-acou-wave` | Wave acoustics: frequency-domain BEM/FEM for the Helmholtz equation |
| `tpt-acou-vibro` | Vibroacoustics: structural-acoustic coupling, radiation efficiency, TL/STC |
| `tpt-acou-aero` | Aeroacoustics: Lighthill analogy, Curle, Ffowcs Williams–Hawkings |
| `tpt-acou-psycho` | Psychoacoustics: weighting filters, loudness, HRTFs, Ambisonics, VBAP |
| `tpt-acou-verify` | Verification harnesses: Kani proofs, proptest strategies, invariant checks |
| `tpt-acoustics` | Feature-gated umbrella crate re-exporting all of the above |

## Non-goals

- No C/C++/Fortran FFI, no external PDE solver processes, no GPU-required paths.
- Not a real-time audio engine: `tpt-acou-psycho` models perception, it does not
  mix or render an audio graph (that is `tpt-audio`'s job).
- Not a CFD solver: aeroacoustic source data is consumed, not produced (that is
  `tpt-physics`' job).

## Building & testing

The build fetches the `tpt-math`, `tpt-fem`, `tpt-dsp` and `tpt-materials`
crates from GitHub at pinned revisions (see `[workspace.dependencies]` in
`Cargo.toml`), so it needs network access the first time.

`core`, `ray`, `wave`, `aero` and `psycho` build `no_std`; `vibro` and `verify`
are std-linked (FEM and proptest tooling). `core`'s optional `materials`
feature bridges to `tpt-mat-database` for standard structural-material
properties.

```sh
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --check
cargo deny check
```

`no_std` targets (`tpt-acou-core`, …) are built for `thumbv6m-none-eabi` in CI:

```sh
cargo build -p tpt-acou-core --target thumbv6m-none-eabi
```

## License

Licensed under either of [MIT](LICENSE-MIT) or [Apache-2.0](LICENSE-APACHE) at
your option.
