# tpt-acoustics

A pure-Rust, AI-native computational acoustics library: geometric room acoustics,
wave acoustics (BEM/FEM for the Helmholtz equation), aeroacoustics, vibroacoustics,
and psychoacoustic modeling — **no C/C++ FFI, no external solver dependencies**.

Dual-licensed **MIT OR Apache-2.0**. © TPT Solutions.
See [spec.txt](spec.txt) for the full design and [todo.md](todo.md) for the phased
build tracker.

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

## Status

v1.0: all four phases of [todo.md](todo.md) are implemented — geometric
acoustics, wave acoustics (BEM/FEM), vibroacoustics, aeroacoustics and
psychoacoustics, with 150+ unit/property/integration tests, Kani proof
harnesses (CI on Linux), `no_std` builds for the solver crates, and
`cargo deny` license enforcement. `core`, `ray`, `wave`, `aero` and
`psycho` build `no_std` (psycho's DSP substrate is `tpt-dsp-core`);
`vibro`/`verify` are std-linked (FEM/proptest tooling), and `core`'s
optional `materials` feature (off for `no_std` builds) bridges to
`tpt-mat-database` for standard structural-material properties.

## Building & testing

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
