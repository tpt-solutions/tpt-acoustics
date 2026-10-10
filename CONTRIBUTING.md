# Contributing to tpt-acoustics

Thanks for your interest in contributing to `tpt-acoustics`! This repository is
part of the TPT Solutions `tpt-*` family and follows its house rules.

## Ground rules

1. **Pure Rust only.** No C/C++ FFI, no external solver processes. Algorithms
   are implemented from scratch on top of the in-house `tpt-math` substrate.
2. **Licensing.** Every dependency must be MIT, Apache-2.0, or a dual
   MIT/Apache-2.0 license — *never* Apache-2.0-only. `cargo deny check` must
   stay clean (`deny.toml` enforces this; CI runs it on every PR).
3. **Physical correctness is testable.** New physics must ship with either
   (a) a comparison against an analytic solution, (b) a property test of a
   physical law (reciprocity, energy conservation, passivity), or (c) a Kani
   proof of a safety invariant. Numerical tolerances must be justified in the
   test name or comment.
4. **Unit safety.** Use the phantom-typed quantities from `tpt-acou-core`
   (`Pressure`, `Intensity`, `SoundPressureLevel`, …) at API boundaries. Raw
   `f64` is fine inside a numerical kernel once units are documented.
5. **`no_std` where the spec asks.** `tpt-acou-core` must stay `no_std + alloc`
   compatible. Don't reach for `std`; use `libm` for transcendentals.

## Checklist for every crate change

- [ ] `cargo fmt` clean
- [ ] `cargo clippy --workspace --all-targets -- -D warnings` clean
- [ ] `cargo test --workspace` green (including new property tests)
- [ ] `cargo deny check` clean
- [ ] `rustdoc` on all public items (`cargo doc --no-deps` warning-free)
- [ ] `no_std` build still green for the crates that promise it
- [ ] Update [todo.md](todo.md) checkboxes when you complete a tracked item

## Commit style

Short imperative subject lines (`ism: guard against zero-order shoebox`),
body paragraphs explaining *why* when the diff doesn't make it obvious.

## Reporting bugs

Open an issue with the smallest reproduction (a failing test is ideal), the
expected physical behaviour with a reference, and the observed value.
