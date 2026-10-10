# Validation

What each result is checked against, the tolerance the test enforces, and
where the check lives. The README has the short version. This page is the
detail. A row with no tolerance is a qualitative or structural check.

"Analytic" means a closed-form result. "Literature" means a published
formula or table. "Property" means a randomised invariant (`proptest`).
"Kani" means a bounded model-checking proof, run in CI on Linux only.

## Geometric acoustics (`tpt-acou-ray`)

| Result | Reference | Tolerance | Test |
| --- | --- | --- | --- |
| Ray-traced T30 of a 10 × 8 × 6 m room | Sabine (literature) | 10 % | `ray::rt60_tracks_sabine_within_tolerance` |
| Energy partition, lossless room | Analytic | 1e-9 | `ray::energy_partition_is_exact_in_lossless_rooms` |
| Energy partition, absorbing and scattering walls | Analytic | 1e-9 | `ray::energy_partition_holds_with_absorption_and_scattering` |
| Energy partition holds for random rooms | Property | 1e-9 | `verify::properties` (lossless and lossy tracing) |
| T20, T30, EDT of an exponential decay | Analytic | 2 % | `iso3382::decay_times_recover_exact_exponential_t60` |
| C50 of an exponential decay | Closed form `10·log10(e^{kt} − 1)` | 0.05 dB | `iso3382::clarity_matches_closed_form_for_exponential_decay` |
| D50 of an exponential decay | Closed form `1 − e^{−kt}` | 1e-3 | `iso3382::definition_matches_closed_form_for_exponential_decay` |
| Direct arrival time and weight | Analytic | 1e-12 | `verify::properties::echogram_direct_arrival_matches_geometry` |
| Image-source count bounded by `(2N+1)³` | Analytic | exact | `verify::proofs::ism_count_never_expects_order_bound` (Kani) |
| Scene tracer agrees with shoebox tracer | Internal consistency | 1e-9 | `ray::scene_tracer_matches_shoebox_partition` |
| Schroeder frequency `2000·√(T/V)` | Worked example | 1e-9 | `iso3382::schroeder_frequency_worked_example` |
| OBJ and STL cubes, tetrahedron, ray hit | Analytic geometry | 1e-6 | `obj::tests`, `stl::tests` |

## Wave acoustics (`tpt-acou-wave`)

| Result | Reference | Tolerance | Test |
| --- | --- | --- | --- |
| Pulsating sphere, exterior Neumann | Closed form `p(r)` | 5 % | `bem::exterior_pulsating_sphere_neumann` |
| Pulsating sphere, exterior Dirichlet | Closed form `p(r)` | 5 % | `bem::exterior_sphere_dirichlet` |
| Interior sphere with constant boundary value | Closed form `j₀(kr)/j₀(ka)` | 0.05 | `bem::interior_sphere_dirichlet` |
| Radiated power of a pulsating sphere | Closed form | see test | `bem::radiated_power_of_pulsating_sphere` |
| Dense and FMM matvecs agree | Internal consistency | see test | `fmm::fmm_matvec_matches_dense` |
| GMRES matches the dense BEM solve | Internal consistency | 1e-7 | `krylov::matches_the_dense_bem_exterior_solve` |
| GMRES restarts, preconditioning, non-convergence | Analytic test systems | 1e-8 | `krylov::tests` |
| Cavity FEM modes | Analytic rigid-cavity frequencies | see test | `fem_cavity::fundamental_mode_of_box_cavity`, `analytic_frequencies_match_textbook` |
| BEM kernel reciprocity; radiation efficiency is passive (σ ≤ 1) | Analytic symmetry and passivity | exact | `verify::properties::bem_single_layer_kernel_reciprocity`, `bem_radiation_efficiency_is_passive` |

## Vibroacoustics (`tpt-acou-vibro`)

| Result | Reference | Tolerance | Test |
| --- | --- | --- | --- |
| Diffuse mass law: doubling of mass gains 6 dB | Mass law | 1 dB | `transmission::mass_law_doubles_per_surface_mass` |
| STC of the reference contour | ASTM E413 contour | exact (52) | `transmission::stc_of_the_reference_contour_is_52` |
| STC of a flat 60 dB spectrum | ASTM E413 criteria | exact (60) | `transmission::stc_of_uniform_60_db_spectrum` |
| ASTM E413 contour shape | Standard | exact | `transmission::stc_contour_matches_astm_e413_shape` |
| Beam first mode of a partition | Analytic Euler–Bernoulli | 5 % | `tests/integration_partition.rs` |
| TL of a partition, FEM → modal → BEM | Mass law at 500 Hz | 8 dB | `tests/integration_partition.rs` |
| STC of the partition vs mass law | Internal consistency | 8 | `tests/integration_partition.rs` |

A published STC worked example is not yet in the test suite. The contour
checks use the ASTM E413 shape and analytic spectra.

## Aeroacoustics (`tpt-acou-aero`)

| Result | Reference | Tolerance | Test |
| --- | --- | --- | --- |
| Stationary monopole peak pressure | Analytic `Q₀ω / (4π a c₀)` | 5 % | `fwh::stationary_monopole_matches_curle_like_scaling` |
| Approaching source louder than receding | Qualitative (Doppler) | 10 % | `fwh::moving_monopole_is_amplified_upstream` |
| Rotating ring tone at the blade-passing frequency | Structure (BPF harmonic) | 4× the off-tone | `fwh::rotating_ring_has_bpf_tone` |
| Stationary ring radiates nothing | Analytic | exact (0) | `fwh::stationary_ring_radiates_nothing` |
| Dipole directivity | Qualitative | 0.15 ratio | `fwh::dipole_directivity_of_rotating_loading` |

No published rotor-noise benchmark is wired in yet. Amplitude checks against
a reference are open; see `todo.md`.

## Psychoacoustics (`tpt-acou-psycho`)

| Result | Reference | Tolerance | Test |
| --- | --- | --- | --- |
| A-weighting at 100 Hz, 1 kHz, 10 kHz | IEC 61672 nominal (−19.1, 0, −2.5 dB) | 0.1 dB | `weighting::a_weighting_spot_values`, `examples/loudness_a_weighting` |
| ITU-R 468 peak at 6.3 kHz | BS.468-4 (secondary sources) | exact | `weighting::itu468_table_values_and_interpolation` |
| Phon and sone anchor, doubling | ISO 532-1 definition | 1e-12 | `loudness::phon_sone_anchor_and_doubling` |
| 1 kHz / 40 dB tone loudness | ISO 532-1 (1 sone) | 5 % | `loudness::zwicker_one_khz_tone_near_one_sone` |
| Spreading slopes | Zwicker's published slopes | 1e-9 | `loudness::spreading_slopes_match_zwicker` |

Open: the ITU-R 468 table is only partly confirmed. The 3.15–10 kHz points
come from secondary sources, not the standard. The loudness model is a
simplified Zwicker approximation, not a certified ISO 532-1 implementation,
and Moore–Glasberg is not implemented.

## Core (`tpt-acou-core`)

| Result | Reference | Tolerance | Test |
| --- | --- | --- | --- |
| SPL round trip through pressure, intensity, power | Definition | 1e-9 | `quantity::spl_round_trips_through_pressure` |
| SPL `try_` constructors reject non-positive input | Definition | exact | `quantity::spl_try_variants_reject_non_positive_input` |

## Numerical proofs (`tpt-acou-verify`)

| Property | Method | Where it runs |
| --- | --- | --- |
| Image-source order bound | Kani | CI, Linux |
| CFL classification matches the rational bound | Kani | CI, Linux |
| Grazing-incidence finiteness | Kani | CI, Linux |
| Reciprocity `H_AB = H_BA` | Property | `cargo test` |
| Energy conservation in lossless rays and BEM | Property | `cargo test` |

The CFL harness checks a float-versus-integer comparison. It does not exercise
a time-domain solver, because the workspace has none.

## Reproducing a check

```sh
cargo test -p tpt-acou-ray iso3382          # ISO 3382 closed forms
cargo test -p tpt-acou-vibro transmission   # STC and mass law
cargo run --release -p tpt-acoustics --example sphere_radiation_bem
```

Each example asserts its reference and exits non-zero on failure.
