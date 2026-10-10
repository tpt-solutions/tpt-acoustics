# Cookbook

A single-page guide to choosing a solver, the units and conventions used
throughout, and a few recipes. Every recipe has a runnable counterpart in
`crates/tpt-acoustics/examples/` or the `tpt-acou` command line.

## Choosing a solver

| Question | Use | Crate | Valid range |
| --- | --- | --- | --- |
| Reverberation time of a box room | Sabine or Eyring formula | `tpt-acou-ray` (`ShoeboxRoom`) | Diffuse field; `α` well below 1 |
| Early reflections of a box room | Image Source Method | `tpt-acou-ray` (`echogram`) | Shoebox geometry; wall phase ignored |
| Late reverberation, arbitrary geometry | Monte Carlo ray tracing | `tpt-acou-ray` (`trace_shoebox`, `trace_scene`) | High frequency; energy only |
| Radiation from a closed surface | Helmholtz BEM | `tpt-acou-wave` (`HelmholtzBem`) | Exterior/interior, away from interior resonances |
| Large BEM systems | GMRES over a matvec | `tpt-acou-wave` (`krylov::gmres`) | Any operator closure; FMM wiring still open |
| Cavity modes of a box | Cavity FEM | `tpt-acou-wave` (`CavityFem`, feature `fem`) | Rigid walls |
| Partition transmission loss | Mass law, STC | `tpt-acou-vibro` (`transmission`) | Mass-controlled range |
| Rotor or fan tones | FW-H rotating ring | `tpt-acou-aero` (`fwh`) | Subsonic, compact source |
| Loudness, weightings | A/B/C, ITU-R 468, Zwicker | `tpt-acou-psycho` | See the validation table in the README |

## Units and conventions

- **Length** in metres, **time** in seconds, **frequency** in hertz.
- **Pressure** in pascals; sound pressure levels in dB re 20 µPa.
- **Intensity** in W/m² (re 1 pW/m²); **power** in W (re 1 pW).
- **Surface mass** `m` in kg/m²; **absorption** `α` and **scattering** `s`
  are dimensionless in `[0, 1]`.
- **Time convention** is `e^{+iωt}`. Outgoing waves carry `e^{−ikr}`.
- **Speed of sound** defaults to 343 m/s in `SPEED_OF_SOUND_REFERENCE`. Use
  `Medium::air(θ)` for temperature-dependent air.
- **Echogram weights** are energies (`rc²/(4πd)²` for a unit source), not
  pressures. `rir_from_echogram` uses random signs, so its waveform is
  energy-correct but not phase-correct.
- **Ray tracer energy** is normalised to 1 injected unit; `detected_*` is a
  tally of receiver passes, separate from the energy partition.

## Glossary

- **Absorption coefficient α**: fraction of incident energy a wall absorbs.
- **BEM**: boundary element method. Solves on the surface, not the volume.
- **C50 / C80 (clarity)**: early-to-late energy ratio at 50 ms (speech) or
  80 ms (music), in dB.
- **D50 (definition)**: early energy fraction within 50 ms.
- **EDT**: early decay time, from the 0 to −10 dB slope, extrapolated to 60 dB.
- **FW-H**: Ffowcs Williams–Hawkings equation for moving sources.
- **ISM**: image source method.
- **RIR**: room impulse response.
- **Schroeder frequency** `f_S = 2000 √(T/V)`: above it the modal density is
  high enough for statistical methods.
- **STC**: sound transmission class, a single-number rating from the
  ASTM E413 contour.
- **T20 / T30**: decay times fitted over the −5 to −25 dB or −5 to −35 dB
  windows, extrapolated to 60 dB.

## Recipes

### Reverberation time of a room

```sh
tpt-acou rt60 10 8 6 0.15 1000
```

Prints Sabine and Eyring at 1 kHz. For a full study with a ray-traced T30
check, see the next recipe.

### A room study from a file

```sh
tpt-acou study templates/room-study.toml
```

The report gives Sabine and Eyring T60, a ray-traced T30 with its difference
from Sabine, the Schroeder frequency, and C80 and D50 from the image-source
response. Copy `templates/room-study.toml` and edit the room, positions and
analysis settings. Unknown keys are rejected.

### Partition and rotor studies from files

```sh
tpt-acou study templates/partition-study.toml   # mass-law TL at the STC bands, STC
tpt-acou study templates/rotor-study.toml       # FW-H blade-passing amplitude
```

Each file names its `kind`. The rotor study rejects a blade tip at or above
Mach 1, where the far-field model is undefined, and it reports a zero
half-harmonic as zero rather than as a ratio of noise.

### Hybrid impulse response

`tpt_acou_ray::rir::hybrid_rir` joins an image-source early part to a
ray-traced late part. The ray tail is scaled so its energy in a window around
the crossover matches the echogram. The gain is returned for inspection.

### Early reflections and clarity in code

```rust
use tpt_acoustics::prelude::*;

let room = ShoeboxRoom::uniform([8.0, 5.0, 3.0], BoundaryMaterial::concrete_painted())?;
let air = Medium::air_20c();
let echogram = room.echogram(Vec3::new(1.5, 1.2, 1.4), Vec3::new(6.2, 3.7, 1.6), 1000.0, 12, &air)?;
let rir = rir_from_echogram(&echogram, 48_000.0, 1)?;
let c80 = clarity_db(&rir, 48_000.0, 80.0);
```

Clarity and definition need late energy. The image-source echogram stops at
its highest arrival, so late-decay metrics such as T30 must come from the ray
tracer.

### Transmission loss and STC of a partition

```sh
tpt-acou stc mass 5
```

For a measured spectrum, pass the 16 one-third-octave values from 125 Hz to
4 kHz: `tpt-acou stc tl 34 37 40 ...`.

### Loading a mesh for ray tracing

```rust
let scene = tpt_acou_ray::obj::scene_from_obj(&std::fs::read_to_string("room.obj")?, material)?;
let scene = tpt_acou_ray::stl::scene_from_stl(&std::fs::read("room.stl")?, material)?;
```

Every face gets the same material. Per-material assignment from `usemtl`
names is not yet available.

### Schroeder frequency

```rust
use tpt_acou_ray::iso3382::schroeder_frequency_hz;

let f_s = schroeder_frequency_hz(0.5, 200.0); // Some(100.0) Hz
```

Below `f_S`, use wave methods (BEM, FEM) or modal analysis; above it, the ray
and diffuse-field methods apply.

### Fitting a room to a target energy

See `examples/optimize_room_gradient.rs`. The image-source energy oscillates
with room length on the scale of the wavelength, so refinement only converges
near the starting point. Scan the length first, then refine.

## Not covered yet

- A mdBook layout. This page is the single-file version.
- Octave-band absorption databases (ISO 354 style).
- Hybrid ISM/ray/wave impulse responses; see `todo.md`.
