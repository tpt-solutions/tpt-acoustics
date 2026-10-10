//! Studies described in TOML.
//!
//! A study file names its `kind` and carries that kind's settings. Three
//! kinds are supported:
//!
//! * `room`: a shoebox room with a source and receiver. The report gives
//!   Sabine and Eyring T60, a ray-traced T30 checked against Sabine, the
//!   Schroeder frequency, and clarity and definition from the image-source
//!   response.
//! * `partition`: a single leaf of given surface mass. The report gives the
//!   diffuse-field mass-law transmission loss at the 16 STC bands and the STC.
//! * `rotor`: a rotor of `B` blades on a ring. The report gives the FW-H
//!   thickness-noise amplitude at the blade-passing frequency and its
//!   half-harmonic, from a sampled pressure history.
//!
//! Templates for each kind are in `templates/`. Unknown keys are rejected.

use serde::Deserialize;

use tpt_acoustics::aero::fwh::{rotating_ring_pressure, Ambient};
use tpt_acoustics::core::medium::Medium;
use tpt_acoustics::core::{AcousticError, Vec3};
use tpt_acoustics::ray::ism::{ShoeboxRoom, SPEED_OF_SOUND_REFERENCE};
use tpt_acoustics::ray::iso3382::{clarity_db, definition, schroeder_frequency_hz};
use tpt_acoustics::ray::material::BoundaryMaterial;
use tpt_acoustics::ray::ray::{rt60_from_curve, schroeder_curve, trace_shoebox, RayTracerConfig};
use tpt_acoustics::ray::rir::rir_from_echogram;
use tpt_acoustics::vibro::transmission::{
    mass_law_diffuse_tl, sound_transmission_class, STC_BAND_HZ,
};

/// A parsed study: one of the supported kinds.
#[derive(Debug, Deserialize, PartialEq)]
#[serde(tag = "kind", rename_all = "lowercase")]
pub enum StudyConfig {
    /// Shoebox room with a source and a receiver.
    Room(RoomStudy),
    /// Single leaf of surface mass, rated for STC.
    Partition(PartitionStudy),
    /// Rotor blade-passing tone from the FW-H ring model.
    Rotor(RotorStudy),
}

/// Settings for a room study.
#[derive(Debug, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct RoomStudy {
    /// The room.
    pub room: RoomSection,
    /// Source position.
    pub source: PointSection,
    /// Receiver position.
    pub receiver: PointSection,
    /// Analysis settings.
    pub analysis: RoomAnalysis,
}

/// Room geometry and wall absorption.
#[derive(Debug, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct RoomSection {
    /// Length, width and height in metres.
    pub dimensions: [f64; 3],
    /// Uniform absorption coefficient of all six walls, in [0, 1].
    pub absorption: f64,
}

/// A position in metres.
#[derive(Debug, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct PointSection {
    /// x, y and z in metres.
    pub position: [f64; 3],
}

/// Room analysis settings.
#[derive(Debug, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct RoomAnalysis {
    /// Frequency for absorption and the reverberation-time formulas, in Hz.
    pub frequency_hz: f64,
    /// Number of ray-tracer particles.
    pub rays: u32,
    /// Highest image-source reflection order for the echogram.
    pub max_order: u32,
    /// RNG seed, so the report is reproducible.
    pub seed: u64,
}

/// Settings for a partition study.
#[derive(Debug, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct PartitionStudy {
    /// The leaf.
    pub partition: PartitionSection,
}

/// Leaf description.
#[derive(Debug, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct PartitionSection {
    /// Surface mass in kg/m².
    pub surface_mass: f64,
}

/// Settings for a rotor study.
#[derive(Debug, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct RotorStudy {
    /// The rotor.
    pub rotor: RotorSection,
    /// Observer position in metres, off the blade path.
    pub observer: PointSection,
    /// Sampling of the pressure history.
    pub analysis: RotorAnalysis,
}

/// Rotor geometry and source strength.
#[derive(Debug, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct RotorSection {
    /// Number of blades.
    pub blades: u32,
    /// Ring radius in metres.
    pub radius_m: f64,
    /// Rotation rate in revolutions per second.
    pub rev_per_s: f64,
    /// Volume strength of each blade, in m³/s.
    pub strength_m3_per_s: f64,
}

/// Pressure-history sampling.
#[derive(Debug, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct RotorAnalysis {
    /// Sample rate in Hz.
    pub sample_rate_hz: f64,
    /// Number of samples. The window should hold an integer number of
    /// blade-passing periods for a clean projection.
    pub samples: usize,
    /// Start time of the window in seconds.
    pub start_s: f64,
}

/// Parses a study file's text.
///
/// # Errors
/// Returns the TOML parser's message for malformed input, unknown keys, or an
/// unknown or missing `kind`.
pub fn parse_study(text: &str) -> Result<StudyConfig, String> {
    toml::from_str(text).map_err(|e| e.to_string())
}

/// Runs a study and returns its report text.
///
/// # Errors
/// Returns the validation or computation error from the library.
pub fn run_study(config: &StudyConfig) -> Result<String, AcousticError> {
    match config {
        StudyConfig::Room(s) => run_room(s),
        StudyConfig::Partition(s) => run_partition(s),
        StudyConfig::Rotor(s) => run_rotor(s),
    }
}

fn run_room(config: &RoomStudy) -> Result<String, AcousticError> {
    let a = &config.analysis;
    if !(0.0..=1.0).contains(&config.room.absorption) {
        return Err(AcousticError::InvalidArgument(
            "absorption must lie in [0, 1]",
        ));
    }
    let material =
        BoundaryMaterial::new("study", flat_table(config.room.absorption), flat_table(0.0))?;
    let room = ShoeboxRoom::uniform(config.room.dimensions, material)?;
    let src = vec3(config.source.position);
    let rcv = vec3(config.receiver.position);
    let air = Medium::air_20c();

    let sabine = room.sabine_rt60(a.frequency_hz);
    let eyring = room.eyring_rt60(a.frequency_hz);
    let volume = room.volume();

    let ray_config = RayTracerConfig {
        n_rays: a.rays,
        time_limit_s: 3.0,
        capture_radius: 0.5,
        min_energy: 1e-9,
        seed: a.seed,
        speed_of_sound: SPEED_OF_SOUND_REFERENCE,
    };
    let traced = trace_shoebox(&room, src, rcv, a.frequency_hz, &ray_config)?;
    let curve = schroeder_curve(&traced.detected_times, &traced.detected_energies, 0.005);
    let ray_t30 = rt60_from_curve(&curve, 0.005);

    let echogram = room.echogram(src, rcv, a.frequency_hz, a.max_order, &air)?;
    let fs = 48_000.0;
    let rir = rir_from_echogram(&echogram, fs, a.seed)?;
    let c80 = clarity_db(&rir, fs, 80.0);
    let d50 = definition(&rir, fs, 50.0);

    let mut out = String::new();
    out.push_str(&format!(
        "room             : {:.2} x {:.2} x {:.2} m\n",
        config.room.dimensions[0], config.room.dimensions[1], config.room.dimensions[2]
    ));
    out.push_str(&format!(
        "absorption       : {:.3}\n",
        config.room.absorption
    ));
    out.push_str(&format!("frequency        : {:.0} Hz\n", a.frequency_hz));
    out.push_str(&format!("direct distance  : {:.3} m\n", src.distance(rcv)));
    out.push_str(&format!("sabine T60       : {sabine:.3} s\n"));
    out.push_str(&format!("eyring T60       : {eyring:.3} s\n"));
    match ray_t30 {
        Some(t30) => {
            let rel = 100.0 * (t30 - sabine).abs() / sabine;
            out.push_str(&format!(
                "ray-traced T30   : {t30:.3} s ({rel:.1} % from Sabine)\n"
            ));
        }
        None => out.push_str("ray-traced T30   : not enough decay to fit (raise rays or time)\n"),
    }
    match schroeder_frequency_hz(sabine, volume) {
        Some(f) => out.push_str(&format!("schroeder freq   : {f:.0} Hz\n")),
        None => out.push_str("schroeder freq   : undefined\n"),
    }
    match c80 {
        Some(c) => out.push_str(&format!(
            "C80              : {c:.2} dB (image sources, order {})\n",
            a.max_order
        )),
        None => out.push_str("C80              : undefined\n"),
    }
    match d50 {
        Some(d) => out.push_str(&format!("D50              : {d:.3}\n")),
        None => out.push_str("D50              : undefined\n"),
    }
    Ok(out)
}

fn run_partition(config: &PartitionStudy) -> Result<String, AcousticError> {
    let mass = config.partition.surface_mass;
    if !(mass.is_finite() && mass > 0.0) {
        return Err(AcousticError::InvalidArgument(
            "surface mass must be positive",
        ));
    }
    let air = Medium::air_20c();
    let spectrum: Vec<(f64, f64)> = STC_BAND_HZ
        .iter()
        .map(|&f| (f, mass_law_diffuse_tl(mass, f, &air)))
        .collect();
    let rating = sound_transmission_class(&spectrum)?;
    let mut out = format!("surface mass    : {mass:.2} kg/m2\n   f [Hz]   TL [dB]\n");
    for (f, tl) in &spectrum {
        out.push_str(&format!("{f:>9.0}   {tl:>7.1}\n"));
    }
    out.push_str(&format!("STC              : {rating}\n"));
    Ok(out)
}

fn run_rotor(config: &RotorStudy) -> Result<String, AcousticError> {
    let r = &config.rotor;
    let a = &config.analysis;
    if r.blades == 0 {
        return Err(AcousticError::InvalidArgument("blades must be at least 1"));
    }
    if !(r.radius_m.is_finite() && r.radius_m > 0.0) {
        return Err(AcousticError::InvalidArgument("radius must be positive"));
    }
    if !(r.rev_per_s.is_finite() && r.rev_per_s > 0.0) {
        return Err(AcousticError::InvalidArgument("rev_per_s must be positive"));
    }
    if !(a.sample_rate_hz.is_finite() && a.sample_rate_hz > 0.0) || a.samples < 2 {
        return Err(AcousticError::InvalidArgument(
            "sample_rate_hz must be positive and samples at least 2",
        ));
    }
    let ambient = Ambient::air_20c();
    let omega = core::f64::consts::TAU * r.rev_per_s;
    // The thickness model is only defined for a subsonic blade tip; at or
    // above Mach 1 the emission-time geometry has no valid far field.
    let tip_mach = r.radius_m * omega / ambient.c0;
    if tip_mach >= 1.0 {
        return Err(AcousticError::InvalidArgument(
            "blade tip Mach number must be below 1 (far field undefined)",
        ));
    }
    let observer = vec3(config.observer.position);
    let series: Vec<f64> = (0..a.samples)
        .map(|i| {
            let t = a.start_s + i as f64 / a.sample_rate_hz;
            rotating_ring_pressure(
                &ambient,
                r.blades,
                r.radius_m,
                omega,
                r.strength_m3_per_s,
                observer,
                t,
            )
        })
        .collect();
    if series.iter().any(|p| !p.is_finite()) {
        return Err(AcousticError::InvalidGeometry(
            "far field undefined: observer on a blade path, or blade tip supersonic",
        ));
    }
    let bpf = f64::from(r.blades) * r.rev_per_s;
    let at_bpf = tone_amplitude(&series, a.sample_rate_hz, bpf);
    let at_half = tone_amplitude(&series, a.sample_rate_hz, bpf / 2.0);
    let mut out = String::new();
    out.push_str(&format!("blade-passing    : {bpf:.1} Hz\n"));
    out.push_str(&format!("amplitude at BPF : {at_bpf:.3e} Pa\n"));
    out.push_str(&format!("amplitude BPF/2  : {at_half:.3e} Pa\n"));
    // The half-harmonic of a symmetric ring is zero in exact arithmetic; a
    // tiny value is floating-point noise, so do not print a ratio of noise.
    if at_half <= 1e-6 * at_bpf {
        out.push_str("BPF / (BPF/2)    : half-harmonic numerically zero\n");
    } else {
        out.push_str(&format!("BPF / (BPF/2)    : {:.1}\n", at_bpf / at_half));
    }
    Ok(out)
}

/// Single-frequency amplitude of `series` at `freq` (direct DFT projection).
fn tone_amplitude(series: &[f64], sample_rate: f64, freq: f64) -> f64 {
    let (mut re, mut im) = (0.0, 0.0);
    for (n, &x) in series.iter().enumerate() {
        let phase = core::f64::consts::TAU * freq * n as f64 / sample_rate;
        re += x * phase.cos();
        im -= x * phase.sin();
    }
    2.0 * (re * re + im * im).sqrt() / series.len() as f64
}

fn vec3(p: [f64; 3]) -> Vec3 {
    Vec3::new(p[0], p[1], p[2])
}

fn flat_table(alpha: f64) -> Vec<(f64, f64)> {
    vec![(125.0, alpha), (4000.0, alpha)]
}

#[cfg(test)]
mod tests {
    use super::*;

    const ROOM: &str = "\
kind = \"room\"

[room]
dimensions = [10.0, 8.0, 6.0]
absorption = 0.15

[source]
position = [2.0, 2.0, 3.0]

[receiver]
position = [8.0, 6.0, 3.0]

[analysis]
frequency_hz = 1000.0
rays = 20000
max_order = 8
seed = 7
";

    const PARTITION: &str = "\
kind = \"partition\"

[partition]
surface_mass = 5.0
";

    const ROTOR: &str = "\
kind = \"rotor\"

[rotor]
blades = 4
radius_m = 0.25
rev_per_s = 50.0
strength_m3_per_s = 1.0

[observer]
position = [2.0, 0.0, 2.0]

[analysis]
sample_rate_hz = 8192.0
samples = 8192
start_s = 0.5
";

    fn report(text: &str) -> String {
        run_study(&parse_study(text).unwrap()).unwrap()
    }

    #[test]
    fn parses_each_kind() {
        assert!(matches!(parse_study(ROOM), Ok(StudyConfig::Room(_))));
        assert!(matches!(
            parse_study(PARTITION),
            Ok(StudyConfig::Partition(_))
        ));
        assert!(matches!(parse_study(ROTOR), Ok(StudyConfig::Rotor(_))));
    }

    #[test]
    fn missing_or_unknown_kind_is_rejected() {
        assert!(parse_study(&ROOM.replace("kind = \"room\"\n", "")).is_err());
        assert!(parse_study(&ROOM.replace("\"room\"", "\"sphere\"")).is_err());
    }

    #[test]
    fn unknown_keys_are_rejected_in_every_kind() {
        assert!(parse_study(&ROOM.replace("seed = 7", "seed = 7\nsead = 3")).is_err());
        assert!(parse_study(
            &PARTITION.replace("surface_mass = 5.0", "surface_mass = 5.0\nmass = 1.0")
        )
        .is_err());
        assert!(parse_study(&ROTOR.replace("blades = 4", "blades = 4\nbalde = 1")).is_err());
    }

    #[test]
    fn room_report_contains_every_metric() {
        let out = report(ROOM);
        for key in [
            "sabine T60",
            "eyring T60",
            "ray-traced T30",
            "schroeder freq",
            "C80",
            "D50",
        ] {
            assert!(out.contains(key), "missing {key}:\n{out}");
        }
    }

    #[test]
    fn ray_t30_agrees_with_sabine_for_the_room_study() {
        let out = report(ROOM);
        let line = out
            .lines()
            .find(|l| l.starts_with("ray-traced T30"))
            .unwrap();
        let pct: f64 = line
            .split('(')
            .nth(1)
            .and_then(|s| s.split(' ').next())
            .unwrap()
            .parse()
            .unwrap();
        assert!(pct < 10.0, "{line}");
    }

    #[test]
    fn partition_report_gives_the_stc_of_a_five_kilogram_leaf() {
        let out = report(PARTITION);
        let rating: i32 = out
            .lines()
            .find_map(|l| l.strip_prefix("STC              : "))
            .unwrap()
            .trim()
            .parse()
            .unwrap();
        assert!((18..=26).contains(&rating), "STC {rating}");
    }

    #[test]
    fn rotor_report_finds_the_blade_passing_tone() {
        let out = report(ROTOR);
        assert!(out.contains("blade-passing    : 200.0 Hz"), "{out}");
        // A symmetric four-blade ring has no half-harmonic, so the report says
        // it is numerically zero rather than printing a ratio of noise.
        assert!(out.contains("half-harmonic numerically zero"), "{out}");
    }

    #[test]
    fn supersonic_tip_is_rejected() {
        // Tip Mach R·Ω / c0 = 0.25 · 2π · 250 / 343 ≈ 1.15.
        let fast = ROTOR.replace("rev_per_s = 50.0", "rev_per_s = 250.0");
        let result = run_study(&parse_study(&fast).unwrap());
        assert!(
            matches!(result, Err(AcousticError::InvalidArgument(_))),
            "{result:?}"
        );
    }

    #[test]
    fn invalid_values_are_rejected() {
        let mut cfg = match parse_study(ROOM).unwrap() {
            StudyConfig::Room(r) => r,
            _ => unreachable!(),
        };
        cfg.room.absorption = 2.0;
        assert!(run_study(&StudyConfig::Room(cfg)).is_err());
        let bad_mass = PARTITION.replace("surface_mass = 5.0", "surface_mass = -1.0");
        assert!(run_study(&parse_study(&bad_mass).unwrap()).is_err());
    }
}
