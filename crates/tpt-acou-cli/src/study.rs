//! Room study from a TOML description.
//!
//! A study file names one shoebox room, a source, a receiver and the analysis
//! settings. [`run_study`] returns a text report with the Sabine and Eyring
//! reverberation times, a ray-traced T30 checked against Sabine, the Schroeder
//! frequency, and clarity and definition from the image-source response. The
//! file format is documented in `templates/room-study.toml`.

use serde::Deserialize;

use tpt_acoustics::core::medium::Medium;
use tpt_acoustics::core::AcousticError;
use tpt_acoustics::ray::ism::{ShoeboxRoom, SPEED_OF_SOUND_REFERENCE};
use tpt_acoustics::ray::iso3382::{clarity_db, definition, schroeder_frequency_hz};
use tpt_acoustics::ray::material::BoundaryMaterial;
use tpt_acoustics::ray::ray::{rt60_from_curve, schroeder_curve, trace_shoebox, RayTracerConfig};
use tpt_acoustics::ray::rir::rir_from_echogram;
use tpt_acoustics::ray::vec3::Vec3;

/// Parsed study file. Unknown keys are rejected so typos fail loudly.
#[derive(Debug, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct StudyConfig {
    /// The room.
    pub room: RoomSection,
    /// Source position.
    pub source: PointSection,
    /// Receiver position.
    pub receiver: PointSection,
    /// Analysis settings.
    pub analysis: AnalysisSection,
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
    /// x, y and z in metres, strictly inside the room.
    pub position: [f64; 3],
}

/// Analysis settings.
#[derive(Debug, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct AnalysisSection {
    /// Frequency for absorption and the reverberation-time formulas, in Hz.
    pub frequency_hz: f64,
    /// Number of ray-tracer particles.
    pub rays: u32,
    /// Highest image-source reflection order for the echogram.
    pub max_order: u32,
    /// RNG seed, so the report is reproducible.
    pub seed: u64,
}

/// Parses a study file's text.
///
/// # Errors
/// Returns the TOML parser's message for malformed input or unknown keys.
pub fn parse_study(text: &str) -> Result<StudyConfig, String> {
    toml::from_str(text).map_err(|e| e.to_string())
}

/// Runs the study and returns the report text.
///
/// # Errors
/// Returns the geometry or validation error from the library as a string.
pub fn run_study(config: &StudyConfig) -> Result<String, AcousticError> {
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
                "ray-traced T30    : {t30:.3} s ({rel:.1} % from Sabine)\n"
            ));
        }
        None => out.push_str("ray-traced T30    : not enough decay to fit (raise rays or time)\n"),
    }
    match schroeder_frequency_hz(sabine, volume) {
        Some(fs_hz) => out.push_str(&format!("schroeder freq   : {fs_hz:.0} Hz\n")),
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

fn vec3(p: [f64; 3]) -> Vec3 {
    Vec3::new(p[0], p[1], p[2])
}

fn flat_table(alpha: f64) -> Vec<(f64, f64)> {
    vec![(125.0, alpha), (4000.0, alpha)]
}

#[cfg(test)]
mod tests {
    use super::*;

    const STUDY: &str = "\
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

    #[test]
    fn parses_the_template_layout() {
        let cfg = parse_study(STUDY).unwrap();
        assert_eq!(cfg.room.dimensions, [10.0, 8.0, 6.0]);
        assert_eq!(cfg.analysis.rays, 20000);
    }

    #[test]
    fn unknown_keys_are_rejected() {
        let bad = STUDY.replace("seed = 7", "seed = 7\nsead = 3");
        assert!(parse_study(&bad).is_err());
    }

    #[test]
    fn report_contains_every_metric() {
        let cfg = parse_study(STUDY).unwrap();
        let report = run_study(&cfg).unwrap();
        for key in [
            "sabine T60",
            "eyring T60",
            "ray-traced T30",
            "schroeder freq",
            "C80",
            "D50",
        ] {
            assert!(report.contains(key), "missing {key}:\n{report}");
        }
    }

    #[test]
    fn ray_t30_agrees_with_sabine_for_the_study_room() {
        let cfg = parse_study(STUDY).unwrap();
        let report = run_study(&cfg).unwrap();
        let line = report
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
    fn bad_absorption_is_rejected() {
        let mut cfg = parse_study(STUDY).unwrap();
        cfg.room.absorption = 2.0;
        assert!(run_study(&cfg).is_err());
    }
}
