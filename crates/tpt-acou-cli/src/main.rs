//! `tpt-acou`: command-line front end for tpt-acoustics.
//!
//! ```text
//! tpt-acou rt60 <length> <width> <height> <alpha> [freq_hz]
//! tpt-acou stc mass <kg_per_m2>
//! tpt-acou stc tl <tl_1> <tl_2> ... <tl_16>      (16 values, 125 Hz … 4 kHz)
//! tpt-acou weight <freq_hz>
//! ```
//!
//! Output goes to standard output. Errors go to standard error with exit
//! status 2 for usage problems and 1 for failed computations.

use std::process::ExitCode;

mod study;
mod wav;

use tpt_acoustics::ray::design;

use tpt_acoustics::core::medium::Medium;
use tpt_acoustics::core::AcousticError;
use tpt_acoustics::psycho::weighting::{
    a_weighting_db, b_weighting_db, c_weighting_db, itu468_weighting_db,
};
use tpt_acoustics::ray::ism::ShoeboxRoom;
use tpt_acoustics::ray::material::BoundaryMaterial;
use tpt_acoustics::vibro::transmission::{
    mass_law_diffuse_tl, sound_transmission_class, STC_BAND_HZ,
};

const USAGE: &str = "usage:
  tpt-acou rt60 <length> <width> <height> <alpha> [freq_hz]
  tpt-acou stc mass <kg_per_m2>
  tpt-acou stc tl <tl_1> ... <tl_16>
  tpt-acou weight <freq_hz>
  tpt-acou study <file.toml>                     (see templates/*.toml)
  tpt-acou design rt60 <lx> <ly> <lz> <target_s>
  tpt-acou uq rt60 <lx> <ly> <lz> <alpha_mean> <alpha_sd> <samples> <seed>
  tpt-acou auralize <dry.wav> <rir.wav> <out.wav>";

/// A usage problem (exit status 2) or a failed computation (exit status 1).
#[derive(Debug, PartialEq)]
enum CliError {
    Usage(String),
    Compute(String),
}

fn parse_number(text: &str, name: &str) -> Result<f64, CliError> {
    let value: f64 = text
        .parse()
        .map_err(|_| CliError::Usage(format!("{name} must be a number, got '{text}'")))?;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(CliError::Usage(format!("{name} must be finite")))
    }
}

fn compute(err: AcousticError) -> CliError {
    CliError::Compute(err.to_string())
}

/// Runs one command and returns its output text.
fn run(args: &[String]) -> Result<String, CliError> {
    match args.first().map(String::as_str) {
        Some("rt60") => rt60(&args[1..]),
        Some("stc") => stc(&args[1..]),
        Some("weight") => weight(&args[1..]),
        Some("study") => study_file(&args[1..]),
        Some("design") => design(&args[1..]),
        Some("uq") => uq(&args[1..]),
        Some("auralize") => auralize(&args[1..]),
        Some(other) => Err(CliError::Usage(format!("unknown command '{other}'"))),
        None => Err(CliError::Usage("no command given".into())),
    }
}

fn rt60(args: &[String]) -> Result<String, CliError> {
    if !(4..=5).contains(&args.len()) {
        return Err(CliError::Usage(
            "rt60 takes <length> <width> <height> <alpha> [freq_hz]".into(),
        ));
    }
    let dims = [
        parse_number(&args[0], "length")?,
        parse_number(&args[1], "width")?,
        parse_number(&args[2], "height")?,
    ];
    let alpha = parse_number(&args[3], "alpha")?;
    if !(0.0..=1.0).contains(&alpha) {
        return Err(CliError::Usage("alpha must lie in [0, 1]".into()));
    }
    let freq = match args.get(4) {
        Some(text) => parse_number(text, "freq_hz")?,
        None => 1000.0,
    };
    // A flat absorption coefficient across the band; scattering is zero.
    let material = BoundaryMaterial::new(
        "uniform",
        vec![(125.0, alpha), (4000.0, alpha)],
        vec![(125.0, 0.0), (4000.0, 0.0)],
    )
    .map_err(compute)?;
    let room = ShoeboxRoom::uniform(dims, material).map_err(compute)?;
    Ok(format!(
        "frequency      : {freq:.0} Hz\nsabine T60     : {:.3} s\neyring T60     : {:.3} s\n",
        room.sabine_rt60(freq),
        room.eyring_rt60(freq),
    ))
}

fn stc(args: &[String]) -> Result<String, CliError> {
    match args.first().map(String::as_str) {
        Some("mass") => {
            if args.len() != 2 {
                return Err(CliError::Usage("stc mass takes <kg_per_m2>".into()));
            }
            let mass = parse_number(&args[1], "surface mass")?;
            if mass <= 0.0 {
                return Err(CliError::Usage("surface mass must be positive".into()));
            }
            let air = Medium::air_20c();
            let spectrum: Vec<(f64, f64)> = STC_BAND_HZ
                .iter()
                .map(|&f| (f, mass_law_diffuse_tl(mass, f, &air)))
                .collect();
            let rating = sound_transmission_class(&spectrum).map_err(compute)?;
            let mut out = String::from("   f [Hz]   TL [dB]\n");
            for (f, tl) in &spectrum {
                out.push_str(&format!("{f:>9.0}   {tl:>7.1}\n"));
            }
            out.push_str(&format!("STC            : {rating}\n"));
            Ok(out)
        }
        Some("tl") => {
            if args.len() != 1 + STC_BAND_HZ.len() {
                return Err(CliError::Usage(format!(
                    "stc tl takes {} transmission-loss values",
                    STC_BAND_HZ.len()
                )));
            }
            let spectrum: Vec<(f64, f64)> = STC_BAND_HZ
                .iter()
                .zip(&args[1..])
                .map(|(&f, text)| Ok((f, parse_number(text, "transmission loss")?)))
                .collect::<Result<_, CliError>>()?;
            let rating = sound_transmission_class(&spectrum).map_err(compute)?;
            Ok(format!("STC            : {rating}\n"))
        }
        _ => Err(CliError::Usage(
            "stc takes 'mass <kg_per_m2>' or 'tl <16 values>'".into(),
        )),
    }
}

fn weight(args: &[String]) -> Result<String, CliError> {
    if args.len() != 1 {
        return Err(CliError::Usage("weight takes <freq_hz>".into()));
    }
    let f = parse_number(&args[0], "freq_hz")?;
    if f <= 0.0 {
        return Err(CliError::Usage("freq_hz must be positive".into()));
    }
    Ok(format!(
        "frequency : {f:.1} Hz\nA         : {:.2} dB\nB         : {:.2} dB\nC         : {:.2} dB\nITU-R 468 : {:.2} dB\n",
        a_weighting_db(f),
        b_weighting_db(f),
        c_weighting_db(f),
        itu468_weighting_db(f),
    ))
}

fn room_dims(args: &[String]) -> Result<[f64; 3], CliError> {
    Ok([
        parse_number(&args[0], "length")?,
        parse_number(&args[1], "width")?,
        parse_number(&args[2], "height")?,
    ])
}

/// `design rt60 <lx> <ly> <lz> <target_s>`: the uniform absorption that gives
/// the target Sabine time, and a check of the result.
fn design(args: &[String]) -> Result<String, CliError> {
    match args.first().map(String::as_str) {
        Some("rt60") if args.len() == 5 => {
            let dims = room_dims(&args[1..4])?;
            let target = parse_number(&args[4], "target T60")?;
            let alpha = design::absorption_for_rt60(dims, target).map_err(compute)?;
            let (check, _, _) = design::sabine_with_gradient(dims, [alpha; 6]).map_err(compute)?;
            Ok(format!(
                "uniform absorption : {alpha:.4}\nsabine T60 check  : {check:.4} s (target {target:.3} s)\n"
            ))
        }
        _ => Err(CliError::Usage(
            "design takes 'rt60 <lx> <ly> <lz> <target_s>'".into(),
        )),
    }
}

/// `uq rt60 <lx> <ly> <lz> <alpha_mean> <alpha_sd> <samples> <seed>`.
fn uq(args: &[String]) -> Result<String, CliError> {
    match args.first().map(String::as_str) {
        Some("rt60") if args.len() == 8 => {
            let dims = room_dims(&args[1..4])?;
            let mean = parse_number(&args[4], "alpha_mean")?;
            let sd = parse_number(&args[5], "alpha_sd")?;
            let samples: usize = args[6]
                .parse()
                .map_err(|_| CliError::Usage("samples must be a whole number".into()))?;
            let seed: u64 = args[7]
                .parse()
                .map_err(|_| CliError::Usage("seed must be a whole number".into()))?;
            let u = design::rt60_uncertainty(dims, mean, sd, samples, seed).map_err(compute)?;
            Ok(format!(
                "sabine T60 mean   : {:.4} s\nsabine T60 std    : {:.4} s\n5th percentile    : {:.4} s\n95th percentile   : {:.4} s\n",
                u.mean_s, u.std_s, u.p05_s, u.p95_s
            ))
        }
        _ => Err(CliError::Usage(
            "uq takes 'rt60 <lx> <ly> <lz> <alpha_mean> <alpha_sd> <samples> <seed>'".into(),
        )),
    }
}

/// `auralize <dry.wav> <rir.wav> <out.wav>`: convolves the dry signal with the
/// impulse response and writes a float WAV, scaled to a peak of 0.99 if the
/// result would clip.
fn auralize(args: &[String]) -> Result<String, CliError> {
    if args.len() != 3 {
        return Err(CliError::Usage(
            "auralize takes <dry.wav> <rir.wav> <out.wav>".into(),
        ));
    }
    let read = |path: &String| -> Result<wav::WavData, CliError> {
        let bytes = std::fs::read(path)
            .map_err(|e| CliError::Usage(format!("cannot read '{path}': {e}")))?;
        wav::parse_wav(&bytes).map_err(|e| CliError::Usage(format!("'{path}': {e}")))
    };
    let dry = read(&args[0])?;
    let rir = read(&args[1])?;
    if dry.sample_rate != rir.sample_rate {
        return Err(CliError::Usage(format!(
            "sample rates differ: {} Hz dry, {} Hz RIR",
            dry.sample_rate, rir.sample_rate
        )));
    }
    let mut out = auralize_signal(&dry.samples, &rir.samples);
    let peak = out.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
    let gain = if peak > 0.99 { 0.99 / peak } else { 1.0 };
    for v in &mut out {
        *v *= gain;
    }
    let bytes = tpt_acoustics::ray::rir::write_wav_f32(&out, dry.sample_rate);
    std::fs::write(&args[2], bytes)
        .map_err(|e| CliError::Usage(format!("cannot write '{}': {e}", args[2])))?;
    Ok(format!(
        "wrote {} samples at {} Hz to {} (gain {gain:.3})\n",
        out.len(),
        dry.sample_rate,
        args[2]
    ))
}

fn auralize_signal(dry: &[f64], rir: &[f64]) -> Vec<f64> {
    tpt_acoustics::ray::rir::auralize(dry, rir)
}

fn study_file(args: &[String]) -> Result<String, CliError> {
    if args.len() != 1 {
        return Err(CliError::Usage("study takes <file.toml>".into()));
    }
    let text = std::fs::read_to_string(&args[0])
        .map_err(|e| CliError::Usage(format!("cannot read '{}': {e}", args[0])))?;
    let config =
        study::parse_study(&text).map_err(|e| CliError::Usage(format!("study file: {e}")))?;
    study::run_study(&config).map_err(compute)
}

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().skip(1).collect();
    match run(&args) {
        Ok(text) => {
            print!("{text}");
            ExitCode::SUCCESS
        }
        Err(CliError::Usage(message)) => {
            eprintln!("error: {message}\n{USAGE}");
            ExitCode::from(2)
        }
        Err(CliError::Compute(message)) => {
            eprintln!("error: {message}");
            ExitCode::from(1)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn argv(text: &str) -> Vec<String> {
        text.split_whitespace().map(str::to_owned).collect()
    }

    #[test]
    fn rt60_reports_sabine_and_eyring() {
        let out = run(&argv("rt60 10 8 6 0.15")).unwrap();
        assert!(out.contains("sabine T60"));
        assert!(out.contains("eyring T60"));
        // Sabine for 10×8×6 m at α = 0.15 is about 1.37 s (see the ray tests).
        assert!(
            out.contains("1.37") || out.contains("1.36") || out.contains("1.38"),
            "{out}"
        );
    }

    #[test]
    fn rt60_rejects_bad_alpha() {
        assert!(matches!(
            run(&argv("rt60 10 8 6 1.5")),
            Err(CliError::Usage(_))
        ));
    }

    #[test]
    fn stc_from_mass_is_in_the_expected_range() {
        let out = run(&argv("stc mass 5")).unwrap();
        let rating: i32 = out
            .lines()
            .find_map(|l| l.strip_prefix("STC            : "))
            .unwrap()
            .trim()
            .parse()
            .unwrap();
        assert!((18..=26).contains(&rating), "STC {rating}");
    }

    #[test]
    fn stc_from_tl_matches_the_analytic_contour() {
        // Spectrum equal to the STC-50 contour rates 52 (see the vibro tests).
        let values = [
            34.0, 37.0, 40.0, 43.0, 46.0, 49.0, 50.0, 51.0, 52.0, 53.0, 54.0, 54.0, 54.0, 54.0,
            54.0, 54.0,
        ];
        let mut args = vec!["stc".to_string(), "tl".to_string()];
        args.extend(values.iter().map(|v| v.to_string()));
        let out = run(&args).unwrap();
        assert!(out.ends_with("STC            : 52\n"), "{out}");
    }

    #[test]
    fn stc_tl_needs_sixteen_values() {
        assert!(matches!(
            run(&argv("stc tl 40 40")),
            Err(CliError::Usage(_))
        ));
    }

    #[test]
    fn weight_at_one_kilohertz_is_zero_for_a() {
        let out = run(&argv("weight 1000")).unwrap();
        assert!(out.contains("A         : 0.00 dB"), "{out}");
    }

    #[test]
    fn unknown_command_is_a_usage_error() {
        assert!(matches!(run(&argv("frobnicate")), Err(CliError::Usage(_))));
        assert!(matches!(run(&[]), Err(CliError::Usage(_))));
    }
}
