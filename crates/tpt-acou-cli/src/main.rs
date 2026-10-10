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
  tpt-acou weight <freq_hz>";

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
