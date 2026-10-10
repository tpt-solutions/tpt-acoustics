//! WebAssembly bindings for tpt-acoustics.
//!
//! The exports take and return plain numbers, so no JavaScript objects cross
//! the boundary. Where a computation can fail, the function returns `NaN`, and
//! the JavaScript caller checks with `Number.isNaN`. The Rust API returns
//! `Result` for the same cases.
//!
//! Build for the web with `wasm-pack build crates/tpt-acou-wasm --target web`,
//! or check the target with
//! `cargo build -p tpt-acou-wasm --target wasm32-unknown-unknown`.

use wasm_bindgen::prelude::*;

use tpt_acoustics::core::medium::Medium;
use tpt_acoustics::psycho::loudness_moore::moore_loudness_tones;
use tpt_acoustics::psycho::weighting::a_weighting_db;
use tpt_acoustics::ray::design::sabine_with_gradient;
use tpt_acoustics::ray::iso3382::schroeder_frequency_hz;
use tpt_acoustics::ray::sti::{sti_from_rir, STI_WEIGHTS_MALE};
use tpt_acoustics::vibro::transmission::{
    mass_law_diffuse_tl, sound_transmission_class, STC_BAND_HZ,
};

/// Sabine reverberation time (s) of a box room with uniform absorption
/// `alpha`, or `NaN` for invalid input.
#[wasm_bindgen]
pub fn sabine_rt60(lx: f64, ly: f64, lz: f64, alpha: f64) -> f64 {
    let alphas = [alpha; 6];
    match sabine_with_gradient([lx, ly, lz], alphas) {
        Ok((t, _, _)) => t,
        Err(_) => f64::NAN,
    }
}

/// Schroeder frequency (Hz) for reverberation time `t60_s` and volume
/// `volume_m3`, or `NaN` for invalid input.
#[wasm_bindgen]
pub fn schroeder_frequency(t60_s: f64, volume_m3: f64) -> f64 {
    schroeder_frequency_hz(t60_s, volume_m3).unwrap_or(f64::NAN)
}

/// Diffuse-field STC of a single leaf of surface mass `kg_per_m2`, or `-1`
/// for invalid input (a rating is never negative).
#[wasm_bindgen]
pub fn stc_from_mass(kg_per_m2: f64) -> i32 {
    if !(kg_per_m2.is_finite() && kg_per_m2 > 0.0) {
        return -1;
    }
    let air = Medium::air_20c();
    let spectrum: Vec<(f64, f64)> = STC_BAND_HZ
        .iter()
        .map(|&f| (f, mass_law_diffuse_tl(kg_per_m2, f, &air)))
        .collect();
    sound_transmission_class(&spectrum).unwrap_or(-1)
}

/// A-weighting correction (dB) at `freq_hz`, or `NaN` for a non-positive
/// frequency.
#[wasm_bindgen]
pub fn a_weighting(freq_hz: f64) -> f64 {
    if freq_hz.is_finite() && freq_hz > 0.0 {
        a_weighting_db(freq_hz)
    } else {
        f64::NAN
    }
}

/// Moore–Glasberg loudness (sones) of one pure tone, or `NaN` for invalid
/// input.
#[wasm_bindgen]
pub fn loudness_of_tone(freq_hz: f64, level_db_spl: f64) -> f64 {
    moore_loudness_tones(&[(freq_hz, level_db_spl)]).unwrap_or(f64::NAN)
}

/// Speech Transmission Index (male weights) of an impulse response sampled at
/// `sample_rate_hz`, or `NaN` if the input is invalid (for example a sample
/// rate at or below 22.6 kHz, or an empty response).
#[wasm_bindgen]
pub fn sti(impulse_response: &[f64], sample_rate_hz: f64) -> f64 {
    sti_from_rir(impulse_response, sample_rate_hz, &STI_WEIGHTS_MALE).unwrap_or(f64::NAN)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sabine_matches_the_library() {
        let t = sabine_rt60(10.0, 8.0, 6.0, 0.15);
        assert!((t - 1.370).abs() < 1e-2, "{t}");
        assert!(sabine_rt60(10.0, 8.0, 6.0, 0.0).is_nan());
    }

    #[test]
    fn invalid_inputs_are_reported_not_panicked() {
        assert!(schroeder_frequency(-1.0, 60.0).is_nan());
        assert_eq!(stc_from_mass(-2.0), -1);
        assert!(a_weighting(0.0).is_nan());
        assert!(loudness_of_tone(0.0, 40.0).is_nan());
        assert!(sti(&[], 48_000.0).is_nan());
    }

    #[test]
    fn anchors_agree_with_the_library() {
        assert!(a_weighting(1000.0).abs() < 0.01);
        assert!((loudness_of_tone(1000.0, 40.0) - 1.0).abs() < 1e-9);
        assert!((schroeder_frequency(0.5, 200.0) - 100.0).abs() < 1e-9);
    }
}
