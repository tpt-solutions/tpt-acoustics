//! Python bindings for tpt-acoustics.
//!
//! Errors raise `ValueError` with the library's message, so Python callers
//! see the same reasons the Rust API reports.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use tpt_acoustics::core::medium::Medium;
use tpt_acoustics::core::AcousticError;
use tpt_acoustics::psycho::loudness_moore::moore_loudness_tones;
use tpt_acoustics::psycho::weighting::a_weighting_db;
use tpt_acoustics::ray::design::sabine_with_gradient;
use tpt_acoustics::ray::iso3382::schroeder_frequency_hz;
use tpt_acoustics::ray::sti::{sti_from_rir, STI_WEIGHTS_MALE};
use tpt_acoustics::vibro::transmission::{
    mass_law_diffuse_tl, sound_transmission_class, STC_BAND_HZ,
};

fn to_py(err: AcousticError) -> PyErr {
    PyValueError::new_err(err.to_string())
}

/// Sabine reverberation time (s) of a box room with uniform absorption.
#[pyfunction]
fn sabine_rt60(dims: [f64; 3], alpha: f64) -> PyResult<f64> {
    sabine_with_gradient(dims, [alpha; 6])
        .map(|(t, _, _)| t)
        .map_err(to_py)
}

/// Schroeder frequency (Hz) for reverberation time and volume.
#[pyfunction]
fn schroeder_frequency(t60_s: f64, volume_m3: f64) -> PyResult<f64> {
    schroeder_frequency_hz(t60_s, volume_m3)
        .ok_or_else(|| PyValueError::new_err("T60 and volume must be positive"))
}

/// Diffuse-field STC of a single leaf of surface mass (kg/m²).
#[pyfunction]
fn stc_from_mass(kg_per_m2: f64) -> PyResult<i32> {
    if !(kg_per_m2.is_finite() && kg_per_m2 > 0.0) {
        return Err(PyValueError::new_err("surface mass must be positive"));
    }
    let air = Medium::air_20c();
    let spectrum: Vec<(f64, f64)> = STC_BAND_HZ
        .iter()
        .map(|&f| (f, mass_law_diffuse_tl(kg_per_m2, f, &air)))
        .collect();
    sound_transmission_class(&spectrum).map_err(to_py)
}

/// A-weighting correction (dB) at a frequency in hertz.
#[pyfunction]
fn a_weighting(freq_hz: f64) -> PyResult<f64> {
    if !(freq_hz.is_finite() && freq_hz > 0.0) {
        return Err(PyValueError::new_err("frequency must be positive"));
    }
    Ok(a_weighting_db(freq_hz))
}

/// Moore–Glasberg loudness (sones) of one pure tone.
#[pyfunction]
fn loudness_of_tone(freq_hz: f64, level_db_spl: f64) -> PyResult<f64> {
    moore_loudness_tones(&[(freq_hz, level_db_spl)]).map_err(to_py)
}

/// Speech Transmission Index (male weights) of an impulse response.
#[pyfunction]
fn sti(impulse_response: Vec<f64>, sample_rate_hz: f64) -> PyResult<f64> {
    sti_from_rir(&impulse_response, sample_rate_hz, &STI_WEIGHTS_MALE).map_err(to_py)
}

#[pymodule]
#[pyo3(name = "tpt_acoustics")]
fn tpt_acoustics_module(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(sabine_rt60, m)?)?;
    m.add_function(wrap_pyfunction!(schroeder_frequency, m)?)?;
    m.add_function(wrap_pyfunction!(stc_from_mass, m)?)?;
    m.add_function(wrap_pyfunction!(a_weighting, m)?)?;
    m.add_function(wrap_pyfunction!(loudness_of_tone, m)?)?;
    m.add_function(wrap_pyfunction!(sti, m)?)?;
    Ok(())
}
