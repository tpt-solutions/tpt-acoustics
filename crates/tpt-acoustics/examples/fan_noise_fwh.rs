//! Blade-passing tone of a four-blade rotor from the rotating-ring FW-H model.
//!
//! Four blades turning at 50 rev/s give a blade-passing frequency of 200 Hz.
//! The pressure history at a point off the rotor plane is projected onto
//! 200 Hz and onto the half-harmonic at 100 Hz. The BPF must dominate.
//!
//! Run with `cargo run --release --example fan_noise_fwh`.

use tpt_acoustics::aero::fwh::{rotating_ring_pressure, Ambient};
use tpt_acoustics::core::Vec3;

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

fn main() {
    let ambient = Ambient::air_20c();
    let blades = 4;
    let omega = core::f64::consts::TAU * 50.0;
    let bpf = f64::from(blades) * 50.0;
    let radius = 0.25;
    let observer = Vec3::new(2.0, 0.0, 2.0);
    let sample_rate = 8192.0;
    let n = 8192;

    let series: Vec<f64> = (0..n)
        .map(|i| {
            let t = 0.5 + i as f64 / sample_rate;
            rotating_ring_pressure(&ambient, blades, radius, omega, 1.0, observer, t)
        })
        .collect();

    let at_bpf = tone_amplitude(&series, sample_rate, bpf);
    let at_half = tone_amplitude(&series, sample_rate, bpf / 2.0);
    println!("amplitude at BPF ({bpf:.0} Hz)     : {at_bpf:.3e} Pa");
    println!(
        "amplitude at BPF/2 ({:.0} Hz)      : {at_half:.3e} Pa",
        bpf / 2.0
    );
    assert!(at_bpf > 4.0 * at_half.max(1e-12), "BPF tone must dominate");
}
