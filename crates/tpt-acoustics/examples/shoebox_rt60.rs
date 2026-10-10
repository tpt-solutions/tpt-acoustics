//! Reverberation time of a shoebox room: ray-traced T30 against Sabine.
//!
//! A 10 × 8 × 6 m room with uniform absorption α = 0.15. Sabine gives
//! T60 ≈ 1.37 s, and the Monte Carlo tracer's T30 should agree within 10 %.
//!
//! Run with `cargo run --release --example shoebox_rt60`.

use tpt_acoustics::core::AcousticError;
use tpt_acoustics::ray::ism::{ShoeboxRoom, SPEED_OF_SOUND_REFERENCE};
use tpt_acoustics::ray::material::BoundaryMaterial;
use tpt_acoustics::ray::ray::{rt60_from_curve, schroeder_curve, trace_shoebox, RayTracerConfig};
use tpt_acoustics::ray::vec3::Vec3;

fn main() -> Result<(), AcousticError> {
    let alpha = 0.15;
    let material = BoundaryMaterial::new(
        "uniform",
        vec![(125.0, alpha), (4000.0, alpha)],
        vec![(125.0, 0.0), (4000.0, 0.0)],
    )?;
    let room = ShoeboxRoom::uniform([10.0, 8.0, 6.0], material)?;

    let sabine = room.sabine_rt60(1000.0);
    let config = RayTracerConfig {
        n_rays: 30_000,
        time_limit_s: 3.0,
        capture_radius: 0.5,
        min_energy: 1e-9,
        seed: 7,
        speed_of_sound: SPEED_OF_SOUND_REFERENCE,
    };
    let result = trace_shoebox(
        &room,
        Vec3::new(2.0, 2.0, 3.0),
        Vec3::new(8.0, 6.0, 3.0),
        1000.0,
        &config,
    )?;
    let curve = schroeder_curve(&result.detected_times, &result.detected_energies, 0.005);
    let t30 = rt60_from_curve(&curve, 0.005).expect("decay long enough for T30");

    let rel = (t30 - sabine).abs() / sabine;
    println!("Sabine T60 at 1 kHz : {sabine:.3} s");
    println!("Ray-traced T30      : {t30:.3} s");
    println!("Relative difference : {:.1} %", 100.0 * rel);
    assert!(
        rel < 0.10,
        "T30 differs from Sabine by {:.1} %",
        100.0 * rel
    );
    Ok(())
}
