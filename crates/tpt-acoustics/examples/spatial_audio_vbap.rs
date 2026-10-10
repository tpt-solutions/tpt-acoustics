//! Vector Base Amplitude Panning (VBAP) on a square horizontal loudspeaker ring.
//!
//! Four speakers at 0°, 90°, 180° and 270° azimuth. A source on a speaker
//! plays from that speaker alone; a source at 45° is shared equally by the
//! two neighbouring speakers. Gains are normalised so that their Euclidean
//! norm is one (Pulkki's constant-power convention).
//!
//! Run with `cargo run --example spatial_audio_vbap`.

use tpt_acoustics::core::{AcousticError, Vec3};
use tpt_acoustics::psycho::vbap::VbapLayout;

/// Unit vector on the horizontal plane at `azimuth_deg` (0° = +x, counter-clockwise).
fn horizontal(azimuth_deg: f64) -> Vec3 {
    let a = azimuth_deg.to_radians();
    Vec3::new(a.cos(), a.sin(), 0.0)
}

fn main() -> Result<(), AcousticError> {
    let speakers: Vec<Vec3> = [0.0, 90.0, 180.0, 270.0]
        .into_iter()
        .map(horizontal)
        .collect();
    let layout = VbapLayout::new_2d(speakers)?;

    for azimuth in [0.0, 45.0, 90.0, 135.0] {
        let gains = layout.gains(horizontal(azimuth));
        let norm = gains.iter().map(|g| g * g).sum::<f64>().sqrt();
        let formatted: Vec<String> = gains.iter().map(|g| format!("{g:.3}")).collect();
        println!(
            "source {azimuth:>5.1}°  gains [{}]  |g| = {norm:.6}",
            formatted.join(", ")
        );
        assert!((norm - 1.0).abs() < 1e-9, "gain norm must be one");
    }

    // The layout sorts speakers by azimuth, so look up indices rather than
    // assuming the input order.
    let index_of = |azimuth: f64| {
        let target = horizontal(azimuth);
        layout
            .speakers()
            .iter()
            .position(|s| (*s - target).norm() < 1e-9)
            .expect("speaker present in layout")
    };

    // On a speaker: that speaker alone at unit gain.
    let on_speaker = layout.gains(horizontal(0.0));
    assert!((on_speaker[index_of(0.0)] - 1.0).abs() < 1e-9);
    let others: f64 = (0..4)
        .filter(|&i| i != index_of(0.0))
        .map(|i| on_speaker[i].abs())
        .sum();
    assert!(others < 1e-9, "other speakers must be silent");

    // Between two: an equal split of the two neighbours.
    let between = layout.gains(horizontal(45.0));
    let expected = core::f64::consts::FRAC_1_SQRT_2;
    assert!((between[index_of(0.0)] - expected).abs() < 1e-9);
    assert!((between[index_of(90.0)] - expected).abs() < 1e-9);
    assert!(between[index_of(180.0)].abs() < 1e-9);
    assert!(between[index_of(270.0)].abs() < 1e-9);
    Ok(())
}
