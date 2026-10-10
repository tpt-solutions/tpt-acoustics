//! Differentiable room design: recover a room length from a target energy.
//!
//! `ism_energy_and_gradient` returns the image-source energy at a receiver
//! and its exact derivative with respect to the room length `L` along x,
//! computed by forward-mode automatic differentiation. This example first
//! checks that derivative against a finite difference. It then runs damped
//! Newton iterations on `E(L) − E*` to recover the length of a known room.
//!
//! The energy oscillates with `L` on the scale of the acoustic wavelength,
//! because image paths interfere. Newton only converges near the starting
//! point, and a target taken from a far-away length is not reachable from 4 m.
//! Global design needs a scan over `L` before refinement.
//!
//! Run with `cargo run --release --example optimize_room_gradient`.

use tpt_acoustics::core::AcousticError;
use tpt_acoustics::ray::differentiable::ism_energy_and_gradient;
use tpt_acoustics::ray::ism::ShoeboxRoom;
use tpt_acoustics::ray::material::BoundaryMaterial;
use tpt_acoustics::ray::vec3::Vec3;

const FREQ_HZ: f64 = 1000.0;
const SPEED: f64 = 343.0;
const MAX_ORDER: u32 = 4;
const WIDTH: f64 = 5.0;
const HEIGHT: f64 = 3.0;

fn room(length: f64) -> Result<ShoeboxRoom, AcousticError> {
    ShoeboxRoom::uniform(
        [length, WIDTH, HEIGHT],
        BoundaryMaterial::concrete_painted(),
    )
}

/// Energy and its exact derivative with respect to the room length.
fn energy(length: f64, src: Vec3, rcv: Vec3) -> Result<(f64, f64), AcousticError> {
    ism_energy_and_gradient(&room(length)?, src, rcv, FREQ_HZ, SPEED, MAX_ORDER, 0)
}

fn main() -> Result<(), AcousticError> {
    let src = Vec3::new(1.0, 2.0, 1.5);
    let rcv = Vec3::new(3.0, 2.5, 1.5);
    // Both points must stay inside the room, so the length cannot drop
    // below the x-coordinate of either of them.
    let min_length = src.x.max(rcv.x) + 0.1;

    // 1. The analytic gradient must match a central finite difference.
    let l0 = 4.0;
    let h = 1e-5;
    let (e0, grad) = energy(l0, src, rcv)?;
    let fd = (energy(l0 + h, src, rcv)?.0 - energy(l0 - h, src, rcv)?.0) / (2.0 * h);
    println!("E({l0}) = {e0:.6e}, dE/dL analytic = {grad:.6e}, finite diff = {fd:.6e}");
    assert!(
        (grad - fd).abs() <= 1e-4 * fd.abs().max(1e-12),
        "gradient mismatch: {grad} vs {fd}"
    );

    // 2. Local recovery. The energy oscillates with length on the scale of the
    //    acoustic wavelength (image paths interfere), so refinement only works
    //    near the starting point. Pick a target that a room 2 cm longer produces
    //    and recover that length from 4 m.
    let known_length = 4.02;
    let target = energy(known_length, src, rcv)?.0;
    let mut length = l0;
    // Damped Newton on the residual r(L) = E(L) − E*, using the exact
    // derivative. Each step is halved until |r| decreases and the length
    // stays above the source and receiver.
    for iter in 0..50 {
        let (e, de_dl) = energy(length, src, rcv)?;
        let residual = e - target;
        if residual.abs() <= 1e-12 * target || de_dl.abs() < 1e-30 {
            break;
        }
        let full_step = residual / de_dl;
        let mut scale = 1.0;
        let mut accepted = false;
        for _ in 0..30 {
            let candidate = length - scale * full_step;
            if candidate > min_length
                && (energy(candidate, src, rcv)?.0 - target).abs() < residual.abs()
            {
                length = candidate;
                accepted = true;
                break;
            }
            scale *= 0.5;
        }
        if !accepted {
            break;
        }
        if iter % 5 == 0 {
            println!(
                "iter {iter:>2}: L = {length:.6} m, |E - E*| = {:.3e}",
                residual.abs()
            );
        }
    }

    let final_energy = energy(length, src, rcv)?.0;
    let rel = (final_energy - target).abs() / target;
    println!("recovered length: {length:.6} m (target room {known_length} m)");
    println!("energy error   : {:.2e} relative", rel);
    assert!(
        (length - known_length).abs() < 1e-4,
        "recovered {length} m, expected {known_length} m"
    );
    assert!(rel < 1e-6, "energy mismatch {rel:.2e}");
    Ok(())
}
