//! Image-source echogram to room impulse response, with ISO 3382 metrics.
//!
//! Builds an energy echogram for an 8 × 5 × 3 m room, turns it into a
//! sampled impulse response at 48 kHz, and reports clarity and definition.
//! The response is also serialised as a float WAV to show the export path.
//! Nothing is written to disk.
//!
//! Reverberation times are deliberately not reported here. The image-source
//! echogram stops at its highest-order arrival (about 0.5 s for this room),
//! so the response has no late tail and a decay fit would be meaningless.
//! Use the ray tracer's late energy for T30 (see the `shoebox_rt60` example).
//!
//! Run with `cargo run --release --example room_impulse_response`.

use tpt_acoustics::core::medium::Medium;
use tpt_acoustics::core::AcousticError;
use tpt_acoustics::ray::ism::ShoeboxRoom;
use tpt_acoustics::ray::iso3382::{clarity_db, definition};
use tpt_acoustics::ray::material::BoundaryMaterial;
use tpt_acoustics::ray::rir::{rir_from_echogram, write_wav_f32};
use tpt_acoustics::ray::vec3::Vec3;

fn main() -> Result<(), AcousticError> {
    let room = ShoeboxRoom::uniform([8.0, 5.0, 3.0], BoundaryMaterial::concrete_painted())?;
    let src = Vec3::new(1.5, 1.2, 1.4);
    let rcv = Vec3::new(6.2, 3.7, 1.6);
    let air = Medium::air_20c();

    let echogram = room.echogram(src, rcv, 1000.0, 20, &air)?;
    let fs = 48_000.0;
    let rir = rir_from_echogram(&echogram, fs, 1)?;

    let c50 = clarity_db(&rir, fs, 50.0).expect("C50");
    let c80 = clarity_db(&rir, fs, 80.0).expect("C80");
    let d50 = definition(&rir, fs, 50.0).expect("D50");

    println!("impulses            : {}", echogram.len());
    println!("RIR length          : {:.3} s", rir.len() as f64 / fs);
    println!("C50                 : {c50:.2} dB");
    println!("C80                 : {c80:.2} dB");
    println!("D50                 : {d50:.3}");

    assert!(c80 > c50, "C80 must exceed C50 for the same response");
    assert!((0.0..=1.0).contains(&d50), "D50 is a ratio");

    let wav = write_wav_f32(&rir, fs as u32);
    println!("WAV bytes           : {}", wav.len());
    assert_eq!(wav.len(), 44 + 4 * rir.len());
    Ok(())
}
