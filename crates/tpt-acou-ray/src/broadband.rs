//! Broadband hybrid impulse response: modal below the Schroeder frequency,
//! image sources above it.
//!
//! The two transfer functions share the free-field Green's-function
//! normalisation, so their spectra join directly. The combined spectrum has
//! Hermitian symmetry, and the inverse FFT gives a real impulse response.
//!
//! * Below `f_S = 2000 √(T/V)` the response is modal, and the rigid-box modal
//!   sum ([`crate::modal::RigidBoxModes`]) is used, with the Sabine damping of
//!   the room at 1 kHz.
//! * Above `f_S` the image-source transfer function is used, with reflection
//!   coefficients from the wall absorption.
//!
//! **Limits.** The high band uses image sources only. The ray-traced late
//! reverberation is not merged here, because the ray tracer is energy-only and
//! has no phase. The crossover is a hard cut at `f_S`. The result is a
//! frequency-domain construction of the Green's function, with no separate
//! level calibration between the two bands.

use alloc::vec;
use alloc::vec::Vec;

use tpt_acou_core::error::AcousticError;
use tpt_acou_core::medium::Medium;
use tpt_math_linalg_complex::Complex;

use crate::fft::{fft_in_place, next_pow2};
use crate::ism::ShoeboxRoom;
use crate::modal::RigidBoxModes;
use crate::vec3::Vec3;

/// Broadband impulse response of `room` between `src` and `rcv`.
///
/// `fs` is the sample rate, and the output has `n_samples` samples, rounded up
/// to a power of two internally and truncated back. The image-source order is
/// `max_order`.
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] for non-positive `fs` or
/// `n_samples`, and propagates geometry errors for positions outside the room.
pub fn hybrid_broadband_ir(
    room: &ShoeboxRoom,
    src: Vec3,
    rcv: Vec3,
    fs: f64,
    n_samples: usize,
    max_order: u32,
    medium: &Medium,
) -> Result<Vec<f64>, AcousticError> {
    if !(fs.is_finite() && fs > 0.0) || n_samples == 0 {
        return Err(AcousticError::InvalidArgument(
            "fs must be positive and n_samples at least 1",
        ));
    }
    let dims = room.dims();
    let volume = dims[0] * dims[1] * dims[2];
    let t60 = room.sabine_rt60(1000.0);
    let f_schroeder = crate::iso3382::schroeder_frequency_hz(t60, volume)
        .ok_or(AcousticError::InvalidArgument(
            "room has no finite reverberation time",
        ))?
        .min(fs / 2.0);

    let n = next_pow2(n_samples);
    let half = n / 2;
    let modes = RigidBoxModes::new(
        dims,
        [src.x, src.y, src.z],
        [rcv.x, rcv.y, rcv.z],
        medium.speed_of_sound,
        2.0 * f_schroeder,
        t60,
    )?;

    // Spectrum on the FFT grid, bins 0..=N/2.
    let mut spectrum = vec![Complex::new(0.0, 0.0); n];
    for k in 1..=half {
        let f = k as f64 * fs / n as f64;
        let value = if f < f_schroeder {
            modes.green_at_hz(f)
        } else {
            room.transfer_function(src, rcv, &[f], max_order, medium)?[0]
        };
        spectrum[k] = value;
        if k != half {
            spectrum[n - k] = value.conj();
        }
    }
    // Inverse transform: h[n] = fs · IFFT(H)[n].
    fft_in_place(&mut spectrum, true);
    let out: Vec<f64> = spectrum.iter().take(n_samples).map(|c| c.re * fs).collect();
    Ok(out)
}

/// The spectrum on the same grid as [`hybrid_broadband_ir`] (bins 0..=N/2 for
/// the padded length), for checks and plots. It uses the same crossover.
///
/// # Errors
/// As [`hybrid_broadband_ir`].
pub fn hybrid_broadband_spectrum(
    room: &ShoeboxRoom,
    src: Vec3,
    rcv: Vec3,
    fs: f64,
    n_samples: usize,
    max_order: u32,
    medium: &Medium,
) -> Result<Vec<Complex>, AcousticError> {
    let ir = hybrid_broadband_ir(room, src, rcv, fs, n_samples, max_order, medium)?;
    let n = next_pow2(n_samples);
    let mut padded: Vec<Complex> = vec![Complex::new(0.0, 0.0); n];
    for (slot, &h) in padded.iter_mut().zip(&ir) {
        *slot = Complex::new(h, 0.0);
    }
    fft_in_place(&mut padded, false);
    Ok(padded.iter().map(|c| *c * (1.0 / fs)).collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::BoundaryMaterial;

    fn room() -> ShoeboxRoom {
        ShoeboxRoom::uniform([5.0, 4.0, 3.0], BoundaryMaterial::concrete_painted()).unwrap()
    }

    #[test]
    fn output_is_real_and_has_the_requested_length() {
        let air = Medium::air_20c();
        let ir = hybrid_broadband_ir(
            &room(),
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(4.0, 3.0, 2.0),
            2_000.0,
            1_000,
            2,
            &air,
        )
        .unwrap();
        assert_eq!(ir.len(), 1_000);
        assert!(ir.iter().all(|v| v.is_finite()));
    }

    #[test]
    fn spectrum_matches_the_two_transfer_functions_in_their_bands() {
        let air = Medium::air_20c();
        let r = room();
        let src = Vec3::new(1.0, 1.0, 1.0);
        let rcv = Vec3::new(4.0, 3.0, 2.0);
        let fs = 2_000.0;
        let n_samples = 1_024;
        let spec = hybrid_broadband_spectrum(&r, src, rcv, fs, n_samples, 3, &air).unwrap();
        let n = next_pow2(n_samples);
        let t60 = r.sabine_rt60(1000.0);
        let f_s = crate::iso3382::schroeder_frequency_hz(t60, 60.0)
            .unwrap()
            .min(fs / 2.0);
        let modes = RigidBoxModes::new(
            r.dims(),
            [src.x, src.y, src.z],
            [rcv.x, rcv.y, rcv.z],
            air.speed_of_sound,
            2.0 * f_s,
            t60,
        )
        .unwrap();
        // Bins well inside each band (away from the crossover bin).
        for k in [1_usize, 3, 7] {
            let f = k as f64 * fs / n as f64;
            assert!(f < f_s, "test bin {k} must lie below f_S = {f_s}");
            let want = modes.green_at_hz(f);
            assert!(
                (spec[k] - want).norm() < 1e-9 * want.norm().max(1e-30),
                "bin {k}"
            );
        }
        for k in [n / 3, n / 2 - 10] {
            let f = k as f64 * fs / n as f64;
            assert!(f > f_s, "test bin {k} must lie above f_S");
            let want = r.transfer_function(src, rcv, &[f], 3, &air).unwrap()[0];
            assert!(
                (spec[k] - want).norm() < 1e-9 * want.norm().max(1e-30),
                "bin {k}"
            );
        }
    }

    #[test]
    fn rejects_invalid_sampling() {
        let air = Medium::air_20c();
        assert!(hybrid_broadband_ir(
            &room(),
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(4.0, 3.0, 2.0),
            0.0,
            10,
            1,
            &air
        )
        .is_err());
        assert!(hybrid_broadband_ir(
            &room(),
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(4.0, 3.0, 2.0),
            1000.0,
            0,
            1,
            &air
        )
        .is_err());
    }
}
