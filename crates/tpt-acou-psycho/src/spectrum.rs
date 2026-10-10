//! FFT-based band levels and weighted levels (on `tpt-dsp-core`).

use alloc::vec;
use alloc::vec::Vec;
use num_complex::Complex;
use tpt_acou_core::quantity::{Pressure, SoundPressureLevel};
use tpt_dsp_core::{fft_inplace, next_power_of_two};

/// One-sided RMS-amplitude spectrum (Pa per bin) of a real signal sampled at
/// `fs`, with a Hann window. The transform length is the next power of two
/// at or above `signal.len()` (zero-padded, per `tpt-dsp-core`'s FFT); the
/// power normalization uses the window's mean square, so the energetic sum
/// over bins conserves the *unwindowed* signal's mean square (ENBW-style).
#[must_use]
pub fn magnitude_spectrum(signal: &[f64], fs: f64) -> (Vec<f64>, Vec<f64>) {
    let n = signal.len();
    if n < 2 || fs <= 0.0 {
        return (vec![], vec![]);
    }
    // Hann window applied in place; track its mean square for the
    // ENBW-style power normalization.
    let mut mean_w2 = 0.0_f64;
    let mut buf: Vec<Complex<f64>> = Vec::with_capacity(n);
    for (i, &x) in signal.iter().enumerate() {
        let w = 0.5 * (1.0 - libm::cos(core::f64::consts::TAU * i as f64 / n as f64));
        mean_w2 += w * w;
        buf.push(Complex::new(x * w, 0.0));
    }
    mean_w2 /= n as f64;

    let np2 = next_power_of_two(n);
    buf.resize(np2, Complex::new(0.0, 0.0));
    let mut scratch = vec![Complex::new(0.0, 0.0); np2];
    fft_inplace(&mut buf, &mut scratch);

    let half = np2 / 2;
    let mut freqs = Vec::with_capacity(half);
    let mut mags = Vec::with_capacity(half);
    for k in 0..half {
        freqs.push(fs * k as f64 / np2 as f64);
        // One-sided bin mean-square: 2|X|^2 / (N * n * mean(w^2)) for
        // interior bins. With zero-padding (N >= n) this keeps the
        // energetic sum over all bins equal to the *unwindowed* signal's
        // mean square for any N/n ratio (Parseval: sum|X|^2 = N*sum|x_w|^2).
        let ms = if k > 0 && k * 2 != np2 {
            2.0 * buf[k].norm_sqr()
        } else {
            buf[k].norm_sqr()
        } / (np2 as f64 * n as f64 * mean_w2);
        mags.push(libm::sqrt(ms));
    }
    (freqs, mags)
}

/// A-weighted sound pressure level of a signal sampled at `fs`: the
/// energetic sum of A-weighted one-sided bin levels.
///
/// Returns `None` for an empty/silent signal.
#[must_use]
pub fn a_weighted_level(signal: &[f64], fs: f64) -> Option<SoundPressureLevel> {
    let (freqs, mags) = magnitude_spectrum(signal, fs);
    let mut levels = Vec::new();
    for (f, m) in freqs.iter().zip(&mags) {
        if *m <= 0.0 {
            continue;
        }
        levels.push(
            SoundPressureLevel::from_pressure_rms(Pressure::new(*m))
                .add_db(crate::weighting::a_weighting_db(*f)),
        );
    }
    SoundPressureLevel::combine(&levels)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sine_level_recovers_rms() {
        // 1 kHz sine at 1 Pa RMS -> ~93.98 dB; the A-weighting at 1 kHz is
        // ~0 dB so the weighted level matches. Zero-padding to a power of
        // two spreads leakage across bins, which the energetic sum folds
        // back in (ENBW-conserving normalization).
        let fs = 48000.0_f64;
        let n = 48000;
        let signal: Vec<f64> = (0..n)
            .map(|i| {
                core::f64::consts::SQRT_2
                    * libm::sin(core::f64::consts::TAU * 1000.0 * i as f64 / fs)
            })
            .collect();
        let lvl = a_weighted_level(&signal, fs).unwrap();
        let want = 20.0 * libm::log10(1.0 / 20.0e-6);
        assert!(
            (lvl.db() - want).abs() < 0.5,
            "got {:.3} want {want:.3}",
            lvl.db()
        );
    }

    #[test]
    fn silent_signal_has_no_level() {
        assert!(a_weighted_level(&[0.0; 1024], 48000.0).is_none());
    }
}
