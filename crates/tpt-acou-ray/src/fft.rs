//! In-place radix-2 complex FFT, used for octave filtering and convolution.
//!
//! Lengths must be powers of two; [`next_pow2`] gives the padded length.
//! The transform pair is unnormalised forward and `1/N`-scaled inverse, so
//! `ifft(fft(x)) = x`.

use alloc::vec;
use alloc::vec::Vec;

use tpt_math_linalg_complex::Complex;

/// Smallest power of two that is at least `n` (and at least 1).
#[must_use]
pub fn next_pow2(n: usize) -> usize {
    n.max(1).next_power_of_two()
}

/// In-place FFT (`inverse = false`) or inverse FFT with `1/N` scaling
/// (`inverse = true`). `x.len()` must be a power of two; other lengths are
/// left untouched.
pub fn fft_in_place(x: &mut [Complex], inverse: bool) {
    let n = x.len();
    if n < 2 || !n.is_power_of_two() {
        return;
    }
    // Bit-reversal permutation.
    let mut j = 0;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j |= bit;
        if i < j {
            x.swap(i, j);
        }
    }
    let sign = if inverse { 1.0 } else { -1.0 };
    let mut len = 2;
    while len <= n {
        let angle = sign * core::f64::consts::TAU / len as f64;
        let wlen = Complex::new(libm::cos(angle), libm::sin(angle));
        for start in (0..n).step_by(len) {
            let mut w = Complex::new(1.0, 0.0);
            for k in 0..len / 2 {
                let a = x[start + k];
                let b = x[start + k + len / 2] * w;
                x[start + k] = a + b;
                x[start + k + len / 2] = a - b;
                w = w * wlen;
            }
        }
        len <<= 1;
    }
    if inverse {
        let scale = 1.0 / n as f64;
        for v in x.iter_mut() {
            *v = *v * scale;
        }
    }
}

/// Linear convolution of two real sequences, via the FFT. The output has
/// `a.len() + b.len() − 1` samples, or is empty if either input is empty.
#[must_use]
pub fn convolve_real(a: &[f64], b: &[f64]) -> Vec<f64> {
    if a.is_empty() || b.is_empty() {
        return Vec::new();
    }
    let out_len = a.len() + b.len() - 1;
    let n = next_pow2(out_len);
    let mut fa = vec![Complex::new(0.0, 0.0); n];
    let mut fb = vec![Complex::new(0.0, 0.0); n];
    for (slot, &v) in fa.iter_mut().zip(a) {
        *slot = Complex::new(v, 0.0);
    }
    for (slot, &v) in fb.iter_mut().zip(b) {
        *slot = Complex::new(v, 0.0);
    }
    fft_in_place(&mut fa, false);
    fft_in_place(&mut fb, false);
    for (x, y) in fa.iter_mut().zip(&fb) {
        *x = *x * *y;
    }
    fft_in_place(&mut fa, true);
    fa.iter().take(out_len).map(|c| c.re).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn c(re: f64) -> Complex {
        Complex::new(re, 0.0)
    }

    #[test]
    fn round_trip_recovers_the_input() {
        let original: Vec<Complex> = (0..16).map(|i| c((i as f64 * 0.7).sin())).collect();
        let mut x = original.clone();
        fft_in_place(&mut x, false);
        fft_in_place(&mut x, true);
        for (a, b) in x.iter().zip(&original) {
            assert!((*a - *b).norm() < 1e-12);
        }
    }

    #[test]
    fn single_tone_lands_in_its_bin() {
        let n = 32;
        let k = 5;
        let mut x: Vec<Complex> = (0..n)
            .map(|i| {
                c(libm::cos(
                    core::f64::consts::TAU * k as f64 * i as f64 / n as f64,
                ))
            })
            .collect();
        fft_in_place(&mut x, false);
        assert!((x[k].norm() - n as f64 / 2.0).abs() < 1e-9);
        assert!(x[k + 1].norm() < 1e-9);
    }

    #[test]
    fn convolution_with_a_delta_is_identity() {
        let a = [1.0, 2.0, 3.0];
        let delta = [0.0, 1.0];
        let y = convolve_real(&a, &delta);
        assert_eq!(y.len(), 4);
        for (got, want) in y.iter().zip([0.0, 1.0, 2.0, 3.0]) {
            assert!((got - want).abs() < 1e-12, "{y:?}");
        }
    }

    #[test]
    fn convolution_matches_the_direct_sum() {
        let a: Vec<f64> = (0..7).map(|i| (i as f64).cos()).collect();
        let b: Vec<f64> = (0..5).map(|i| 1.0 / (1.0 + i as f64)).collect();
        let fast = convolve_real(&a, &b);
        for (n, got) in fast.iter().enumerate() {
            let mut direct = 0.0;
            for (i, &ai) in a.iter().enumerate() {
                if n >= i && n - i < b.len() {
                    direct += ai * b[n - i];
                }
            }
            assert!((got - direct).abs() < 1e-10, "n={n}: {got} vs {direct}");
        }
    }

    #[test]
    fn next_pow2_rounds_up() {
        assert_eq!(next_pow2(0), 1);
        assert_eq!(next_pow2(5), 8);
        assert_eq!(next_pow2(8), 8);
    }
}
