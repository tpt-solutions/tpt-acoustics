//! Ambisonics: SN3D B-format encoding to 2nd order and velocity decoding.
//!
//! Channel ordering is ACN/N3D index order `(#, index)` with **SN3D**
//! normalization (max-N with the N3D factors dropped), the AmbiX convention:
//! channels `W, Y, Z, X, V2, V3, ...` carry the real spherical harmonics
//! `Y_nm` (SN3D), with `W = Y_00 = 1` (not 1/sqrt(2) — that is FuMa).
//!
//! Decoding uses the **velocity** (least-squares) decoder: speaker gains
//! solve `A g = h` direction-matched in the least-squares sense, computed
//! with the real dense solver from `tpt-math-linalg-dense`.

use alloc::vec;
use alloc::vec::Vec;

use tpt_acou_core::error::AcousticError;
use tpt_acou_core::Vec3;
use tpt_math_linalg_dense::DMatrix;

/// Number of channels for 2nd-order ambisonics (ACN 0..8).
pub const CHANNELS_ORDER_2: usize = 9;

/// Real spherical harmonic `Y_nm` (SN3D normalization) for direction
/// `d` (unit vector). Returns the value and the flat ACN channel index
/// `(n, m)` with `n` the order and `m in [-n, n]`.
///
/// Channel index = `n*n + n + m` (ACN).
#[must_use]
pub fn sn3d_harmonic(n: usize, m: i32, d: Vec3) -> f64 {
    debug_assert!(m.unsigned_abs() as usize <= n);
    let (x, y, z) = (d.x, d.y, d.z);
    match (n, m) {
        (0, 0) => 1.0,
        (1, -1) => y,
        (1, 0) => z,
        (1, 1) => x,
        // SN3D AmbiX set (all channels bounded by ~1 on the unit sphere).
        (2, -2) => libm::sqrt(3.0) * x * y,
        (2, -1) => libm::sqrt(3.0) * y * z,
        (2, 0) => (3.0 * z * z - 1.0) / 2.0,
        (2, 1) => libm::sqrt(3.0) * x * z,
        (2, 2) => (libm::sqrt(3.0) / 2.0) * (x * x - y * y),
        _ => 0.0,
    }
}

/// Encodes a unit-amplitude plane-wave source from direction `d` into a
/// 2nd-order SN3D B-format vector (ACN order).
#[must_use]
pub fn encode_order2(d: Vec3) -> Vec<f64> {
    let mut out = vec![0.0_f64; CHANNELS_ORDER_2];
    for n in 0..=2 {
        for m in -(n as i32)..=(n as i32) {
            let ch = (n * n + n) as i32 + m;
            out[ch as usize] = sn3d_harmonic(n, m, d);
        }
    }
    out
}

/// A loudspeaker direction with unit position.
#[derive(Debug, Clone, Copy)]
pub struct Speaker {
    /// Unit direction to the speaker.
    pub direction: Vec3,
}

impl Speaker {
    /// Creates a speaker.
    ///
    /// # Errors
    /// Returns [`AcousticError::InvalidArgument`] for a non-unit direction.
    pub fn new(direction: Vec3) -> Result<Self, AcousticError> {
        if (direction.norm() - 1.0).abs() > 1e-6 {
            return Err(AcousticError::InvalidArgument(
                "speaker direction must be unit length",
            ));
        }
        Ok(Self { direction })
    }
}

/// Velocity (least-squares) decoder: `g = (A^T A)^{-1} A^T h`, where `A` is
/// the `speakers x channels` harmonic matrix and `h` the encoded source.
/// Rigid (imaginary) modes are not constrained — this is the classic
/// velocity decoder, adequate below the spatial aliasing limit.
///
/// Returns one gain per speaker.
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] for fewer speakers than
/// channels (underdetermined) and [`AcousticError::Singular`] if the
/// speaker layout is degenerate.
pub fn decode_velocity(speakers: &[Speaker], encoded: &[f64]) -> Result<Vec<f64>, AcousticError> {
    let s = speakers.len();
    if s < encoded.len() {
        return Err(AcousticError::InvalidArgument(
            "velocity decoding needs at least one speaker per B-format channel",
        ));
    }
    // Pure harmonic matrix: A[i][ch] = Y_ch(speaker_i); the encoded source
    // enters only through the right-hand side of A^T g = encoded.
    let a = DMatrix::from_fn(s, encoded.len(), |i, ch| {
        sn3d_harmonic(order_of(ch), m_of(ch), speakers[i].direction)
    });
    // Normal equations: (A^T A) g = A^T h_vec — but here "h" is implicit:
    // the velocity decoder matches the SUM over speakers to the encoded
    // vector: minimize ||A^T g - encoded^T||_2 (reproducing the encoded
    // direction from the speaker array). Solve (A A^T) is wrong; the
    // standard form: g = A (A^T A)^{-1} encoded? Use g = (A^T)^+ encoded:
    // least-norm solution of A^T g = encoded.
    let ata_inv = pseudo_inverse_ata(&a)?;
    // g = A (A^T A)^{-1} encoded: (A^T g) matches encoded exactly when
    // determined.
    let mut g = vec![0.0_f64; s];
    for i in 0..s {
        let mut acc = 0.0;
        for ch in 0..encoded.len() {
            // (A^T A)^{-1} e : row ch of the inverse times e.
            let mut row_dot_e = 0.0;
            for k in 0..encoded.len() {
                row_dot_e += ata_inv[ch][k] * encoded[k];
            }
            acc += sn3d_harmonic(order_of(ch), m_of(ch), speakers[i].direction) * row_dot_e;
        }
        g[i] = acc;
    }
    Ok(g)
}

fn order_of(ch: usize) -> usize {
    libm::sqrt(ch as f64) as usize
}

fn m_of(ch: usize) -> i32 {
    let n = order_of(ch);
    ch as i32 - (n * n + n) as i32
}

/// `(A^T A)^{-1}` via the dense real solver (LU) from `tpt-math-linalg-dense`.
fn pseudo_inverse_ata(a: &DMatrix) -> Result<Vec<Vec<f64>>, AcousticError> {
    let cols = a.ncols();
    let ata = DMatrix::from_fn(cols, cols, |i, j| {
        (0..a.nrows()).map(|r| a[(r, i)] * a[(r, j)]).sum::<f64>()
    });
    let mut inv = vec![vec![0.0_f64; cols]; cols];
    for col in 0..cols {
        let mut e = vec![0.0_f64; cols];
        e[col] = 1.0;
        let x = ata
            .solve(&tpt_math_linalg_dense::DVector::from_vec(e))
            .map_err(|_| AcousticError::Singular("degenerate speaker layout"))?;
        for r in 0..cols {
            inv[r][col] = x[r];
        }
    }
    Ok(inv)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn w_channel_is_one() {
        assert!(approx(
            encode_order2(Vec3::new(0.0, 0.0, 1.0))[0],
            1.0,
            1e-12
        ));
    }

    #[test]
    fn first_order_channels_match_direction() {
        // Source at +x: X channel (ACN 3) = 1, Y (ACN 1) = 0.
        let e = encode_order2(Vec3::new(1.0, 0.0, 0.0));
        assert!(approx(e[1], 0.0, 1e-12)); // Y
        assert!(approx(e[3], 1.0, 1e-12)); // X
        assert!(approx(e[2], 0.0, 1e-12)); // Z
                                           // Source at +z: Z channel = 1.
        let e = encode_order2(Vec3::new(0.0, 0.0, 1.0));
        assert!(approx(e[2], 1.0, 1e-12));
    }

    #[test]
    fn second_order_channels_match_ambix_table() {
        // ACN4 (2,-2) at the 45-degree horizontal diagonal: sqrt(3)*xy = sqrt(3)/2.
        let d45 = Vec3::new(1.0, 1.0, 0.0).normalized();
        let e = encode_order2(d45);
        assert!(approx(e[4], libm::sqrt(3.0) / 2.0, 1e-9), "ACN4 = {}", e[4]);
        // ACN8 (2,2) vanishes on that diagonal (x^2 = y^2).
        assert!(approx(e[8], 0.0, 1e-12));
        // ACN6 (2,0) at zenith: (3-1)/2 = 1.
        let e = encode_order2(Vec3::new(0.0, 0.0, 1.0));
        assert!(approx(e[6], 1.0, 1e-12));
    }

    #[test]
    fn velocity_decoder_pins_a_matching_speaker() {
        // A tetrahedral-ish layout: the decoder should reproduce the encoded
        // direction through the array (A^T g = h exactly, determined case
        // with 9 speakers for 9 channels is heavy; use 1st-order content on
        // 4 speakers via a truncated encode).
        let sp = [
            Speaker::new(Vec3::new(1.0, 0.0, 0.0)).unwrap(),
            Speaker::new(Vec3::new(-1.0, 0.0, 0.0)).unwrap(),
            Speaker::new(Vec3::new(0.0, 1.0, 0.0)).unwrap(),
            Speaker::new(Vec3::new(0.0, -1.0, 0.0)).unwrap(),
            Speaker::new(Vec3::new(0.0, 0.0, 1.0)).unwrap(),
            Speaker::new(Vec3::new(0.0, 0.0, -1.0)).unwrap(),
        ];
        // Encode only the 1st-order channels (ACN 0..3).
        let src = Vec3::new(0.6, 0.0, 0.8).normalized();
        let mut encoded = vec![1.0_f64, src.y, src.z, src.x];
        // Velocity decoder on 6 speakers, 4 channels: g should activate the
        // +z and +x speakers most.
        let g = decode_velocity(&sp, &encoded).unwrap();
        assert_eq!(g.len(), 6);
        let pos_total = g[0] + g[4]; // +x, +z speakers
        let neg_total = g[1] + g[2] + g[3] + g[5];
        assert!(pos_total > neg_total, "gains {g:?}");
        encoded[0] = encoded[0].max(1.0); // keep W normalized
    }

    #[test]
    fn rejects_insufficient_speakers() {
        let sp = [Speaker::new(Vec3::new(1.0, 0.0, 0.0)).unwrap()];
        assert!(decode_velocity(&sp, &encode_order2(Vec3::new(1.0, 0.0, 0.0))).is_err());
    }

    #[test]
    fn rejects_bad_speaker_direction() {
        assert!(Speaker::new(Vec3::new(2.0, 0.0, 0.0)).is_err());
    }
}
