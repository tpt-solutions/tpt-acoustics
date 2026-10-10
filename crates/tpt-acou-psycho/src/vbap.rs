//! Vector Base Amplitude Panning (Ville Pulkki, 1997).
//!
//! A source direction is panned across the loudspeaker *base* (a pair in 2-D,
//! a triplet in 3-D) that contains it: the gains come from expressing the
//! source direction in the inverted base of unit speaker vectors,
//! `g = B^{-1} p`, normalized so `|g| = 1` and clipped at negative gains.

use alloc::vec;
use alloc::vec::Vec;

use tpt_acou_core::error::AcousticError;
use tpt_acou_core::Vec3;
use tpt_math_linalg_dense::{DMatrix, DVector};

/// A 2-D or 3-D loudspeaker layout with precomputed bases.
#[derive(Debug, Clone)]
pub struct VbapLayout {
    speakers: Vec<Vec3>,
    /// 3-D triplets (or 2-D pairs) as speaker indices per base.
    bases: Vec<[usize; 3]>,
    is_2d: bool,
}

impl VbapLayout {
    /// Builds a 3-D layout from unit speaker directions, computing the
    /// valid base triplets (speakers on a common great-circle arc produce
    /// degenerate bases and are skipped).
    ///
    /// # Errors
    /// Returns [`AcousticError::InvalidArgument`] for non-unit directions.
    pub fn new_3d(speakers: Vec<Vec3>) -> Result<Self, AcousticError> {
        for s in &speakers {
            if (s.norm() - 1.0).abs() > 1e-6 {
                return Err(AcousticError::InvalidArgument(
                    "speaker directions must be unit vectors",
                ));
            }
        }
        let n = speakers.len();
        let mut bases = Vec::new();
        for i in 0..n {
            for j in (i + 1)..n {
                for k in (j + 1)..n {
                    // Base must positively span the 3-D space.
                    let b = [speakers[i], speakers[j], speakers[k]];
                    let det = b[0].dot(b[1].cross(b[2]));
                    if det.abs() < 1e-9 {
                        continue;
                    }
                    // Only bases whose Voronoi region contains no other
                    // speaker (cheap check: all speakers inside the base's
                    // spherical triangle are exactly these three).
                    if is_valid_base(&speakers, i, j, k) {
                        bases.push([i, j, k]);
                    }
                }
            }
        }
        Ok(Self {
            speakers,
            bases,
            is_2d: false,
        })
    }

    /// Builds a 2-D layout from unit speaker directions in the xy plane
    /// (z ignored), with adjacent-pair bases.
    ///
    /// # Errors
    /// Returns [`AcousticError::InvalidArgument`] for non-unit directions.
    pub fn new_2d(mut speakers: Vec<Vec3>) -> Result<Self, AcousticError> {
        for s in &speakers {
            if (s.norm() - 1.0).abs() > 1e-6 {
                return Err(AcousticError::InvalidArgument(
                    "speaker directions must be unit vectors",
                ));
            }
        }
        speakers.sort_by(|a, b| {
            let aa = libm::atan2(a.y, a.x);
            let bb = libm::atan2(b.y, b.x);
            aa.partial_cmp(&bb).unwrap_or(core::cmp::Ordering::Equal)
        });
        let n = speakers.len();
        let mut bases = Vec::new();
        for i in 0..n {
            bases.push([i, (i + 1) % n, 0]);
        }
        Ok(Self {
            speakers,
            bases,
            is_2d: true,
        })
    }

    /// Speaker directions.
    #[must_use]
    pub fn speakers(&self) -> &[Vec3] {
        &self.speakers
    }

    /// VBAP gains for a unit source direction: one gain per speaker, largest
    /// first-come normalization `|g| = 1` (Pulkki's convention), zero for
    /// speakers outside the active base.
    #[must_use]
    pub fn gains(&self, direction: Vec3) -> Vec<f64> {
        let d = direction.normalized();
        let mut best_gains = vec![0.0_f64; self.speakers.len()];
        let mut best_norm = -1.0_f64;
        for base in &self.bases {
            let (g0, g1, g2) = if self.is_2d {
                let (a, b) = (self.speakers[base[0]], self.speakers[base[1]]);
                // Solve [a b] [g0 g1]^T = d in 2-D.
                let det = a.x * b.y - a.y * b.x;
                if det.abs() < 1e-12 {
                    continue;
                }
                (
                    (d.x * b.y - d.y * b.x) / det,
                    (a.x * d.y - a.y * d.x) / det,
                    0.0,
                )
            } else {
                let bmat = DMatrix::from_vec(
                    3,
                    3,
                    vec![
                        self.speakers[base[0]].x,
                        self.speakers[base[0]].y,
                        self.speakers[base[0]].z,
                        self.speakers[base[1]].x,
                        self.speakers[base[1]].y,
                        self.speakers[base[1]].z,
                        self.speakers[base[2]].x,
                        self.speakers[base[2]].y,
                        self.speakers[base[2]].z,
                    ],
                );
                let dv = DVector::from_vec(vec![d.x, d.y, d.z]);
                let Some(g) = (|| {
                    let sol = bmat.solve(&dv).ok()?;
                    Some((sol[0], sol[1], sol[2]))
                })() else {
                    continue;
                };
                g
            };
            if g0 < -1e-9 || g1 < -1e-9 || g2 < -1e-9 {
                continue; // source outside this base
            }
            let norm = libm::sqrt(g0 * g0 + g1 * g1 + g2 * g2);
            if norm > best_norm && norm > 1e-12 {
                best_norm = norm;
                best_gains = vec![0.0; self.speakers.len()];
                best_gains[base[0]] = g0 / norm;
                best_gains[base[1]] = g1 / norm;
                if !self.is_2d {
                    best_gains[base[2]] = g2 / norm;
                }
            }
        }
        best_gains
    }
}

/// Cheap validity check: no other speaker lies inside the spherical triangle
/// spanned by the base (dot products within the cone bounds).
fn is_valid_base(speakers: &[Vec3], i: usize, j: usize, k: usize) -> bool {
    let (a, b, c) = (speakers[i], speakers[j], speakers[k]);
    for (idx, s) in speakers.iter().enumerate() {
        if idx == i || idx == j || idx == k {
            continue;
        }
        // Inside the spherical triangle (plain cone intersection test).
        let inside_ab = b.cross(a).dot(*s) < 0.0;
        let inside_bc = c.cross(b).dot(*s) < 0.0;
        let inside_ca = a.cross(c).dot(*s) < 0.0;
        if inside_ab && inside_bc && inside_ca {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn stereo_panning_on_speaker_gives_unity() {
        let layout = VbapLayout::new_2d(alloc::vec![
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ])
        .unwrap();
        let g = layout.gains(Vec3::new(1.0, 0.0, 0.0));
        assert!(approx(g[0], 1.0, 1e-9), "g {g:?}");
        assert!(g[1].abs() < 1e-9);
    }

    #[test]
    fn stereo_panning_midway_is_equal_power_pair() {
        let layout = VbapLayout::new_2d(alloc::vec![
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ])
        .unwrap();
        let d = Vec3::new(1.0, 1.0, 0.0).normalized();
        let g = layout.gains(d);
        // Equal gains on both speakers, normalized |g| = 1.
        assert!(approx(g[0], g[1], 1e-9));
        assert!(approx(g[0], core::f64::consts::FRAC_1_SQRT_2, 1e-9));
    }

    #[test]
    fn threedee_triplet_on_speaker() {
        // Octahedral layout: ±x, ±y, ±z.
        let layout = VbapLayout::new_3d(alloc::vec![
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(-1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, -1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(0.0, 0.0, -1.0),
        ])
        .unwrap();
        // Source straight up +z: gains should concentrate on +z (and its
        // triplet partners sharing energy), never on -z.
        let g = layout.gains(Vec3::new(0.0, 0.0, 1.0));
        assert!(g[5].abs() < 1e-12, "g {g:?}");
        assert!(g[4] > 0.0 || g[0] + g[2] > 0.0);
    }

    #[test]
    fn gains_are_normalized() {
        let layout = VbapLayout::new_3d(alloc::vec![
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
            Vec3::new(-1.0, -1.0, -1.0).normalized(),
        ])
        .unwrap();
        for dir in [
            Vec3::new(0.5, 0.5, 0.5),
            Vec3::new(0.8, 0.1, 0.3),
            Vec3::new(-0.5, -0.6, -0.7),
        ] {
            let g = layout.gains(dir.normalized());
            let n: f64 = libm::sqrt(g.iter().map(|x| x * x).sum::<f64>());
            assert!(approx(n, 1.0, 1e-9) || n < 1e-12, "|g| = {n}");
        }
    }

    #[test]
    fn rejects_nonunit_speakers() {
        assert!(VbapLayout::new_2d(alloc::vec![Vec3::new(2.0, 0.0, 0.0)]).is_err());
        assert!(VbapLayout::new_3d(alloc::vec![Vec3::new(2.0, 0.0, 0.0)]).is_err());
    }
}
