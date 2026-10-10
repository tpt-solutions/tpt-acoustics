//! Head-Related Transfer Functions (HRTFs) and interaural time differences.
//!
//! * [`spherical_head_itd`] — Woodworth's formula: the interaural time
//!   difference for a rigid sphere of radius `a` (the "snowman head" model),
//!   `ITD(theta) = a/c * (theta + sin(theta))` for the azimuth `theta` in
//!   radians (0 = front, pi/2 = side).
//! * [`HrtfGrid`] — a set of HRIRs on a direction grid with
//!   nearest-neighbour + spherical-interpolation retrieval.
//!
//! Pinna spectral cues are carried by the measured HRIRs; this module does
//! not synthesize them.

use alloc::vec;
use alloc::vec::Vec;

use tpt_acou_core::error::AcousticError;
use tpt_acou_core::Vec3;

/// Interaural time difference (seconds) of the rigid-sphere (Woodworth)
/// model for azimuth `theta` (0 = front, ±pi = back) at ear distance
/// `radius` and sound speed `c`.
#[must_use]
pub fn spherical_head_itd(theta: f64, radius: f64, c: f64) -> f64 {
    radius / c * (theta + libm::sin(theta))
}

/// One measured head-related impulse response plus its direction.
#[derive(Debug, Clone)]
pub struct Hrir {
    /// Unit direction vector from the head center to the source.
    pub direction: Vec3,
    /// Impulse response samples for the left ear.
    pub left: Vec<f64>,
    /// Impulse response samples for the right ear.
    pub right: Vec<f64>,
    /// Sample rate (Hz).
    pub sample_rate: f64,
}

/// A collection of HRIRs with nearest-direction retrieval.
#[derive(Debug, Clone, Default)]
pub struct HrtfGrid {
    hrirs: Vec<Hrir>,
}

impl HrtfGrid {
    /// Builds a grid from measurements.
    ///
    /// # Errors
    /// Returns [`AcousticError::InvalidArgument`] if any direction is not
    /// (approximately) unit length or any impulse-response pair has
    /// different lengths.
    pub fn new(hrirs: Vec<Hrir>) -> Result<Self, AcousticError> {
        for h in &hrirs {
            if (h.direction.norm() - 1.0).abs() > 1e-6 {
                return Err(AcousticError::InvalidArgument(
                    "HRIR directions must be unit vectors",
                ));
            }
            if h.left.len() != h.right.len() {
                return Err(AcousticError::InvalidArgument(
                    "left/right HRIRs must have equal length",
                ));
            }
        }
        Ok(Self { hrirs })
    }

    /// Number of measurement directions.
    #[must_use]
    pub fn len(&self) -> usize {
        self.hrirs.len()
    }

    /// True if the grid is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.hrirs.is_empty()
    }

    /// Nearest measured direction's HRIRs (linear search; grids are small).
    #[must_use]
    pub fn nearest(&self, direction: Vec3) -> Option<&Hrir> {
        let d = direction.normalized();
        self.hrirs.iter().min_by(|a, b| {
            let da = 1.0 - a.direction.dot(d);
            let db = 1.0 - b.direction.dot(d);
            da.partial_cmp(&db).unwrap_or(core::cmp::Ordering::Equal)
        })
    }

    /// Left/right impulse responses interpolated between the two nearest
    /// measurements (`HRIR = (1-t)*A + t*B` weighted by the angle ratio —
    /// a first-order directional interpolation adequate for grids with
    /// ~15-degree spacing).
    #[must_use]
    pub fn interpolated(&self, direction: Vec3) -> Option<(Vec<f64>, Vec<f64>)> {
        let d = direction.normalized();
        if self.hrirs.is_empty() {
            return None;
        }
        // Two nearest.
        let mut best = [0_usize, 0_usize];
        let mut best_d = [f64::INFINITY, f64::INFINITY];
        for (i, h) in self.hrirs.iter().enumerate() {
            let dist = 1.0 - h.direction.dot(d);
            if dist < best_d[0] {
                best_d[1] = best_d[0];
                best[1] = best[0];
                best_d[0] = dist;
                best[0] = i;
            } else if dist < best_d[1] {
                best_d[1] = dist;
                best[1] = i;
            }
        }
        let (a, b) = (&self.hrirs[best[0]], &self.hrirs[best[1]]);
        let cos_a = a.direction.dot(d);
        let cos_b = b.direction.dot(d);
        let t = if cos_a + cos_b <= 1e-12 {
            0.5
        } else {
            (cos_b / (cos_a + cos_b)).clamp(0.0, 1.0)
        };
        let n = a.left.len();
        let mut left = vec![0.0_f64; n];
        let mut right = vec![0.0_f64; n];
        for i in 0..n {
            left[i] = (1.0 - t) * a.left[i] + t * b.left[i];
            right[i] = (1.0 - t) * a.right[i] + t * b.right[i];
        }
        Some((left, right))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn woodworth_itd_at_side_is_radius_over_c() {
        let a = 0.0875; // 8.75 cm head radius
        let c = 343.0;
        // theta = pi/2: ITD = a/c (pi/2 + 1) ~ 650 microseconds for a
        // typical head radius.
        let itd = spherical_head_itd(core::f64::consts::FRAC_PI_2, a, c);
        assert!((itd - a / c * (core::f64::consts::FRAC_PI_2 + 1.0)).abs() < 1e-12);
        // Front: zero. (The arc model gives (a/c)*pi at the rear — a known
        // quirk of the Woodworth approximation, kept honest here.)
        assert!(spherical_head_itd(0.0, a, c).abs() < 1e-15);
        assert!(
            (spherical_head_itd(core::f64::consts::PI, a, c) - core::f64::consts::PI * a / c).abs()
                < 1e-12
        );
    }

    #[test]
    fn grid_nearest_and_interpolation() {
        let mk = |dir: Vec3, gain: f64| Hrir {
            direction: dir,
            left: alloc::vec![gain; 8],
            right: alloc::vec![gain * 2.0; 8],
            sample_rate: 44100.0,
        };
        let grid = HrtfGrid::new(alloc::vec![
            mk(Vec3::new(0.0, 0.0, 1.0), 1.0), // front
            mk(Vec3::new(1.0, 0.0, 0.0), 0.5), // right
        ])
        .unwrap();
        assert_eq!(grid.len(), 2);
        // Exactly at a measurement: that measurement.
        let (l, r) = grid.interpolated(Vec3::new(0.0, 0.0, 1.0)).unwrap();
        assert!((l[0] - 1.0).abs() < 1e-12 && (r[0] - 2.0).abs() < 1e-12);
        // Midway (45 degrees): blend ~0.75/0.25-ish, bounded by the ends.
        let diag = (Vec3::new(1.0, 0.0, 1.0)).normalized();
        let (l, _) = grid.interpolated(diag).unwrap();
        assert!(l[0] > 0.5 && l[0] < 1.0, "l[0] = {}", l[0]);
        // Nearest match for the exact direction.
        assert!((grid.nearest(Vec3::new(1.0, 0.0, 0.0)).unwrap().left[0] - 0.5).abs() < 1e-12);
    }

    #[test]
    fn rejects_bad_measurements() {
        let bad = Hrir {
            direction: Vec3::new(2.0, 0.0, 0.0),
            left: alloc::vec![0.0; 4],
            right: alloc::vec![0.0; 5],
            sample_rate: 44100.0,
        };
        assert!(HrtfGrid::new(alloc::vec![bad]).is_err());
    }
}
