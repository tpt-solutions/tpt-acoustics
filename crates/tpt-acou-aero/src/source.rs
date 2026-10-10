//! Source-side data types for aeroacoustic kernels.
//!
//! These are plain numeric containers (no CFD crate dependency): the
//! `tpt-physics` CFD adapter will *fill* these from unsteady flow fields once
//! that integration lands (todo.md "Integration & Synergies").

use alloc::vec::Vec;

use tpt_acou_core::Vec3;

/// A time history sampled at a uniform rate.
#[derive(Debug, Clone, Default)]
pub struct TimeSignal {
    /// Sample values.
    pub samples: Vec<f64>,
    /// Sample rate (Hz).
    pub sample_rate: f64,
}

impl TimeSignal {
    /// Builds a signal from samples and a sample rate.
    ///
    /// # Errors
    /// Returns [`tpt_acou_core::AcousticError::NonPositive`] for a
    /// non-positive sample rate.
    pub fn new(samples: Vec<f64>, sample_rate: f64) -> Result<Self, tpt_acou_core::AcousticError> {
        if sample_rate <= 0.0 {
            return Err(tpt_acou_core::AcousticError::NonPositive("sample rate"));
        }
        Ok(Self {
            samples,
            sample_rate,
        })
    }

    /// Time at sample index `i`.
    #[must_use]
    pub fn time(&self, i: usize) -> f64 {
        i as f64 / self.sample_rate
    }

    /// Linear interpolation of the signal at time `t` (zero outside the
    /// sampled window).
    #[must_use]
    pub fn at(&self, t: f64) -> f64 {
        let x = t * self.sample_rate;
        let i = libm::floor(x) as usize;
        if i + 1 >= self.samples.len() {
            return 0.0;
        }
        let f = x - i as f64;
        self.samples[i] * (1.0 - f) + self.samples[i + 1] * f
    }
}

/// A compact acoustic source whose strength is a scalar time history
/// (monopole: volume acceleration; dipole: force component projection).
#[derive(Debug, Clone)]
pub struct CompactSource {
    /// Source position.
    pub position: Vec3,
    /// Unit directivity vector (for a dipole: the force direction; for a
    /// monopole: unused, set to [`Vec3::ZERO`]).
    pub direction: Vec3,
    /// Strength history: monopole `Q(t)` (volume injection, m³/s) or dipole
    /// `F(t)` (force magnitude along `direction`, N).
    pub strength: TimeSignal,
}

/// Unsteady surface-force history for Curle's dipole: the total force `F(t)`
/// exerted by the fluid on the body (N), with components in world axes.
#[derive(Debug, Clone, Default)]
pub struct UnsteadyForce {
    /// Force components vs time, sampled at `sample_rate`.
    pub fx: TimeSignal,
    /// See [`UnsteadyForce::fx`].
    pub fy: TimeSignal,
    /// See [`UnsteadyForce::fx`].
    pub fz: TimeSignal,
    /// Common sample rate (Hz).
    pub sample_rate: f64,
}

/// One quadrupole cell of a Lighthill source distribution: a compact volume
/// element with an averaged Lighthill-stress tensor `T_ij` history.
#[derive(Debug, Clone)]
pub struct QuadrupoleCell {
    /// Cell centroid.
    pub centroid: Vec3,
    /// Cell volume (m³).
    pub volume: f64,
    /// Averaged Lighthill stress `T_rr` along the radiation direction,
    /// sampled uniformly (Pa).
    pub t_rr: TimeSignal,
}

#[cfg(test)]
mod tests {
    use super::*;
    use alloc::vec;

    #[test]
    fn signal_interpolation() {
        let s = TimeSignal::new(vec![0.0, 1.0, 2.0, 3.0], 100.0).unwrap();
        assert!((s.at(0.005) - 0.5).abs() < 1e-12);
        assert!((s.at(0.02) - 2.0).abs() < 1e-12);
        assert_eq!(s.at(1.0), 0.0); // outside window
        assert!(TimeSignal::new(vec![], -1.0).is_err());
    }
}
