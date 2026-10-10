//! Frequency-dependent absorption and scattering coefficients at boundaries.
//!
//! Coefficients are stored as `(frequency_hz, value)` breakpoints sorted by
//! frequency; evaluation between breakpoints is linear in `log(frequency)`
//! (linear in level), constant outside the tabulated range.

use alloc::vec::Vec;

/// Boundary surface behaviour: absorption `alpha(f)` in [0, 1] (fraction of
/// incident *energy* not reflected) and scattering `s(f)` in [0, 1] (fraction
/// of reflected energy redirected diffusely, Lambert-distributed).
#[derive(Debug, Clone, PartialEq)]
pub struct BoundaryMaterial {
    name: &'static str,
    absorption: Vec<(f64, f64)>,
    scattering: Vec<(f64, f64)>,
}

impl BoundaryMaterial {
    /// Builds a material from absorption and scattering breakpoint tables.
    ///
    /// # Errors
    /// Returns [`tpt_acou_core::AcousticError::InvalidArgument`] if a table is
    /// empty, unsorted, or contains values outside [0, 1].
    pub fn new(
        name: &'static str,
        absorption: Vec<(f64, f64)>,
        scattering: Vec<(f64, f64)>,
    ) -> Result<Self, tpt_acou_core::AcousticError> {
        Self::check_table(&absorption, "absorption")?;
        Self::check_table(&scattering, "scattering")?;
        Ok(Self {
            name,
            absorption,
            scattering,
        })
    }

    fn check_table(
        table: &[(f64, f64)],
        what: &'static str,
    ) -> Result<(), tpt_acou_core::AcousticError> {
        use tpt_acou_core::AcousticError;
        if table.is_empty() {
            return Err(AcousticError::InvalidArgument(what));
        }
        if table.windows(2).any(|w| w[0].0 >= w[1].0) {
            return Err(AcousticError::InvalidArgument(what));
        }
        if table.iter().any(|&(_, v)| !(0.0..=1.0).contains(&v)) {
            return Err(AcousticError::InvalidArgument(what));
        }
        Ok(())
    }

    /// The material name.
    #[must_use]
    pub const fn name(&self) -> &'static str {
        self.name
    }

    /// Interpolated absorption coefficient `alpha(f)` (clamped to the
    /// tabulated range).
    #[must_use]
    pub fn absorption_at(&self, freq_hz: f64) -> f64 {
        interp(&self.absorption, freq_hz)
    }

    /// Interpolated scattering coefficient `s(f)` (clamped to the tabulated
    /// range).
    #[must_use]
    pub fn scattering_at(&self, freq_hz: f64) -> f64 {
        interp(&self.scattering, freq_hz)
    }

    /// Amplitude reflection factor magnitude `|r| = sqrt(1 − alpha)` for a
    /// real reflection coefficient (energy conservation: `|r|² = 1 − alpha`).
    #[must_use]
    pub fn amplitude_reflection(&self, freq_hz: f64) -> f64 {
        libm::sqrt((1.0 - self.absorption_at(freq_hz)).max(0.0))
    }

    /// Hard, rigid wall: no absorption, no scattering (frequency-independent).
    #[must_use]
    pub fn rigid() -> Self {
        Self {
            name: "rigid",
            absorption: alloc::vec![(20.0, 0.0), (20_000.0, 0.0)],
            scattering: alloc::vec![(20.0, 0.0), (20_000.0, 0.0)],
        }
    }

    /// Painted concrete: low absorption rising slightly at HF, mild
    /// scattering. Representative mid-band values from common auditorium
    /// tables (not a certified measurement).
    #[must_use]
    pub fn concrete_painted() -> Self {
        Self {
            name: "concrete_painted",
            absorption: alloc::vec![
                (125.0, 0.01),
                (250.0, 0.01),
                (500.0, 0.02),
                (1000.0, 0.02),
                (2000.0, 0.02),
                (4000.0, 0.03)
            ],
            scattering: alloc::vec![(125.0, 0.02), (4000.0, 0.08)],
        }
    }

    /// Heavy curtain draped (pleated) against a wall: strong HF absorption.
    #[must_use]
    pub fn curtain_heavy() -> Self {
        Self {
            name: "curtain_heavy",
            absorption: alloc::vec![
                (125.0, 0.07),
                (250.0, 0.31),
                (500.0, 0.49),
                (1000.0, 0.75),
                (2000.0, 0.70),
                (4000.0, 0.60)
            ],
            scattering: alloc::vec![(125.0, 0.15), (4000.0, 0.35)],
        }
    }

    /// Single pane of window glass.
    #[must_use]
    pub fn glass() -> Self {
        Self {
            name: "glass",
            absorption: alloc::vec![
                (125.0, 0.35),
                (250.0, 0.25),
                (500.0, 0.18),
                (1000.0, 0.12),
                (2000.0, 0.07),
                (4000.0, 0.04)
            ],
            scattering: alloc::vec![(125.0, 0.01), (4000.0, 0.02)],
        }
    }

    /// Typical absorber panel (porous, ~50 mm on wall).
    #[must_use]
    pub fn absorber_panel() -> Self {
        Self {
            name: "absorber_panel",
            absorption: alloc::vec![
                (125.0, 0.25),
                (250.0, 0.60),
                (500.0, 0.90),
                (1000.0, 0.95),
                (2000.0, 0.95),
                (4000.0, 0.90)
            ],
            scattering: alloc::vec![(125.0, 0.30), (4000.0, 0.50)],
        }
    }
}

/// Log-frequency linear interpolation, clamped at both ends.
fn interp(table: &[(f64, f64)], freq_hz: f64) -> f64 {
    debug_assert!(!table.is_empty());
    let f = freq_hz.max(f64::MIN_POSITIVE);
    if f <= table[0].0 {
        return table[0].1;
    }
    let last = table.len() - 1;
    if f >= table[last].0 {
        return table[last].1;
    }
    // Find the bracketing pair (linear scan; tables are tiny).
    let mut i = 0;
    while i + 1 < table.len() && table[i + 1].0 < f {
        i += 1;
    }
    let (f0, v0) = table[i];
    let (f1, v1) = table[i + 1];
    let t = (libm::log(f) - libm::log(f0)) / (libm::log(f1) - libm::log(f0));
    v0 + (v1 - v0) * t
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn interpolation_between_and_beyond_breakpoints() {
        let m = BoundaryMaterial::glass();
        // At a breakpoint: exact value.
        assert!((m.absorption_at(500.0) - 0.18).abs() < 1e-12);
        // Below the first breakpoint: clamped.
        assert!((m.absorption_at(10.0) - 0.35).abs() < 1e-12);
        // Above the last: clamped.
        assert!((m.absorption_at(30_000.0) - 0.04).abs() < 1e-12);
        // Between 500 and 1000 the geometric midpoint maps to t = 0.5 in
        // log-frequency, i.e. the arithmetic mean of the two values.
        let mid = libm::sqrt(500.0 * 1000.0);
        let a = m.absorption_at(mid);
        assert!((a - 0.5 * (0.18 + 0.12)).abs() < 1e-12);
    }

    #[test]
    fn amplitude_reflection_conserves_energy() {
        let m = BoundaryMaterial::absorber_panel();
        let r = m.amplitude_reflection(1000.0);
        let alpha = m.absorption_at(1000.0);
        assert!((r * r - (1.0 - alpha)).abs() < 1e-12);
    }

    #[test]
    fn rigid_has_no_losses() {
        let m = BoundaryMaterial::rigid();
        assert_eq!(m.absorption_at(1000.0), 0.0);
        assert_eq!(m.scattering_at(1000.0), 0.0);
        assert!((m.amplitude_reflection(1000.0) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn rejects_bad_tables() {
        use alloc::vec;
        assert!(BoundaryMaterial::new("x", vec![], vec![(100.0, 0.0)]).is_err());
        assert!(
            BoundaryMaterial::new("x", vec![(1000.0, 0.5), (100.0, 0.2)], vec![(100.0, 0.0)])
                .is_err()
        );
        assert!(BoundaryMaterial::new("x", vec![(100.0, 1.5)], vec![(100.0, 0.0)]).is_err());
        assert!(BoundaryMaterial::new("x", vec![(100.0, 0.2)], vec![(100.0, 0.0)]).is_ok());
    }
}
