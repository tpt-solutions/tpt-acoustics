//! Complex-valued acoustic impedance and admittance.
//!
//! The *specific acoustic impedance* `z = p/u` (Pa·s/m, "rayl") is the
//! frequency-domain bridge between pressure and particle velocity at a
//! boundary. Normal-incidence reflection follows
//! `r = (z' − 1)/(z' + 1)` for the impedance normalized by the medium's
//! characteristic impedance `rho·c`, and the absorption coefficient is
//! `alpha = 1 − |r|²`.

use tpt_math_linalg_complex::Complex;

use crate::error::AcousticError;
use crate::medium::Medium;

/// Specific acoustic impedance `z = p/u` in Pa·s/m (complex, frequency
/// dependent).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpecificImpedance(pub Complex);

/// Specific acoustic admittance `y = u/p = 1/z` in m/(Pa·s).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SpecificAdmittance(pub Complex);

impl SpecificImpedance {
    /// Builds an impedance from a complex value in Pa·s/m.
    ///
    /// # Errors
    /// Returns [`AcousticError::InvalidArgument`] for the zero impedance
    /// (no admittance exists).
    pub fn new(z: Complex) -> Result<Self, AcousticError> {
        if z.re == 0.0 && z.im == 0.0 {
            return Err(AcousticError::InvalidArgument("impedance must be nonzero"));
        }
        Ok(Self(z))
    }

    /// Builds a purely real impedance in Pa·s/m.
    ///
    /// # Errors
    /// See [`SpecificImpedance::new`].
    pub fn real(re: f64) -> Result<Self, AcousticError> {
        Self::new(Complex::new(re, 0.0))
    }

    /// The perfectly rigid (hard-wall) limit: `z -> infinity`.
    #[must_use]
    pub fn rigid() -> Self {
        Self(Complex::new(f64::INFINITY, 0.0))
    }

    /// The pressure-release (soft) limit: `z -> 0` (represented as the
    /// smallest normal double; the admittance diverges).
    #[must_use]
    pub fn pressure_release() -> Self {
        Self(Complex::new(f64::MIN_POSITIVE, 0.0))
    }

    /// The impedance normalized by the medium's characteristic impedance
    /// `rho·c`: `z' = z/(rho·c)` (dimensionless).
    #[must_use]
    pub fn normalized(&self, medium: &Medium) -> Complex {
        let rho_c = medium.characteristic_impedance();
        Complex::new(self.0.re / rho_c, self.0.im / rho_c)
    }

    /// Normal-incidence pressure reflection factor `r = (z' − 1)/(z' + 1)`.
    /// Rigid walls reflect in phase (`r = 1`); pressure-release walls reflect
    /// with phase reversal (`r = −1`).
    #[must_use]
    pub fn reflection_factor_normal(&self, medium: &Medium) -> Complex {
        let z = self.normalized(medium);
        // The rigid limit z' -> inf makes (z'-1)/(z'+1) an inf/inf form;
        // the physical reflection factor there is exactly +1.
        if !z.re.is_finite() || !z.im.is_finite() {
            return Complex::new(1.0, 0.0);
        }
        // (z − 1)/(z + 1), componentwise to avoid relying on extra Complex ops.
        let num = Complex::new(z.re - 1.0, z.im);
        let den = Complex::new(z.re + 1.0, z.im);
        complex_div(num, den)
    }

    /// Normal-incidence absorption coefficient `alpha = 1 − |r|²` in [0, 1].
    #[must_use]
    pub fn absorption_normal(&self, medium: &Medium) -> f64 {
        let r = self.reflection_factor_normal(medium);
        1.0 - (r.re * r.re + r.im * r.im)
    }

    /// A real (locally reacting) impedance whose normal-incidence absorption
    /// equals `alpha` in [0, 1): `z' = (1 + sqrt(1−alpha))/(1 − sqrt(1−alpha))`
    /// normalized, i.e. the physically consistent inverse of
    /// [`SpecificImpedance::absorption_normal`].
    ///
    /// # Errors
    /// Returns [`AcousticError::NonPositive`] unless `0 <= alpha < 1`
    /// (`alpha = 1` would need `z = rho·c` exactly — use
    /// [`SpecificImpedance::real`] with the medium's characteristic
    /// impedance for that case).
    pub fn from_absorption(alpha: f64, medium: &Medium) -> Result<Self, AcousticError> {
        if !(0.0..1.0).contains(&alpha) {
            return Err(AcousticError::NonPositive("absorption must lie in [0, 1)"));
        }
        let s = libm::sqrt(1.0 - alpha);
        // alpha = 0 diverges to the rigid limit (infinite impedance), which
        // `reflection_factor_normal` handles explicitly.
        let z_norm = (1.0 + s) / (1.0 - s);
        Ok(Self(Complex::new(
            z_norm * medium.characteristic_impedance(),
            0.0,
        )))
    }

    /// The reciprocal admittance.
    #[must_use]
    pub fn admittance(&self) -> SpecificAdmittance {
        SpecificAdmittance(complex_div(Complex::new(1.0, 0.0), self.0))
    }
}

impl SpecificAdmittance {
    /// The reciprocal impedance.
    #[must_use]
    pub fn impedance(&self) -> SpecificImpedance {
        SpecificImpedance(complex_div(Complex::new(1.0, 0.0), self.0))
    }
}

/// Componentwise complex division `(a.re + i a.im)/(b.re + i b.im)`.
fn complex_div(a: Complex, b: Complex) -> Complex {
    let den = b.re * b.re + b.im * b.im;
    Complex::new(
        (a.re * b.re + a.im * b.im) / den,
        (a.im * b.re - a.re * b.im) / den,
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    const AIR: Medium = Medium {
        density: 1.2041,
        speed_of_sound: 343.21,
    };
    const TOL: f64 = 1e-12;

    #[test]
    fn matched_impedance_absorbs_completely() {
        let z = SpecificImpedance::real(AIR.characteristic_impedance()).unwrap();
        assert!(z.absorption_normal(&AIR) < 1e-9 || (z.absorption_normal(&AIR) - 1.0).abs() < 1e-9);
        assert!((z.absorption_normal(&AIR) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn rigid_wall_reflects_in_phase() {
        let r = SpecificImpedance::rigid().reflection_factor_normal(&AIR);
        assert!((r.re - 1.0).abs() < TOL);
        assert!(r.im.abs() < TOL);
    }

    #[test]
    fn pressure_release_wall_reflects_out_of_phase() {
        let r = SpecificImpedance::pressure_release().reflection_factor_normal(&AIR);
        assert!((r.re + 1.0).abs() < 1e-6, "r = {r:?}");
    }

    #[test]
    fn absorption_impedance_round_trip() {
        for &alpha in &[0.0, 0.05, 0.3, 0.7, 0.95, 0.999] {
            let z = SpecificImpedance::from_absorption(alpha, &AIR).unwrap();
            let alpha_back = z.absorption_normal(&AIR);
            assert!(
                (alpha_back - alpha).abs() < 1e-9,
                "alpha {alpha} -> {alpha_back}"
            );
        }
    }

    #[test]
    fn admittance_is_reciprocal_of_impedance() {
        let z = SpecificImpedance(Complex::new(300.0, 150.0));
        let y = z.admittance();
        let z_back = y.impedance();
        assert!((z_back.0.re - 300.0).abs() < TOL);
        assert!((z_back.0.im - 150.0).abs() < TOL);
    }

    #[test]
    fn rejects_zero_impedance() {
        assert!(SpecificImpedance::new(Complex::new(0.0, 0.0)).is_err());
        assert!(SpecificImpedance::real(0.0).is_err());
    }
}
