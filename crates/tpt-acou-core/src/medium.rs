//! Acoustic media: fluids and solids with their derived acoustic constants.
//!
//! This module is the native replacement for the (not-yet-existing) external
//! `tpt-materials` dependency: it ships a small, tested table of standard
//! media plus the derived quantities the rest of the library needs
//! (characteristic impedance, bulk modulus, wavenumbers, plate stiffness).

use libm::sqrt;

use crate::error::AcousticError;
use crate::wave::{angular_frequency, Frequency, Wavenumber};
use tpt_math_linalg_complex::Complex;

/// A homogeneous, quiescent acoustic fluid (lossless except where an
/// attenuation coefficient is supplied explicitly).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Medium {
    /// Density `rho` in kg/m³.
    pub density: f64,
    /// Small-signal speed of sound `c` in m/s.
    pub speed_of_sound: f64,
}

impl Medium {
    /// Creates a fluid from density (kg/m³) and speed of sound (m/s).
    ///
    /// # Errors
    /// Returns [`AcousticError::NonPositive`] if `density` or
    /// `speed_of_sound` is not strictly positive.
    pub fn new(density: f64, speed_of_sound: f64) -> Result<Self, AcousticError> {
        if density <= 0.0 {
            return Err(AcousticError::NonPositive("density"));
        }
        if speed_of_sound <= 0.0 {
            return Err(AcousticError::NonPositive("speed of sound"));
        }
        Ok(Self {
            density,
            speed_of_sound,
        })
    }

    /// Standard dry air at 20 °C and 101.325 kPa
    /// (`rho ≈ 1.204 kg/m³`, `c ≈ 343.2 m/s`).
    #[must_use]
    pub fn air_20c() -> Self {
        Self::air(20.0)
    }

    /// Dry air at `temperature_c` °C and 101.325 kPa, using the standard
    /// approximations `c = 331.3·sqrt(1 + T/273.15)` m/s and the ideal-gas
    /// density `rho = p/(R_specific·T)` with `R_specific = 287.058 J/(kg·K)`.
    #[must_use]
    pub fn air(temperature_c: f64) -> Self {
        let t_k = temperature_c + 273.15;
        let c = 331.3 * sqrt(1.0 + temperature_c / 273.15);
        let density = 101_325.0 / (287.058 * t_k);
        Self {
            density,
            speed_of_sound: c,
        }
    }

    /// Fresh water at 20 °C (`rho ≈ 998 kg/m³`, `c ≈ 1481 m/s`).
    #[must_use]
    pub fn water_20c() -> Self {
        Self {
            density: 998.0,
            speed_of_sound: 1481.0,
        }
    }

    /// Standard seawater at 20 °C, atmospheric pressure
    /// (`rho ≈ 1025 kg/m³`, `c ≈ 1500 m/s`).
    #[must_use]
    pub fn seawater_20c() -> Self {
        Self {
            density: 1025.0,
            speed_of_sound: 1500.0,
        }
    }

    /// Adiabatic bulk modulus `K = rho·c²` in Pa.
    #[must_use]
    pub fn bulk_modulus(&self) -> f64 {
        self.density * self.speed_of_sound * self.speed_of_sound
    }

    /// Characteristic (specific) acoustic impedance `rho·c` in Pa·s/m.
    #[must_use]
    pub fn characteristic_impedance(&self) -> f64 {
        self.density * self.speed_of_sound
    }

    /// Real wavenumber `k = omega/c = 2πf/c` in rad/m.
    ///
    /// # Errors
    /// Returns [`AcousticError::NonPositive`] for `f <= 0`.
    pub fn wavenumber(&self, f: Frequency) -> Result<Wavenumber, AcousticError> {
        if f.hz() <= 0.0 {
            return Err(AcousticError::NonPositive("frequency"));
        }
        let k = Wavenumber::new(angular_frequency(f) / self.speed_of_sound)?;
        Ok(k)
    }

    /// Complex wavenumber `k = omega/c − i·alpha` for a medium with viscous /
    /// relaxation attenuation `alpha` in nepers per metre (amplitude decay
    /// `exp(−alpha·x)`).
    ///
    /// # Errors
    /// Returns [`AcousticError::NonPositive`] for `f <= 0` or `alpha < 0`.
    pub fn complex_wavenumber(
        &self,
        f: Frequency,
        alpha_np_per_m: f64,
    ) -> Result<Complex, AcousticError> {
        if alpha_np_per_m < 0.0 {
            return Err(AcousticError::NonPositive("attenuation"));
        }
        let k_real = self.wavenumber(f)?;
        Ok(Complex::new(k_real.rad_per_m(), -alpha_np_per_m))
    }
}

/// A homogeneous, isotropic linear-elastic solid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Solid {
    /// Density `rho` in kg/m³.
    pub density: f64,
    /// Young's modulus `E` in Pa.
    pub youngs_modulus: f64,
    /// Poisson's ratio `nu` (dimensionless, `−1 < nu < 0.5`).
    pub poisson_ratio: f64,
}

impl Solid {
    /// Creates a solid from density (kg/m³), Young's modulus (Pa) and
    /// Poisson's ratio.
    ///
    /// # Errors
    /// Returns [`AcousticError::NonPositive`] for non-positive density or
    /// modulus, and [`AcousticError::InvalidArgument`] for Poisson's ratio
    /// outside `−1 < nu < 0.5`.
    pub fn new(
        density: f64,
        youngs_modulus: f64,
        poisson_ratio: f64,
    ) -> Result<Self, AcousticError> {
        if density <= 0.0 {
            return Err(AcousticError::NonPositive("density"));
        }
        if youngs_modulus <= 0.0 {
            return Err(AcousticError::NonPositive("Young's modulus"));
        }
        if !(poisson_ratio > -1.0 && poisson_ratio < 0.5) {
            return Err(AcousticError::InvalidArgument(
                "Poisson's ratio must lie in (-1, 0.5)",
            ));
        }
        Ok(Self {
            density,
            youngs_modulus,
            poisson_ratio,
        })
    }

    /// Bending (flexural) rigidity of a thin plate of thickness `h`:
    /// `D = E·h³ / (12·(1 − nu²))` in N·m.
    ///
    /// # Errors
    /// Returns [`AcousticError::NonPositive`] for `h <= 0`.
    pub fn plate_bending_stiffness(&self, h: f64) -> Result<f64, AcousticError> {
        if h <= 0.0 {
            return Err(AcousticError::NonPositive("plate thickness"));
        }
        Ok(self.youngs_modulus * h * h * h
            / (12.0 * (1.0 - self.poisson_ratio * self.poisson_ratio)))
    }

    /// Thin-plate flexural (bending) wave speed at angular frequency
    /// `omega`: `c_p = (D·omega² / (rho·h))^{1/4}` in m/s (dispersive).
    ///
    /// # Errors
    /// Returns [`AcousticError::NonPositive`] for `h <= 0` or `omega <= 0`.
    pub fn flexural_wave_speed(&self, h: f64, omega: f64) -> Result<f64, AcousticError> {
        if h <= 0.0 {
            return Err(AcousticError::NonPositive("plate thickness"));
        }
        if omega <= 0.0 {
            return Err(AcousticError::NonPositive("angular frequency"));
        }
        let d = self.plate_bending_stiffness(h)?;
        Ok(libm::pow(d * omega * omega / (self.density * h), 0.25))
    }

    /// Critical coincidence frequency of a plate of thickness `h` radiating
    /// into the fluid `fluid`: `f_c = c²/(2π)·sqrt(rho_s/(D))` with the mass
    /// surface density `rho_s = rho·h` in Hz (above which plate radiation is
    /// efficient).
    ///
    /// # Errors
    /// Returns [`AcousticError::NonPositive`] for `h <= 0`.
    pub fn coincidence_frequency(&self, h: f64, fluid: &Medium) -> Result<f64, AcousticError> {
        if h <= 0.0 {
            return Err(AcousticError::NonPositive("plate thickness"));
        }
        let d = self.plate_bending_stiffness(h)?;
        let rho_s = self.density * h;
        Ok(fluid.speed_of_sound * fluid.speed_of_sound / core::f64::consts::TAU * sqrt(rho_s / d))
    }
}

#[cfg(feature = "materials")]
impl Solid {
    /// Builds a solid from a `tpt-mat-database` record (the
    /// github.com/tpt-solutions/tpt-materials integration named in the
    /// spec): Young's modulus converts GPa -> Pa; density and Poisson's
    /// ratio carry over.
    ///
    /// # Errors
    /// Propagates [`Solid::new`] validation errors.
    pub fn from_material_record(
        record: &tpt_mat_database::MaterialRecord,
    ) -> Result<Self, AcousticError> {
        Self::new(
            record.mechanical.density_kg_per_m3,
            record.mechanical.youngs_modulus_gpa * 1.0e9,
            record.mechanical.poissons_ratio,
        )
    }

    /// Looks up a builtin database material by (substring, case-sensitive)
    /// name — e.g. `"Al 6061-T6"` or `"AISI 304"` — and converts it.
    ///
    /// # Errors
    /// Returns [`AcousticError::InvalidArgument`] when no builtin record
    /// matches, and propagates [`Solid::from_material_record`] errors.
    pub fn from_builtin_material(name: &str) -> Result<Self, AcousticError> {
        let db = tpt_mat_database::MaterialsDatabase::load_builtin();
        let record =
            db.search_by_name(name)
                .into_iter()
                .next()
                .ok_or(AcousticError::InvalidArgument(
                    "no builtin material record matches",
                ))?;
        Self::from_material_record(record)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::wave::Frequency;

    const TOL: f64 = 1e-9;

    #[test]
    fn air_20c_standard_values() {
        let air = Medium::air_20c();
        assert!(
            (air.speed_of_sound - 343.21).abs() < 0.05,
            "c = {}",
            air.speed_of_sound
        );
        assert!(
            (air.density - 1.2041).abs() < 0.001,
            "rho = {}",
            air.density
        );
        assert!((air.characteristic_impedance() - 413.3).abs() < 0.5);
        // Adiabatic bulk modulus of air ~ 1.42e5 Pa.
        assert!((air.bulk_modulus() - 1.42e5).abs() < 1.0e3);
    }

    #[test]
    fn water_reference_values() {
        let water = Medium::water_20c();
        assert!((water.characteristic_impedance() - 998.0 * 1481.0).abs() < 1e-6);
        assert!((water.bulk_modulus() - 998.0 * 1481.0 * 1481.0).abs() < 1.0);
    }

    #[test]
    fn wavenumber_matches_two_pi_over_lambda() {
        let air = Medium::air_20c();
        let k = air.wavenumber(Frequency::hertz(1000.0).unwrap()).unwrap();
        assert!((k.rad_per_m() - core::f64::consts::TAU / 0.343_21).abs() < 1e-3);
    }

    #[test]
    fn complex_wavenumber_attenuates() {
        let air = Medium::air_20c();
        let k = air
            .complex_wavenumber(Frequency::hertz(1000.0).unwrap(), 0.5)
            .unwrap();
        assert!(
            (k.re
                - air
                    .wavenumber(Frequency::hertz(1000.0).unwrap())
                    .unwrap()
                    .rad_per_m())
            .abs()
                < TOL
        );
        assert!((k.im + 0.5).abs() < TOL);
    }

    #[test]
    fn rejects_nonpositive_medium() {
        assert!(Medium::new(0.0, 343.0).is_err());
        assert!(Medium::new(1.2, -1.0).is_err());
        assert!(Frequency::hertz(0.0).is_err());
        assert!(Frequency::hertz(-100.0).is_err());
    }

    #[test]
    fn plate_stiffness_of_steel() {
        let steel = Solid::new(7850.0, 2.0e11, 0.3).unwrap();
        let d = steel.plate_bending_stiffness(0.01).unwrap();
        // D = 2e11 * 1e-6 / (12 * 0.91)
        assert!((d - 2.0e11 * 1.0e-6 / (12.0 * 0.91)).abs() < 1.0);
    }

    #[test]
    fn coincidence_frequency_order_of_magnitude() {
        // 1 mm steel plate in air has f_c ≈ 12 kHz (Cremer/Heckl reference).
        let steel = Solid::new(7850.0, 2.0e11, 0.3).unwrap();
        let fc = steel
            .coincidence_frequency(0.001, &Medium::air_20c())
            .unwrap();
        assert!((fc - 12_000.0).abs() < 1_500.0, "fc = {fc}");
    }

    #[cfg(feature = "materials")]
    #[test]
    fn builtin_aluminium_record_converts() {
        let al = Solid::from_builtin_material("Al 6061-T6").unwrap();
        assert!((al.density - 2700.0).abs() < 1.0, "rho = {}", al.density);
        assert!(
            (al.youngs_modulus - 68.9e9).abs() < 1.0e9,
            "E = {}",
            al.youngs_modulus
        );
        assert!((al.poisson_ratio - 0.33).abs() < 0.01);
    }

    #[cfg(feature = "materials")]
    #[test]
    fn unknown_material_is_rejected() {
        assert!(Solid::from_builtin_material("unobtanium").is_err());
    }

    #[test]
    fn rejects_bad_poisson() {
        assert!(Solid::new(7850.0, 2.0e11, 0.5).is_err());
        assert!(Solid::new(7850.0, 2.0e11, -1.5).is_err());
    }
}
