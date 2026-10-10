//! Harmonic wave primitives and the governing equations of linear acoustics.
//!
//! Time convention: `exp(+i·omega·t)`. An outgoing spherical wave therefore
//! carries the factor `exp(−i·k·r)` and the free-field Green's function of
//! the Helmholtz operator `(lap + k²)` is
//!
//! ```text
//! G(r) = exp(−i·k·r) / (4·pi·r)
//! ```

use libm::sqrt;
use tpt_math_linalg_complex::Complex;

/// Frequency in hertz.
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Frequency(f64);

impl Frequency {
    /// Creates a frequency in hertz.
    ///
    /// # Errors
    /// Returns [`crate::AcousticError::NonPositive`] for `value <= 0`.
    pub fn hertz(value: f64) -> Result<Self, crate::AcousticError> {
        if value <= 0.0 {
            return Err(crate::AcousticError::NonPositive("frequency"));
        }
        Ok(Self(value))
    }

    /// Creates a frequency in kilohertz.
    ///
    /// # Errors
    /// See [`Frequency::hertz`].
    pub fn kilohertz(value: f64) -> Result<Self, crate::AcousticError> {
        Self::hertz(value * 1000.0)
    }

    /// The value in hertz.
    #[must_use]
    pub const fn hz(&self) -> f64 {
        self.0
    }

    /// Angular frequency `omega = 2·pi·f` in rad/s.
    #[must_use]
    pub const fn angular(&self) -> f64 {
        core::f64::consts::TAU * self.0
    }

    /// Period `T = 1/f` in seconds.
    #[must_use]
    pub fn period(&self) -> f64 {
        1.0 / self.0
    }
}

/// Wavenumber in radians per metre (magnitude of the wave vector).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Wavenumber(f64);

impl Wavenumber {
    /// Creates a wavenumber in rad/m.
    ///
    /// # Errors
    /// Returns [`crate::AcousticError::NonPositive`] for `value < 0`.
    pub fn new(value: f64) -> Result<Self, crate::AcousticError> {
        if value < 0.0 {
            return Err(crate::AcousticError::NonPositive("wavenumber"));
        }
        Ok(Self(value))
    }

    /// The value in rad/m.
    #[must_use]
    pub const fn rad_per_m(&self) -> f64 {
        self.0
    }

    /// Wavenumber from frequency and speed: `k = 2·pi·f / c`.
    ///
    /// # Errors
    /// Returns [`crate::AcousticError::NonPositive`] for `f <= 0` or `c <= 0`.
    pub fn from_frequency(f: Frequency, c: f64) -> Result<Self, crate::AcousticError> {
        if c <= 0.0 {
            return Err(crate::AcousticError::NonPositive("speed of sound"));
        }
        Self::new(f.angular() / c)
    }
}

/// Angular frequency `omega = 2·pi·f` in rad/s.
#[must_use]
pub const fn angular_frequency(f: Frequency) -> f64 {
    f.angular()
}

/// Complex wavenumber of a lightly attenuated harmonic wave:
/// `k = omega/c − i·alpha` (amplitude decay `exp(−alpha·x)` under the
/// `exp(+i·omega·t)` convention).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ComplexWavenumber(pub Complex);

impl ComplexWavenumber {
    /// Builds a complex wavenumber from real part (rad/m) and (negative of)
    /// amplitude attenuation (Np/m).
    #[must_use]
    pub fn new(real_rad_per_m: f64, minus_imag_np_per_m: f64) -> Self {
        Self(Complex::new(real_rad_per_m, -minus_imag_np_per_m))
    }

    /// Propagates a complex amplitude over distance `x`:
    /// `a(x) = a·exp(−i·k·x)`. With `k = k_r − i·alpha` this is
    /// `a·exp(−alpha·x)·exp(−i·k_r·x)`: decay times phase lag.
    #[must_use]
    pub fn propagate(&self, amplitude: Complex, x: f64) -> Complex {
        // -i k x with k = k_r + i*k_im expands to (k_im - i*k_r) x, so the
        // magnitude factor is exp(k_im * x) = exp(-alpha * x).
        let decay = libm::exp(self.0.im * x);
        amplitude
            * Complex::new(
                decay * libm::cos(-self.0.re * x),
                decay * libm::sin(-self.0.re * x),
            )
    }
}

/// Plane-wave pressure amplitude at axial distance `x`:
/// `p(x) = p0·exp(−i·k·x)`.
#[must_use]
pub fn plane_wave_pressure(p0: f64, k: &ComplexWavenumber, x: f64) -> Complex {
    k.propagate(Complex::new(p0, 0.0), x)
}

/// Free-field (free-space) Green's function of the Helmholtz operator,
/// `G(r) = exp(−i·k·r)/(4·pi·r)` for `r > 0`.
///
/// # Errors
/// Returns [`crate::AcousticError::NonPositive`] for `r <= 0`.
pub fn free_field_green(k: f64, r: f64) -> Result<Complex, crate::AcousticError> {
    if r <= 0.0 {
        return Err(crate::AcousticError::NonPositive("Green's function radius"));
    }
    let magnitude = 1.0 / (4.0 * core::f64::consts::PI * r);
    Ok(Complex::from_polar(magnitude, -k * r))
}

/// Residual of the Helmholtz equation `(lap + k²)p = 0` for a sampled field.
#[must_use]
pub fn helmholtz_residual(laplacian_p: Complex, p: Complex, k: f64) -> Complex {
    laplacian_p + p * (k * k)
}

/// Residual of the linearized Euler momentum equation
/// `rho·(du/dt) + grad p = 0` in one Cartesian direction.
#[must_use]
pub fn euler_momentum_residual(du_dt: f64, dp_dx: f64, rho: f64) -> f64 {
    rho * du_dt + dp_dx
}

/// Residual of the acoustic wave equation `(1/c²)·(d²p/dt²) − lap p = 0`.
#[must_use]
pub fn wave_equation_residual(d2p_dt2: f64, laplacian_p: f64, c: f64) -> f64 {
    d2p_dt2 / (c * c) - laplacian_p
}

/// Speed of sound in an ideal gas: `c = sqrt(gamma·R_specific·T)`.
#[must_use]
pub fn speed_of_sound_ideal_gas(gamma: f64, r_specific: f64, temperature_k: f64) -> f64 {
    sqrt(gamma * r_specific * temperature_k)
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-12;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn frequency_conversions() {
        let f = Frequency::hertz(1000.0).unwrap();
        assert!(approx(f.angular(), core::f64::consts::TAU * 1000.0, TOL));
        assert!(approx(f.period(), 1.0e-3, TOL));
        let f2 = Frequency::kilohertz(2.0).unwrap();
        assert!(approx(f2.hz(), 2000.0, TOL));
        assert!(Frequency::hertz(0.0).is_err());
    }

    #[test]
    fn wavenumber_from_frequency() {
        let k = Wavenumber::from_frequency(Frequency::hertz(343.0).unwrap(), 343.0).unwrap();
        assert!(approx(k.rad_per_m(), core::f64::consts::TAU, 1e-9));
    }

    #[test]
    fn green_function_magnitude_and_phase() {
        let g = free_field_green(2.0, 3.0).unwrap();
        assert!(approx(
            g.norm(),
            1.0 / (4.0 * core::f64::consts::PI * 3.0),
            1e-15
        ));
        // Phase should be -k*r = -6 rad (mod 2*pi).
        let phase = g.arg();
        assert!(
            approx(phase, -6.0_f64.rem_euclid(core::f64::consts::TAU), 1e-9)
                || approx(phase, -6.0 + core::f64::consts::TAU, 1e-9)
        );
        assert!(free_field_green(1.0, 0.0).is_err());
    }

    #[test]
    fn plane_wave_with_attenuation_decays() {
        let k = ComplexWavenumber::new(1.0, 0.5);
        let p0 = plane_wave_pressure(1.0, &k, 0.0);
        assert!(approx(p0.re, 1.0, TOL));
        let p1 = plane_wave_pressure(1.0, &k, 2.0);
        // Amplitude decays by exp(-0.5*2) = e^-1.
        assert!(approx(p1.norm(), libm::exp(-1.0), 1e-12));
    }

    #[test]
    fn residuals_vanish_for_analytic_fields() {
        // Helmholtz: p = exp(-i k x) satisfies lap p + k^2 p = 0.
        let k = 3.0;
        let x = 0.7;
        let p = plane_wave_pressure(1.0, &ComplexWavenumber::new(k, 0.0), x);
        // lap p = -k^2 p exactly for the plane wave.
        let lap = p * (-(k * k));
        assert!(helmholtz_residual(lap, p, k).norm() < 1e-9);
        // Euler: a plane wave has rho·du/dt = -dp/dx exactly, so the
        // residual vanishes term by term.
        let rho = 1.2;
        assert!(euler_momentum_residual(1.0, -rho * 1.0, rho).abs() < TOL);
        // Wave equation: p = sin(t) in nondimensional units with c = 1:
        // d2p/dt2 = -sin t, lap p = -sin t -> residual 0.
        assert!(wave_equation_residual(-1.0, -1.0, 1.0).abs() < TOL);
    }

    #[test]
    fn ideal_gas_speed_of_sound() {
        // Air: gamma = 1.4, R = 287.058, T = 293.15 K -> ~343.0 m/s.
        let c = speed_of_sound_ideal_gas(1.4, 287.058, 293.15);
        assert!((c - 343.0).abs() < 0.5, "c = {c}");
    }
}
