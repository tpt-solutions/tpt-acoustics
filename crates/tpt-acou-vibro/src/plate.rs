//! Analytic simply-supported thin-plate modal model (Kirchhoff–Love).
//!
//! Mode `(m, n)` of a simply supported rectangular plate:
//!
//! ```text
//! phi_mn(x, y) = sin(m*pi*x/Lx) * sin(n*pi*y/Ly)
//! omega_mn     = pi^2 * sqrt(D / (rho*h)) * ((m/Lx)^2 + (n/Ly)^2)
//! modal mass   = rho*h*Lx*Ly/4
//! ```
//!
//! with the bending stiffness `D = E h^3 / (12 (1 - nu^2))` from
//! `tpt_acou_core::medium::Solid`.

use tpt_acou_core::error::AcousticError;
use tpt_acou_core::medium::Solid;

/// A simply supported rectangular thin plate.
#[derive(Debug, Clone)]
pub struct SimplySupportedPlate {
    /// Plate extent in x (m).
    pub lx: f64,
    /// Plate extent in y (m).
    pub ly: f64,
    /// Plate thickness (m).
    pub thickness: f64,
    /// Structural material.
    pub material: Solid,
}

impl SimplySupportedPlate {
    /// Creates a plate.
    ///
    /// # Errors
    /// Returns [`AcousticError::NonPositive`] for non-positive dimensions.
    pub fn new(lx: f64, ly: f64, thickness: f64, material: Solid) -> Result<Self, AcousticError> {
        if lx <= 0.0 || ly <= 0.0 {
            return Err(AcousticError::NonPositive("plate extent"));
        }
        if thickness <= 0.0 {
            return Err(AcousticError::NonPositive("plate thickness"));
        }
        Ok(Self {
            lx,
            ly,
            thickness,
            material,
        })
    }

    /// Bending stiffness `D = E h^3 / (12 (1 - nu^2))` in N·m.
    #[must_use]
    pub fn bending_stiffness(&self) -> f64 {
        self.material
            .plate_bending_stiffness(self.thickness)
            .expect("thickness validated positive")
    }

    /// Mass surface density `rho*h` in kg/m².
    #[must_use]
    pub fn mass_per_area(&self) -> f64 {
        self.material.density * self.thickness
    }

    /// Natural (angular) frequency `omega_mn` in rad/s for mode `(m, n)`
    /// (both ≥ 1 for a simply supported plate).
    ///
    /// # Errors
    /// Returns [`AcousticError::NonPositive`] for mode indices < 1.
    pub fn mode_omega(&self, m: u32, n: u32) -> Result<f64, AcousticError> {
        if m == 0 || n == 0 {
            return Err(AcousticError::NonPositive("simply supported mode index"));
        }
        let d = self.bending_stiffness();
        let rho_h = self.mass_per_area();
        let mf = f64::from(m) / self.lx;
        let nf = f64::from(n) / self.ly;
        let k2 = core::f64::consts::PI * core::f64::consts::PI * (mf * mf + nf * nf);
        Ok(k2 * libm::sqrt(d / rho_h))
    }

    /// Modal mass of mode `(m, n)`: `rho*h*Lx*Ly/4` (independent of `m, n`
    /// for the sine basis).
    #[must_use]
    pub fn modal_mass(&self) -> f64 {
        self.mass_per_area() * self.lx * self.ly / 4.0
    }

    /// Mode shape value at `(x, y)` (unnormalized sine product).
    #[must_use]
    pub fn mode_shape(&self, m: u32, n: u32, x: f64, y: f64) -> f64 {
        libm::sin(core::f64::consts::PI * f64::from(m) * x / self.lx)
            * libm::sin(core::f64::consts::PI * f64::from(n) * y / self.ly)
    }

    /// Harmonic point mobility (velocity per unit point force) at `(x, y)`
    /// for excitation at `(x0, y0)` with damping ratio `zeta`, via modal
    /// superposition over the given modes.
    ///
    /// Returns the complex ratio `v/F` in m/(s·N) using the `exp(+iωt)`
    /// convention (driving-point and transfer mobilities both come from the
    /// same modal sum).
    ///
    /// # Errors
    /// Propagates [`SimplySupportedPlate::mode_omega`] errors.
    pub fn point_mobility(
        &self,
        modes: &[(u32, u32)],
        omega: f64,
        zeta: f64,
        (x, y): (f64, f64),
        (x0, y0): (f64, f64),
    ) -> Result<Complex, AcousticError> {
        let mut acc = Complex::new(0.0, 0.0);
        let i_omega = Complex::new(0.0, omega);
        for &(m, n) in modes {
            let wn = self.mode_omega(m, n)?;
            let phi = self.mode_shape(m, n, x, y);
            let phi0 = self.mode_shape(m, n, x0, y0);
            let mn = self.modal_mass();
            // q_n = F*phi0 / (mn*(wn^2 - w^2 + 2i zeta wn w)); v = i w q_n.
            let denom = Complex::new(mn * (wn * wn - omega * omega), mn * 2.0 * zeta * wn * omega);
            acc = acc + i_omega * (phi * phi0) / denom;
        }
        Ok(acc)
    }
}

use tpt_math_linalg_complex::Complex;
