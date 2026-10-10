//! Phantom-typed acoustic quantities and sound pressure levels.
//!
//! A [`Quantity`] carries its physical unit in a zero-sized type parameter, so
//! adding a [`Pressure`] to a [`ParticleVelocity`] is a *compile* error, while
//! multiplying them yields an [`Intensity`] — the product the physics allows.
//!
//! ```
//! use tpt_acou_core::quantity::{Mps, Pa, ParticleVelocity, Pressure, Wm2};
//!
//! let p = Pressure::<f64, Pa>::new(1.0);              // 1 Pa
//! let u = ParticleVelocity::<f64, Mps>::new(1.0 / 413.0); // ~2.4 mm/s in air
//! let intensity: tpt_acou_core::Intensity<f64, Wm2> = p * u;
//! assert!((intensity.value() - 1.0 / 413.0).abs() < 1e-12);
//! ```
//!
//! The second type parameter of the aliases is the *reference* marker (see
//! [`Pa`]); it exists so that quantities measured against different zero
//! points cannot be silently mixed.

use core::marker::PhantomData;
use core::ops::{Add, Div, Mul, Neg, Sub};

use crate::error::AcousticError;
use crate::wave::Frequency;

/// Marker unit: SI pascal (Pa) for acoustic pressure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Pa;

/// Marker unit: metres per second (m/s) for particle velocity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Mps;

/// Marker unit: watt (W) for sound power.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Watt;

/// Marker unit: watts per square metre (W/m²) for sound intensity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Wm2;

/// Marker unit: pascal-second per metre (Pa·s/m, "rayl") for specific acoustic
/// impedance.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Rayl;

/// A physical scalar tagged with a phantom unit type.
///
/// The stored value is in the base SI unit named by `U` (Pa, m/s, W, W/m²,
/// Pa·s/m). Arithmetic preserves the unit; incompatible units never meet.
#[derive(Debug, Clone, Copy)]
pub struct Quantity<T = f64, U = Pa>(T, PhantomData<fn() -> U>);

impl<T, U> Quantity<T, U> {
    /// Wraps a raw value (in the base SI unit of `U`).
    #[must_use]
    pub const fn new(value: T) -> Self {
        Self(value, PhantomData)
    }

    /// The raw value in the base SI unit of `U`.
    #[must_use]
    pub const fn value(&self) -> &T {
        &self.0
    }

    /// Consumes the quantity, returning the raw value.
    #[must_use]
    pub fn into_value(self) -> T {
        self.0
    }
}

impl<T: Copy, U> Quantity<T, U> {
    /// The raw value by copy.
    #[must_use]
    pub fn get(&self) -> T {
        self.0
    }
}

impl<T: Add<Output = T>, U> Add for Quantity<T, U> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        Self(self.0 + rhs.0, PhantomData)
    }
}

impl<T: Sub<Output = T>, U> Sub for Quantity<T, U> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        Self(self.0 - rhs.0, PhantomData)
    }
}

impl<T: Neg<Output = T>, U> Neg for Quantity<T, U> {
    type Output = Self;
    fn neg(self) -> Self {
        Self(-self.0, PhantomData)
    }
}

impl<T: Mul<Output = T>, U> Mul<T> for Quantity<T, U> {
    type Output = Self;
    fn mul(self, rhs: T) -> Self {
        Self(self.0 * rhs, PhantomData)
    }
}

impl<T: Div<Output = T>, U> Div<T> for Quantity<T, U> {
    type Output = Self;
    fn div(self, rhs: T) -> Self {
        Self(self.0 / rhs, PhantomData)
    }
}

/// Physically valid unit products. `Quantity` pairs only implement `Mul` /
/// `Div` where the physics defines the result, so dimensionally meaningless
/// products never type-check.
impl<T: Mul<Output = T>> Mul<Quantity<T, Mps>> for Quantity<T, Pa> {
    type Output = Quantity<T, Wm2>;
    fn mul(self, rhs: Quantity<T, Mps>) -> Quantity<T, Wm2> {
        Quantity::new(self.0 * rhs.0)
    }
}

impl<T: Div<Output = T>> Div<Quantity<T, Mps>> for Quantity<T, Pa> {
    type Output = Quantity<T, Rayl>;
    fn div(self, rhs: Quantity<T, Mps>) -> Quantity<T, Rayl> {
        Quantity::new(self.0 / rhs.0)
    }
}

/// Acoustic pressure. Base unit: pascal.
pub type Pressure<T = f64, U = Pa> = Quantity<T, U>;

/// Particle (medium) velocity. Base unit: m/s.
pub type ParticleVelocity<T = f64, U = Mps> = Quantity<T, U>;

/// Sound power. Base unit: watt.
pub type Power<T = f64, U = Watt> = Quantity<T, U>;

/// Sound intensity (energy flux). Base unit: W/m².
pub type Intensity<T = f64, U = Wm2> = Quantity<T, U>;

/// Specific acoustic impedance p/u. Base unit: Pa·s/m.
pub type SpecificImpedanceScalar<T = f64, U = Rayl> = Quantity<T, U>;

/// Speed of sound (kept in the [`Mps`] dimension, distinct alias for clarity).
pub type Speed<T = f64> = Quantity<T, Mps>;

/// Standard reference sound pressure: 20 µPa (ANSI S1.1 / IEC 61672), the
/// approximate threshold of human hearing at 1 kHz in air.
pub const REFERENCE_PRESSURE_PA: f64 = 20.0e-6;

/// Standard reference sound intensity: 1 pW/m², the free-progressive-wave
/// intensity corresponding to [`REFERENCE_PRESSURE_PA`] in air at 20 °C.
pub const REFERENCE_INTENSITY_WM2: f64 = 1.0e-12;

/// Standard reference sound power: 1 pW (re 1 pW), used with the decibel
/// sound power level Lw.
pub const REFERENCE_POWER_W: f64 = 1.0e-12;

/// Sound pressure level in decibels, always relative to
/// [`REFERENCE_PRESSURE_PA`] (20 µPa) for sound-in-air use.
///
/// SPL is deliberately *not* a [`Quantity`]: the decibel is a logarithmic
/// unit, so linear arithmetic (`Add`, `Mul`) would be wrong. Energetic
/// combination goes through [`SoundPressureLevel::combine`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SoundPressureLevel(f64);

impl SoundPressureLevel {
    /// Computes SPL from an RMS pressure: `Lp = 20·log10(p_rms / 20 µPa)` dB.
    ///
    /// # Panics
    /// Panics if `p_rms` is not strictly positive (a log-scale level of a
    /// non-positive pressure is undefined).
    #[must_use]
    #[track_caller]
    pub fn from_pressure_rms(p_rms: Pressure) -> Self {
        assert!(p_rms.get() > 0.0, "SPL requires a positive RMS pressure");
        Self(20.0 * libm::log10(p_rms.get() / REFERENCE_PRESSURE_PA))
    }

    /// Fallible [`from_pressure_rms`](Self::from_pressure_rms).
    ///
    /// # Errors
    /// Returns [`AcousticError::InvalidArgument`] unless `p_rms` is positive
    /// and finite.
    pub fn try_from_pressure_rms(p_rms: Pressure) -> Result<Self, AcousticError> {
        let p = p_rms.get();
        if !(p > 0.0 && p.is_finite()) {
            return Err(AcousticError::InvalidArgument(
                "SPL requires a positive finite RMS pressure",
            ));
        }
        Ok(Self(20.0 * libm::log10(p / REFERENCE_PRESSURE_PA)))
    }

    /// The RMS pressure corresponding to this level.
    #[must_use]
    pub fn to_pressure_rms(&self) -> Pressure {
        Pressure::new(REFERENCE_PRESSURE_PA * libm::pow(10.0, self.0 / 20.0))
    }

    /// Computes a sound intensity level (re [`REFERENCE_INTENSITY_WM2`]) from
    /// an intensity: `Li = 10·log10(I / 1 pW/m²)` dB.
    ///
    /// # Panics
    /// Panics if `intensity` is not strictly positive.
    #[must_use]
    #[track_caller]
    pub fn from_intensity(intensity: Intensity) -> Self {
        assert!(intensity.get() > 0.0, "intensity level requires I > 0");
        Self(10.0 * libm::log10(intensity.get() / REFERENCE_INTENSITY_WM2))
    }

    /// Fallible [`from_intensity`](Self::from_intensity).
    ///
    /// # Errors
    /// Returns [`AcousticError::InvalidArgument`] unless `intensity` is
    /// positive and finite.
    pub fn try_from_intensity(intensity: Intensity) -> Result<Self, AcousticError> {
        let i = intensity.get();
        if !(i > 0.0 && i.is_finite()) {
            return Err(AcousticError::InvalidArgument(
                "intensity level requires a positive finite I",
            ));
        }
        Ok(Self(10.0 * libm::log10(i / REFERENCE_INTENSITY_WM2)))
    }

    /// The intensity corresponding to this level (interpreted as an intensity
    /// level).
    #[must_use]
    pub fn to_intensity(&self) -> Intensity {
        Intensity::new(REFERENCE_INTENSITY_WM2 * libm::pow(10.0, self.0 / 10.0))
    }

    /// Computes a sound power level `Lw = 10·log10(W / 1 pW)` dB.
    ///
    /// # Panics
    /// Panics if `power` is not strictly positive.
    #[must_use]
    #[track_caller]
    pub fn from_power(power: Power) -> Self {
        assert!(power.get() > 0.0, "power level requires W > 0");
        Self(10.0 * libm::log10(power.get() / REFERENCE_POWER_W))
    }

    /// Fallible [`from_power`](Self::from_power).
    ///
    /// # Errors
    /// Returns [`AcousticError::InvalidArgument`] unless `power` is positive
    /// and finite.
    pub fn try_from_power(power: Power) -> Result<Self, AcousticError> {
        let w = power.get();
        if !(w > 0.0 && w.is_finite()) {
            return Err(AcousticError::InvalidArgument(
                "power level requires a positive finite W",
            ));
        }
        Ok(Self(10.0 * libm::log10(w / REFERENCE_POWER_W)))
    }

    /// The sound power corresponding to this level (interpreted as a power
    /// level).
    #[must_use]
    pub fn to_power(&self) -> Power {
        Power::new(REFERENCE_POWER_W * libm::pow(10.0, self.0 / 10.0))
    }

    /// Energetically combines independent levels:
    /// `L = 10·log10(Σ 10^(Li/10))`.
    ///
    /// Returns `None` if `levels` is empty.
    #[must_use]
    pub fn combine(levels: &[SoundPressureLevel]) -> Option<Self> {
        let mut total = 0.0;
        for l in levels {
            total += libm::pow(10.0, l.0 / 10.0);
        }
        if levels.is_empty() {
            None
        } else {
            Some(Self(10.0 * libm::log10(total)))
        }
    }

    /// Energetic mean of levels: `L = 10·log10((1/N)·Σ 10^(Li/10))`.
    ///
    /// Returns `None` if `levels` is empty.
    #[must_use]
    pub fn average(levels: &[SoundPressureLevel]) -> Option<Self> {
        let sum = Self::combine(levels)?;
        Some(sum.add_db(-10.0 * libm::log10(levels.len() as f64)))
    }

    /// Adds an energetic level difference (e.g. `+3.0103` dB for doubling).
    #[must_use]
    pub fn add_db(&self, db: f64) -> Self {
        Self(self.0 + db)
    }

    /// The raw value in decibels.
    #[must_use]
    pub const fn db(&self) -> f64 {
        self.0
    }

    /// A-weighted level approximation hook: applies a level delta at a given
    /// frequency, in decibels. The authoritative A-weighting curve lives in
    /// `tpt-acou-psycho`; this helper lets level bookkeeping carry a frequency
    /// dependent correction without importing that crate.
    #[must_use]
    pub fn weighted(&self, delta_db: f64) -> Self {
        self.add_db(delta_db)
    }
}

/// Plane-progressive-wave intensity from an RMS pressure in a medium of
/// characteristic impedance `rho_c`: `I = p² / (ρ c)`.
#[must_use]
pub fn intensity_from_pressure(p_rms: Pressure, rho_c: f64) -> Intensity {
    let p = p_rms.get();
    Intensity::new(p * p / rho_c)
}

/// RMS pressure of a plane-progressive wave of given intensity in a medium of
/// characteristic impedance `rho_c`: `p = sqrt(I·ρc)`.
#[must_use]
pub fn pressure_from_intensity(intensity: Intensity, rho_c: f64) -> Pressure {
    Pressure::new(libm::sqrt(intensity.get() * rho_c))
}

/// Particle velocity amplitude of a plane wave: `u = p / (ρ c)`.
#[must_use]
pub fn plane_wave_particle_velocity(p_rms: Pressure, rho_c: f64) -> ParticleVelocity {
    ParticleVelocity::new(p_rms.get() / rho_c)
}

/// Acoustic power radiated through an area `a` with uniform intensity `i`:
/// `W = I·A`.
#[must_use]
pub fn power_from_uniform_intensity(intensity: Intensity, area_m2: f64) -> Power {
    Power::new(intensity.get() * area_m2)
}

/// Wavelength of a harmonic wave: `λ = c / f`.
#[must_use]
pub fn wavelength(speed_of_sound: f64, frequency: Frequency) -> f64 {
    debug_assert!(frequency.hz() > 0.0);
    speed_of_sound / frequency.hz()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::wave::Frequency;

    const TOL: f64 = 1e-9;

    #[test]
    fn one_pascal_is_about_94_db() {
        let spl = SoundPressureLevel::from_pressure_rms(Pressure::new(1.0));
        assert!((spl.db() - 20.0 * libm::log10(1.0 / 20.0e-6)).abs() < TOL);
        assert!((spl.db() - 93.979_400_086_720_37).abs() < 1e-9);
    }

    #[test]
    fn spl_try_variants_reject_non_positive_input() {
        assert!(SoundPressureLevel::try_from_pressure_rms(Pressure::new(0.0)).is_err());
        assert!(SoundPressureLevel::try_from_pressure_rms(Pressure::new(-1.0)).is_err());
        assert!(SoundPressureLevel::try_from_pressure_rms(Pressure::new(f64::NAN)).is_err());
        assert!(SoundPressureLevel::try_from_intensity(Intensity::new(0.0)).is_err());
        assert!(SoundPressureLevel::try_from_power(Power::new(-1.0)).is_err());
        // Valid input matches the panicking constructor.
        let ok = SoundPressureLevel::try_from_pressure_rms(Pressure::new(1.0)).unwrap();
        assert_eq!(
            ok,
            SoundPressureLevel::from_pressure_rms(Pressure::new(1.0))
        );
    }

    #[test]
    fn spl_round_trips_through_pressure() {
        let spl = SoundPressureLevel::from_pressure_rms(Pressure::new(0.05));
        let p = spl.to_pressure_rms();
        assert!((p.get() - 0.05).abs() < 1e-12);
    }

    #[test]
    fn two_equal_levels_add_three_db() {
        let a = SoundPressureLevel::from_pressure_rms(Pressure::new(0.1));
        let b = SoundPressureLevel::from_pressure_rms(Pressure::new(0.1));
        let sum = SoundPressureLevel::combine(&[a, b]).unwrap();
        assert!((sum.db() - (a.db() + 10.0 * libm::log10(2.0))).abs() < TOL);
    }

    #[test]
    fn average_of_two_levels_is_minus_three_db() {
        let a = SoundPressureLevel::from_pressure_rms(Pressure::new(0.1));
        let avg = SoundPressureLevel::average(&[a, a]).unwrap();
        assert!((avg.db() - a.db()).abs() < TOL);
    }

    #[test]
    fn intensity_pressure_round_trip() {
        let rho_c = 413.2;
        let p = Pressure::new(2.5);
        let i = intensity_from_pressure(p, rho_c);
        let p2 = pressure_from_intensity(i, rho_c);
        assert!((p.get() - p2.get()).abs() < TOL);
    }

    #[test]
    fn plane_wave_velocity_matches_impedance() {
        let u = plane_wave_particle_velocity(Pressure::new(1.0), 413.2);
        assert!((u.get() - 1.0 / 413.2).abs() < TOL);
    }

    #[test]
    fn quantity_units_compose() {
        // p * u -> intensity; p / u -> impedance; mixing is a compile error.
        let p = Pressure::new(1.0);
        let u = ParticleVelocity::new(2.0);
        let i: Intensity = p * u;
        assert_eq!(i.get(), 2.0);
        let z: SpecificImpedanceScalar = p / u;
        assert_eq!(z.get(), 0.5);
    }

    #[test]
    fn power_levels_round_trip() {
        let lw = SoundPressureLevel::from_power(Power::new(1.0)); // re 1 pW
        assert!((lw.db() - 120.0).abs() < TOL);
        assert!((lw.to_power().get() - 1.0).abs() < TOL);
    }

    #[test]
    fn wavelength_is_c_over_f() {
        let lam = wavelength(343.0, Frequency::hertz(1000.0).unwrap());
        assert!((lam - 0.343).abs() < TOL);
    }
}
