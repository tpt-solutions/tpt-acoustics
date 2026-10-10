//! Minimal 3-vector algebra for geometric acoustics.
//!
//! In-house by design: the ray kernels only need dot/cross/normalize and
//! reflections, and keeping the type local keeps the crate `no_std`, `unsafe`-
//! free and dependency-light. (Conversions to/from `tpt-math-geometry` types
//! live with the consumers that need them.)

use core::ops::{Add, AddAssign, Div, Mul, Neg, Sub};

use libm::sqrt;

/// A 3-component double-precision vector (metres unless documented otherwise).
#[derive(Debug, Clone, Copy, PartialEq, Default)]
pub struct Vec3 {
    /// x component.
    pub x: f64,
    /// y component.
    pub y: f64,
    /// z component.
    pub z: f64,
}

impl Vec3 {
    /// The zero vector.
    pub const ZERO: Vec3 = Vec3 {
        x: 0.0,
        y: 0.0,
        z: 0.0,
    };

    /// Builds a vector from components.
    #[must_use]
    pub const fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }

    /// Euclidean dot product.
    #[must_use]
    pub fn dot(self, o: Self) -> f64 {
        self.x * o.x + self.y * o.y + self.z * o.z
    }

    /// Euclidean cross product `self × o`.
    #[must_use]
    pub fn cross(self, o: Self) -> Self {
        Self {
            x: self.y * o.z - self.z * o.y,
            y: self.z * o.x - self.x * o.z,
            z: self.x * o.y - self.y * o.x,
        }
    }

    /// Squared norm.
    #[must_use]
    pub fn norm_sqr(self) -> f64 {
        self.dot(self)
    }

    /// Euclidean norm.
    #[must_use]
    pub fn norm(self) -> f64 {
        sqrt(self.norm_sqr())
    }

    /// Unit vector in the same direction, or [`Vec3::ZERO`] for the zero
    /// vector (no panic, safe for Kani harnesses).
    #[must_use]
    pub fn normalized(self) -> Self {
        let n = self.norm();
        if n == 0.0 {
            Self::ZERO
        } else {
            self / n
        }
    }

    /// Specular reflection of a direction `d` about a surface normal `n`
    /// (`n` must be unit length): `d − 2·(d·n)·n`.
    #[must_use]
    pub fn reflect(self, n: Self) -> Self {
        self - n * (2.0 * self.dot(n))
    }

    /// Componentwise minimum.
    #[must_use]
    pub fn min_components(self, o: Self) -> Self {
        Self {
            x: self.x.min(o.x),
            y: self.y.min(o.y),
            z: self.z.min(o.z),
        }
    }

    /// Componentwise maximum.
    #[must_use]
    pub fn max_components(self, o: Self) -> Self {
        Self {
            x: self.x.max(o.x),
            y: self.y.max(o.y),
            z: self.z.max(o.z),
        }
    }

    /// Distance between two points.
    #[must_use]
    pub fn distance(self, o: Self) -> f64 {
        (self - o).norm()
    }

    /// Linear interpolation `self + t·(o − self)`.
    #[must_use]
    pub fn lerp(self, o: Self, t: f64) -> Self {
        self + (o - self) * t
    }
}

impl Add for Vec3 {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self {
            x: self.x + o.x,
            y: self.y + o.y,
            z: self.z + o.z,
        }
    }
}

impl AddAssign for Vec3 {
    fn add_assign(&mut self, o: Self) {
        *self = *self + o;
    }
}

impl Sub for Vec3 {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self {
            x: self.x - o.x,
            y: self.y - o.y,
            z: self.z - o.z,
        }
    }
}

impl Neg for Vec3 {
    type Output = Self;
    fn neg(self) -> Self {
        Self {
            x: -self.x,
            y: -self.y,
            z: -self.z,
        }
    }
}

impl Mul<f64> for Vec3 {
    type Output = Self;
    fn mul(self, s: f64) -> Self {
        Self {
            x: self.x * s,
            y: self.y * s,
            z: self.z * s,
        }
    }
}

impl Mul<Vec3> for f64 {
    type Output = Vec3;
    fn mul(self, v: Vec3) -> Vec3 {
        v * self
    }
}

impl Div<f64> for Vec3 {
    type Output = Self;
    fn div(self, s: f64) -> Self {
        Self {
            x: self.x / s,
            y: self.y / s,
            z: self.z / s,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const TOL: f64 = 1e-12;

    #[test]
    fn dot_cross_identities() {
        let x = Vec3::new(1.0, 0.0, 0.0);
        let y = Vec3::new(0.0, 1.0, 0.0);
        let z = x.cross(y);
        assert!((z.x - 0.0).abs() < TOL && (z.y - 0.0).abs() < TOL && (z.z - 1.0).abs() < TOL);
        assert!(x.dot(y).abs() < TOL);
        assert!((x.dot(x) - 1.0).abs() < TOL);
    }

    #[test]
    fn reflection_of_axis_vector() {
        let d = Vec3::new(1.0, -1.0, 0.0);
        let n = Vec3::new(0.0, 1.0, 0.0);
        let r = d.reflect(n);
        assert!((r.x - 1.0).abs() < TOL && (r.y - 1.0).abs() < TOL);
    }

    #[test]
    fn normalized_zero_vector_is_safe() {
        assert_eq!(Vec3::ZERO.normalized(), Vec3::ZERO);
    }

    #[test]
    fn lerp_endpoints() {
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(2.0, 4.0, 6.0);
        assert_eq!(a.lerp(b, 0.0), a);
        assert_eq!(a.lerp(b, 1.0), b);
        let m = a.lerp(b, 0.5);
        assert!((m.y - 2.0).abs() < TOL);
    }
}
