//! Differentiable ray-tracing hooks (forward-mode AD over room geometry).
//!
//! The ISM pressure sum is an analytic function of the room dimensions, so a
//! one-parameter forward-mode dual number propagates the exact derivative
//! `∂|H|²/∂L_axis` through the whole sum — the gradient an inverse-design
//! agent needs to fit a room's geometry to a target response by gradient
//! descent (spec.txt §2, "differentiable for neural acoustic rendering and
//! inverse design").
//!
//! # Replacing the two-mirror lattice
//!
//! Along axis `a` with length `L`, the image in mirror room `k` sits at
//! `k·L + s` (even `k`) or `(k+1)·L − s` (odd `k`) — both *linear* in `L`,
//! with derivatives `k` and `k+1`. Distances, the complex exponential and
//! the `1/(4πd)` factor then flow through dual arithmetic unchanged.

use crate::ism::ShoeboxRoom;
use crate::vec3::Vec3;

/// Forward-mode dual number `v + d·ε` (`ε² = 0`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Dual {
    /// Value.
    pub v: f64,
    /// Derivative w.r.t. the single active parameter.
    pub d: f64,
}

impl Dual {
    /// Constant (derivative 0).
    #[must_use]
    pub const fn c(v: f64) -> Self {
        Self { v, d: 0.0 }
    }

    /// Active variable (derivative 1).
    #[must_use]
    pub const fn var(v: f64) -> Self {
        Self { v, d: 1.0 }
    }

    /// Dual with an explicit seed derivative.
    #[must_use]
    pub const fn new(v: f64, d: f64) -> Self {
        Self { v, d }
    }

    /// Square root (chain rule `d = d(v)/(2√v)`).
    #[must_use]
    pub fn sqrt(self) -> Self {
        let s = libm::sqrt(self.v);
        Self::new(s, self.d / (2.0 * s))
    }
}

impl core::ops::Add for Dual {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(self.v + o.v, self.d + o.d)
    }
}

impl core::ops::Sub for Dual {
    type Output = Self;
    fn sub(self, o: Self) -> Self {
        Self::new(self.v - o.v, self.d - o.d)
    }
}

impl core::ops::Mul for Dual {
    type Output = Self;
    fn mul(self, o: Self) -> Self {
        Self::new(self.v * o.v, self.d * o.v + self.v * o.d)
    }
}

impl core::ops::Div for Dual {
    type Output = Self;
    fn div(self, o: Self) -> Self {
        Self::new(self.v / o.v, (self.d * o.v - self.v * o.d) / (o.v * o.v))
    }
}

impl core::ops::Neg for Dual {
    type Output = Self;
    fn neg(self) -> Self {
        Self::new(-self.v, -self.d)
    }
}

/// Complex dual number `re + i·im`, both parts duals.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ComplexDual {
    /// Real part.
    pub re: Dual,
    /// Imaginary part.
    pub im: Dual,
}

impl ComplexDual {
    /// Builds from dual parts.
    #[must_use]
    pub const fn new(re: Dual, im: Dual) -> Self {
        Self { re, im }
    }

    /// Complex exponential `exp(a + i·b) = e^a·(cos b + i·sin b)`.
    #[must_use]
    pub fn exp(self) -> Self {
        let ea = libm::exp(self.re.v);
        Self {
            re: Dual::new(
                ea * libm::cos(self.im.v),
                ea * (self.re.d * libm::cos(self.im.v) - self.im.d * libm::sin(self.im.v)),
            ),
            im: Dual::new(
                ea * libm::sin(self.im.v),
                ea * (self.re.d * libm::sin(self.im.v) + self.im.d * libm::cos(self.im.v)),
            ),
        }
    }

    /// Squared magnitude `|z|² = re² + im²` (dual).
    #[must_use]
    pub fn norm_sqr(self) -> Dual {
        self.re * self.re + self.im * self.im
    }
}

impl core::ops::Add for ComplexDual {
    type Output = Self;
    fn add(self, o: Self) -> Self {
        Self::new(self.re + o.re, self.im + o.im)
    }
}

impl core::ops::Mul<Dual> for ComplexDual {
    type Output = Self;
    fn mul(self, s: Dual) -> Self {
        Self::new(self.re * s, self.im * s)
    }
}

impl core::ops::Div<Dual> for ComplexDual {
    type Output = Self;
    fn div(self, s: Dual) -> Self {
        // z / s = z * conj(s) / |s|², s real-valued dual.
        let inv = Dual::c(1.0) / s;
        self * inv
    }
}

/// ISM pressure sum with the room length `L_axis` treated as the active
/// dual parameter. Returns `(|H|², ∂|H|²/∂L_axis)` exactly (to floating
/// point), for the image set up to `max_order`.
///
/// Walls contribute their real amplitude factor `sqrt(1 − α(f))` per hit
/// (independent of geometry); propagation is the free-field Green's
/// function with the `exp(+iωt)` convention.
///
/// # Errors
/// Returns [`tpt_acou_core::AcousticError::InvalidGeometry`] unless both
/// points lie strictly inside the room.
pub fn ism_energy_and_gradient(
    room: &ShoeboxRoom,
    src: Vec3,
    rcv: Vec3,
    freq_hz: f64,
    speed_of_sound: f64,
    max_order: u32,
    axis: usize,
) -> Result<(f64, f64), tpt_acou_core::AcousticError> {
    use tpt_acou_core::AcousticError;
    if !room.contains(src) || !room.contains(rcv) {
        return Err(AcousticError::InvalidGeometry(
            "gradient evaluation requires source and receiver strictly inside",
        ));
    }
    let dims = room.dims();
    let k = core::f64::consts::TAU * freq_hz / speed_of_sound;
    let images = room.image_sources(src, max_order)?;

    let mut acc = ComplexDual::new(Dual::c(0.0), Dual::c(0.0));
    for img in &images {
        let d_geom = img.position.distance(rcv);
        if d_geom < 1e-9 {
            continue;
        }
        // Image coordinate along the active axis, linear in L.
        let coord = axis_coord_dual(src, dims[axis], rcv, img, axis);
        let pos_axes = [img.position.x, img.position.y, img.position.z];
        let rcv_axes = [rcv.x, rcv.y, rcv.z];
        let mut comp = [Dual::c(0.0); 3];
        for (a, item) in comp.iter_mut().enumerate() {
            *item = if a == axis {
                coord
            } else {
                Dual::c(pos_axes[a] - rcv_axes[a])
            };
        }
        // Distance² = Σ comp² (component duals carry ∂/∂L).
        let dist2 = comp[0] * comp[0] + comp[1] * comp[1] + comp[2] * comp[2];
        let dist = dist2.sqrt();
        // rc(f): product of per-hit amplitude reflection factors.
        let rc: f64 = img
            .walls
            .iter()
            .map(|&w| room.wall(w).amplitude_reflection(freq_hz))
            .product();
        // G = exp(-i k d) / (4 π d)
        let phase = ComplexDual::new(Dual::c(0.0), Dual::c(-k) * dist);
        let green = phase.exp() / (Dual::c(4.0 * core::f64::consts::PI) * dist);
        acc = acc + green * Dual::c(rc);
    }
    let energy = acc.norm_sqr();
    Ok((energy.v, energy.d))
}

/// Dual-valued image coordinate offset `position_axis − rcv_axis` along
/// `axis`, carrying `∂/∂L` as its derivative. Mirrors the ISM lattice: even
/// mirror room `k` gives `k·L + s` (derivative `k`), odd gives
/// `(k+1)·L − s` (derivative `k+1`).
fn axis_coord_dual(src: Vec3, l: f64, rcv: Vec3, img: &crate::ism::ImagePath, axis: usize) -> Dual {
    let src_axes = [src.x, src.y, src.z];
    let rcv_axes = [rcv.x, rcv.y, rcv.z];
    let s = src_axes[axis];
    let pos = [img.position.x, img.position.y, img.position.z][axis];
    // Recover the mirror-room index from the mirrored position: even chain
    // pos = k*L + s => k = (pos-s)/L; odd chain pos = (k+1)*L - s.
    let k_even = libm::round((pos - s) / l);
    let k1 = libm::round((pos + s) / l);
    let even_pos = k_even * l + s;
    if (even_pos - pos).abs() < 1e-9 * l.max(1.0) {
        Dual::new(pos - rcv_axes[axis], k_even)
    } else {
        Dual::new(pos - rcv_axes[axis], k1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::BoundaryMaterial;

    fn room() -> ShoeboxRoom {
        ShoeboxRoom::new(
            [5.0, 4.0, 3.0],
            [
                BoundaryMaterial::concrete_painted(),
                BoundaryMaterial::concrete_painted(),
                BoundaryMaterial::curtain_heavy(),
                BoundaryMaterial::concrete_painted(),
                BoundaryMaterial::glass(),
                BoundaryMaterial::rigid(),
            ],
        )
        .unwrap()
    }

    /// Energy as a plain function of the axis length, for finite differences.
    fn energy_at(room: &ShoeboxRoom, src: Vec3, rcv: Vec3, axis: usize, l: f64) -> f64 {
        let mut dims = room.dims();
        dims[axis] = l;
        let stretched = ShoeboxRoom::new(
            dims,
            [
                room.wall(0).clone(),
                room.wall(1).clone(),
                room.wall(2).clone(),
                room.wall(3).clone(),
                room.wall(4).clone(),
                room.wall(5).clone(),
            ],
        )
        .unwrap();
        let (e, _) = ism_energy_and_gradient(
            &stretched,
            src,
            rcv,
            500.0,
            crate::ism::SPEED_OF_SOUND_REFERENCE,
            2,
            axis,
        )
        .unwrap();
        e
    }

    #[test]
    fn gradient_matches_finite_differences() {
        let room = room();
        let src = Vec3::new(1.2, 0.9, 1.1);
        let rcv = Vec3::new(3.4, 2.6, 1.9);
        for axis in 0..3 {
            let (e, g) = ism_energy_and_gradient(
                &room,
                src,
                rcv,
                500.0,
                crate::ism::SPEED_OF_SOUND_REFERENCE,
                2,
                axis,
            )
            .unwrap();
            assert!(e > 0.0);
            let h = 1e-6;
            let fd = (energy_at(&room, src, rcv, axis, room.dims()[axis] + h)
                - energy_at(&room, src, rcv, axis, room.dims()[axis] - h))
                / (2.0 * h);
            let scale = g.abs().max(fd.abs()).max(1e-30);
            assert!(
                (g - fd).abs() / scale < 1e-4,
                "axis {axis}: dual {g} vs fd {fd}"
            );
        }
    }

    #[test]
    fn dual_arithmetic_chain_rule() {
        // d/dx sqrt(x²+1) at x=2 is 2/sqrt(5).
        let x = Dual::var(2.0);
        let y = (x * x + Dual::c(1.0)).sqrt();
        assert!((y.v - libm::sqrt(5.0)).abs() < 1e-12);
        assert!((y.d - 2.0 / libm::sqrt(5.0)).abs() < 1e-12);
    }

    #[test]
    fn complex_dual_exp_matches_reference() {
        // exp(i·π/3) = 0.5 + i·(√3/2); with respect to θ, derivative is
        // i·exp(iθ).
        let theta = Dual::var(core::f64::consts::FRAC_PI_3);
        let z = ComplexDual::new(Dual::c(0.0), theta).exp();
        assert!((z.re.v - 0.5).abs() < 1e-12);
        assert!((z.im.v - libm::sqrt(3.0) / 2.0).abs() < 1e-12);
        // d/dθ exp(iθ) = i·exp(iθ): d(re) = -im, d(im) = re.
        assert!((z.re.d + z.im.v).abs() < 1e-12);
        assert!((z.im.d - z.re.v).abs() < 1e-12);
    }

    #[test]
    fn rejects_outside_geometry() {
        let room = room();
        assert!(ism_energy_and_gradient(
            &room,
            Vec3::new(-1.0, 1.0, 1.0),
            Vec3::new(2.0, 2.0, 1.5),
            500.0,
            343.0,
            1,
            0
        )
        .is_err());
    }
}
