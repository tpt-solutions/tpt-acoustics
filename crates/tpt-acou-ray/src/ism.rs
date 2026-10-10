//! Image Source Method (ISM) for convex polyhedral (shoebox) rooms.
//!
//! For a rectangular room `[0,Lx]×[0,Ly]×[0,Lz]` the exact early-reflection
//! series is a finite lattice of mirrored sources: per axis, mirror room
//! `k ∈ ℤ` holds the image at `k·L + (s if k even else L−s)` with exactly
//! `|k|` wall hits along that axis, ordered `Low, High, Low, …` for negative
//! `k` and `High, Low, High, …` for positive `k`. A total reflection order
//! cap `|kx|+|ky|+|kz| ≤ max_order` bounds the enumeration (the count never
//! exceeds `(2·max_order+1)³` — a bound we assert and Kani-prove).
//!
//! Energy is handled at pressure-amplitude level: each hit contributes the
//! real amplitude factor `|r| = sqrt(1 − alpha_wall(f))` (so `|r|² = 1 − α`
//! conserves energy), and propagation uses the free-field Green's function.
//! Phase reflection (pressure-release sign flips) is deliberately not
//! modeled; walls are treated as locally reacting with real, positive
//! reflection factors.

use alloc::vec;
use alloc::vec::Vec;

use tpt_acou_core::error::AcousticError;
use tpt_acou_core::medium::Medium;
use tpt_acou_core::wave::free_field_green;
use tpt_math_linalg_complex::Complex;

use crate::material::BoundaryMaterial;
use crate::vec3::Vec3;

/// Wall indices of a shoebox room: `-X, +X, -Y, +Y, -Z, +Z`.
pub const WALL_LOW_X: usize = 0;
/// See [`WALL_LOW_X`].
pub const WALL_HIGH_X: usize = 1;
/// See [`WALL_LOW_X`].
pub const WALL_LOW_Y: usize = 2;
/// See [`WALL_LOW_X`].
pub const WALL_HIGH_Y: usize = 3;
/// See [`WALL_LOW_X`].
pub const WALL_LOW_Z: usize = 4;
/// See [`WALL_LOW_X`].
pub const WALL_HIGH_Z: usize = 5;

/// A rectangular (shoebox) room with per-wall materials.
#[derive(Debug, Clone)]
pub struct ShoeboxRoom {
    dims: [f64; 3],
    walls: [BoundaryMaterial; 6],
}

impl ShoeboxRoom {
    /// Creates a shoebox room from dimensions (metres) and the six wall
    /// materials in [`WALL_LOW_X`]…[`WALL_HIGH_Z`] order.
    ///
    /// # Errors
    /// Returns [`AcousticError::NonPositive`] for any non-positive dimension.
    pub fn new(dims: [f64; 3], walls: [BoundaryMaterial; 6]) -> Result<Self, AcousticError> {
        for &d in dims.iter() {
            if d <= 0.0 {
                return Err(AcousticError::NonPositive("room dimension"));
            }
        }
        Ok(Self { dims, walls })
    }

    /// A shoebox with identical material on all six walls.
    ///
    /// # Errors
    /// See [`ShoeboxRoom::new`].
    pub fn uniform(dims: [f64; 3], wall: BoundaryMaterial) -> Result<Self, AcousticError> {
        Self::new(
            dims,
            [
                wall.clone(),
                wall.clone(),
                wall.clone(),
                wall.clone(),
                wall.clone(),
                wall,
            ],
        )
    }

    /// Room dimensions `[Lx, Ly, Lz]` in metres.
    #[must_use]
    pub const fn dims(&self) -> [f64; 3] {
        self.dims
    }

    /// The wall material for a wall index (`0..6`).
    #[must_use]
    pub fn wall(&self, index: usize) -> &BoundaryMaterial {
        &self.walls[index]
    }

    /// Volume in m³.
    #[must_use]
    pub fn volume(&self) -> f64 {
        self.dims[0] * self.dims[1] * self.dims[2]
    }

    /// Total boundary surface area in m².
    #[must_use]
    pub fn surface_area(&self) -> f64 {
        let [lx, ly, lz] = self.dims;
        2.0 * (lx * ly + ly * lz + lx * lz)
    }

    /// Strict interior test (boundary points are *outside* — ISM and ray
    /// tracing require source and receiver strictly inside).
    #[must_use]
    pub fn contains(&self, p: Vec3) -> bool {
        p.x > 0.0
            && p.x < self.dims[0]
            && p.y > 0.0
            && p.y < self.dims[1]
            && p.z > 0.0
            && p.z < self.dims[2]
    }

    /// Total absorption `A = Σ α_i(f)·S_i` in m² sabins at `freq_hz`.
    #[must_use]
    pub fn total_absorption(&self, freq_hz: f64) -> f64 {
        let [lx, ly, lz] = self.dims;
        let areas = [ly * lz, ly * lz, lx * lz, lx * lz, lx * ly, lx * ly];
        areas
            .iter()
            .zip(self.walls.iter())
            .map(|(&s, w)| w.absorption_at(freq_hz) * s)
            .sum()
    }

    /// Mean absorption coefficient `ᾱ = A / S`.
    #[must_use]
    pub fn mean_absorption(&self, freq_hz: f64) -> f64 {
        self.total_absorption(freq_hz) / self.surface_area()
    }

    /// Sabine reverberation time `T = 0.161·V / A` in seconds.
    #[must_use]
    pub fn sabine_rt60(&self, freq_hz: f64) -> f64 {
        0.161 * self.volume() / self.total_absorption(freq_hz)
    }

    /// Eyring reverberation time
    /// `T = 0.161·V / (−S·ln(1 − ᾱ))` in seconds (more accurate for high
    /// absorption; returns infinity for perfectly rigid rooms).
    #[must_use]
    pub fn eyring_rt60(&self, freq_hz: f64) -> f64 {
        let a_bar = self.mean_absorption(freq_hz);
        if a_bar >= 1.0 {
            return f64::INFINITY;
        }
        0.161 * self.volume() / (-self.surface_area() * libm::log(1.0 - a_bar))
    }

    /// All image sources with total reflection order `≤ max_order`, each with
    /// its per-frequency-independent path data.
    ///
    /// # Errors
    /// Returns [`AcousticError::InvalidGeometry`] unless `src` lies strictly
    /// inside the room.
    pub fn image_sources(
        &self,
        src: Vec3,
        max_order: u32,
    ) -> Result<Vec<ImagePath>, AcousticError> {
        if !self.contains(src) {
            return Err(AcousticError::InvalidGeometry(
                "ISM source must lie strictly inside the room",
            ));
        }
        let mut out = Vec::new();
        // Per-axis mirror-room indices are bounded by max_order on each axis
        // because the total order |kx|+|ky|+|kz| caps each |k*|.
        let n = max_order as i64;
        for kz in -n..=n {
            for ky in -n..=n {
                for kx in -n..=n {
                    let order = (kx.abs() + ky.abs() + kz.abs()) as u32;
                    if order > max_order {
                        continue;
                    }
                    let position = Vec3::new(
                        self.axis_image_position(src.x, self.dims[0], kx),
                        self.axis_image_position(src.y, self.dims[1], ky),
                        self.axis_image_position(src.z, self.dims[2], kz),
                    );
                    let mut walls = Vec::new();
                    self.axis_walls(kx, WALL_LOW_X, WALL_HIGH_X, &mut walls);
                    self.axis_walls(ky, WALL_LOW_Y, WALL_HIGH_Y, &mut walls);
                    self.axis_walls(kz, WALL_LOW_Z, WALL_HIGH_Z, &mut walls);
                    debug_assert_eq!(walls.len(), order as usize);
                    out.push(ImagePath {
                        position,
                        order,
                        walls,
                    });
                }
            }
        }
        Ok(out)
    }

    /// Image coordinate along one axis for mirror-room index `k`.
    fn axis_image_position(&self, s: f64, l: f64, k: i64) -> f64 {
        if k % 2 == 0 {
            k as f64 * l + s
        } else {
            k as f64 * l + (l - s)
        }
    }

    /// Wall-hit sequence along one axis for mirror-room index `k`, in
    /// traversal order. Length equals `|k|`.
    fn axis_walls(&self, k: i64, low: usize, high: usize, out: &mut Vec<usize>) {
        if k == 0 {
            return;
        }
        if k > 0 {
            // Boundaries at l, 2l, ..., k*l: odd multiples map to the high
            // wall, even multiples to the low wall.
            for m in 1..=k {
                out.push(if m % 2 == 1 { high } else { low });
            }
        } else {
            // Boundaries at 0, -l, ..., (k+1)*l: traversal order is
            // Low, High, Low, High, ... (|k| hits).
            for m in 0..(-k) {
                out.push(if m % 2 == 0 { low } else { high });
            }
        }
    }

    /// Frequency-domain transfer function (ISM sum)
    /// `H(f) = Σ_img rc(f)·G(k(f), d)`, where `rc(f)` is the product of
    /// wall amplitude reflection factors along the image path. Frequencies
    /// are taken from `freqs_hz` (plain values for kernel simplicity).
    ///
    /// Images closer than `min_distance` to the receiver (including the
    /// degenerate `src == rcv` direct path) are skipped — the free-field
    /// Green's function is singular there.
    ///
    /// # Errors
    /// Propagates [`ShoeboxRoom::image_sources`] and
    /// [`tpt_acou_core::wave::free_field_green`] errors.
    pub fn transfer_function(
        &self,
        src: Vec3,
        rcv: Vec3,
        freqs_hz: &[f64],
        max_order: u32,
        medium: &Medium,
    ) -> Result<Vec<Complex>, AcousticError> {
        if !self.contains(rcv) {
            return Err(AcousticError::InvalidGeometry(
                "ISM receiver must lie strictly inside the room",
            ));
        }
        let images = self.image_sources(src, max_order)?;
        let mut out = vec![Complex::new(0.0, 0.0); freqs_hz.len()];
        for f_idx in 0..freqs_hz.len() {
            let f = freqs_hz[f_idx];
            let k = medium
                .wavenumber(tpt_acou_core::Frequency::hertz(f)?)?
                .rad_per_m();
            let mut sum = Complex::new(0.0, 0.0);
            for img in &images {
                let d = img.position.distance(rcv);
                if d < 1e-9 {
                    continue;
                }
                let rc: f64 = img
                    .walls
                    .iter()
                    .map(|&w| self.walls[w].amplitude_reflection(f))
                    .product();
                let g = free_field_green(k, d)?;
                sum = sum + g * rc;
            }
            out[f_idx] = sum;
        }
        Ok(out)
    }

    /// Energetic echogram: `(arrival_time_s, energy_weight)` pairs for all
    /// images up to `max_order`, where the weight is `rc² / (4π d)²` —
    /// proportional to received intensity. Includes the direct sound
    /// (`order 0`).
    ///
    /// # Errors
    /// Propagates [`ShoeboxRoom::image_sources`].
    pub fn echogram(
        &self,
        src: Vec3,
        rcv: Vec3,
        freq_hz: f64,
        max_order: u32,
    ) -> Result<Vec<(f64, f64)>, AcousticError> {
        if !self.contains(rcv) {
            return Err(AcousticError::InvalidGeometry(
                "ISM receiver must lie strictly inside the room",
            ));
        }
        let images = self.image_sources(src, max_order)?;
        let mut out = Vec::with_capacity(images.len());
        for img in &images {
            let d = img.position.distance(rcv);
            if d < 1e-9 {
                continue;
            }
            let rc: f64 = img
                .walls
                .iter()
                .map(|&w| self.walls[w].amplitude_reflection(freq_hz))
                .product();
            let denom = 4.0 * core::f64::consts::PI * d;
            out.push((d / SPEED_OF_SOUND_REFERENCE, rc * rc / (denom * denom)));
        }
        out.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(core::cmp::Ordering::Equal));
        Ok(out)
    }
}

/// Reference speed of sound for echogram timing (m/s). Timing callers that
/// need another medium can rescale times by `c_ref / c`.
pub const SPEED_OF_SOUND_REFERENCE: f64 = 343.0;

/// One image source: mirrored position, total reflection order, and the wall
/// indices hit in traversal order (length equals `order`).
#[derive(Debug, Clone, PartialEq)]
pub struct ImagePath {
    /// Mirrored source position (metres).
    pub position: Vec3,
    /// Total reflection order `|kx|+|ky|+|kz|`.
    pub order: u32,
    /// Wall indices hit in traversal order.
    pub walls: Vec<usize>,
}

/// Upper bound on the number of enumerated lattice points for a given order
/// cap: `(2·max_order + 1)³`. Used by tests and the Kani proof that the
/// enumeration never exceeds its bound.
#[must_use]
pub fn image_count_upper_bound(max_order: u32) -> u64 {
    let n = u64::from(max_order) * 2 + 1;
    n * n * n
}

/// Exact number of lattice points `(kx, ky, kz)` with
/// `|kx| + |ky| + |kz| ≤ max_order` — the image count of a full ISM
/// enumeration before any singularity filtering. This is the arithmetic the
/// order-bounded loop performs; the Kani proof in `tpt-acou-verify` checks it
/// against [`image_count_upper_bound`] for symbolic orders.
#[must_use]
pub fn count_images_with_order_at_most(max_order: u32) -> u64 {
    let n = max_order as i64;
    let mut count = 0_u64;
    let mut kz = -n;
    while kz <= n {
        let mut ky = -n;
        while ky <= n {
            let mut kx = -n;
            while kx <= n {
                if kx.abs() + ky.abs() + kz.abs() <= n {
                    count += 1;
                }
                kx += 1;
            }
            ky += 1;
        }
        kz += 1;
    }
    count
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

    fn src() -> Vec3 {
        Vec3::new(1.0, 1.0, 1.5)
    }

    #[test]
    fn lattice_count_matches_bruteforce_and_bound() {
        for n in 0..6_u32 {
            let count = count_images_with_order_at_most(n);
            assert!(count <= image_count_upper_bound(n), "n={n}");
        }
        assert_eq!(count_images_with_order_at_most(0), 1);
        // Shell orders: order ≤1 -> 1 + 6 = 7.
        assert_eq!(count_images_with_order_at_most(1), 7);
        // order ≤2 -> 1 + 6 + 18 = 25.
        assert_eq!(count_images_with_order_at_most(2), 25);
    }

    #[test]
    fn order_zero_is_direct_sound_only() {
        let r = room();
        let imgs = r.image_sources(src(), 0).unwrap();
        assert_eq!(imgs.len(), 1);
        assert_eq!(imgs[0].order, 0);
        assert!(imgs[0].walls.is_empty());
        assert_eq!(imgs[0].position, src());
    }

    #[test]
    fn count_is_bounded_and_orders_are_capped() {
        let r = room();
        let max_order = 3;
        let imgs = r.image_sources(src(), max_order).unwrap();
        assert!(imgs.len() as u64 <= image_count_upper_bound(max_order));
        assert!(imgs.iter().all(|i| i.order <= max_order));
        assert!(imgs.iter().all(|i| i.walls.len() == i.order as usize));
    }

    #[test]
    fn first_order_images_mirrored_across_walls() {
        let r = room();
        let s = src();
        let imgs = r.image_sources(s, 1).unwrap();
        // Expected first-order positions: single-axis mirrors.
        let expected = [
            Vec3::new(-s.x, s.y, s.z),       // low x
            Vec3::new(10.0 - s.x, s.y, s.z), // high x: 2*Lx - s
            Vec3::new(s.x, -s.y, s.z),       // low y
            Vec3::new(s.x, 8.0 - s.y, s.z),  // high y
            Vec3::new(s.x, s.y, -s.z),       // low z
            Vec3::new(s.x, s.y, 6.0 - s.z),  // high z
        ];
        for e in expected {
            assert!(
                imgs.iter()
                    .any(|i| i.order == 1 && (i.position - e).norm() < 1e-12),
                "missing image at {e:?}"
            );
        }
    }

    #[test]
    fn wall_sequences_match_geometry() {
        let r = room();
        let s = src();
        let imgs = r.image_sources(s, 2).unwrap();
        // Image at (2Lx - s.x) has ONE high-x hit: 2*1*Lx - s with k=+1.
        let img = imgs
            .iter()
            .find(|i| i.order == 1 && (i.position.x - (10.0 - s.x)).abs() < 1e-12)
            .unwrap();
        assert_eq!(img.walls, vec![WALL_HIGH_X]);
        // Image at (2Lx + s.x) has k=+2: hits High then Low.
        let img2 = imgs
            .iter()
            .find(|i| i.order == 2 && (i.position.x - (10.0 + s.x)).abs() < 1e-12)
            .unwrap();
        assert_eq!(img2.walls, vec![WALL_HIGH_X, WALL_LOW_X]);
        // Image at (-s.x) has k=-1: single Low hit.
        let img3 = imgs
            .iter()
            .find(|i| i.order == 1 && (i.position.x + s.x).abs() < 1e-12)
            .unwrap();
        assert_eq!(img3.walls, vec![WALL_LOW_X]);
    }

    #[test]
    fn transfer_function_direct_term_matches_free_field() {
        // With rigid walls and max_order = 0 the ISM sum reduces to the
        // free-field Green's function.
        let r = ShoeboxRoom::uniform([5.0, 4.0, 3.0], BoundaryMaterial::rigid()).unwrap();
        let s = src();
        let g = Vec3::new(3.0, 2.0, 2.0);
        let medium = Medium::air_20c();
        let h = r.transfer_function(s, g, &[1000.0], 0, &medium).unwrap();
        let d = s.distance(g);
        let expected = free_field_green(
            medium
                .wavenumber(tpt_acou_core::Frequency::hertz(1000.0).unwrap())
                .unwrap()
                .rad_per_m(),
            d,
        )
        .unwrap();
        assert!((h[0] - expected).norm() < 1e-15);
    }

    #[test]
    fn echogram_arrival_times_sorted_direct_first() {
        let r = room();
        let s = src();
        let g = Vec3::new(3.0, 2.5, 2.0);
        let echo = r.echogram(s, g, 1000.0, 2).unwrap();
        assert_eq!(echo[0].0, s.distance(g) / SPEED_OF_SOUND_REFERENCE);
        for w in echo.windows(2) {
            assert!(w[0].0 <= w[1].0 + 1e-15);
            assert!(w[1].1 >= 0.0);
        }
    }

    #[test]
    fn sabine_and_eyring_agree_for_small_alpha() {
        let r = room();
        let t_sab = r.sabine_rt60(1000.0);
        let t_ey = r.eyring_rt60(1000.0);
        let a_bar = r.mean_absorption(1000.0);
        assert!(a_bar < 0.2);
        assert!((t_sab - t_ey).abs() / t_sab < 0.1, "sab {t_sab} eyr {t_ey}");
        assert!(t_sab > 0.0 && t_sab < 10.0);
    }

    #[test]
    fn rejects_outside_points() {
        let r = room();
        assert!(r.image_sources(Vec3::new(-1.0, 1.0, 1.0), 2).is_err());
        assert!(r.image_sources(src(), 2).is_ok());
        assert!(r
            .transfer_function(
                src(),
                Vec3::new(0.0, 1.0, 1.0),
                &[100.0],
                1,
                &Medium::air_20c()
            )
            .is_err());
    }

    #[test]
    fn transfer_magnitude_respects_triangle_bound() {
        // |H| <= sum_images 1/(4*pi*d): reflected terms may cancel the direct
        // sound pointwise, but the magnitude can never exceed the sum of the
        // path magnitudes.
        let r = ShoeboxRoom::uniform([5.0, 4.0, 3.0], BoundaryMaterial::rigid()).unwrap();
        let s = src();
        let g = Vec3::new(3.0, 2.0, 2.0);
        let medium = Medium::air_20c();
        let freq = 700.0;
        let h = r.transfer_function(s, g, &[freq], 3, &medium).unwrap();
        let mut bound = 0.0;
        for img in r.image_sources(s, 3).unwrap() {
            let d = img.position.distance(g);
            if d < 1e-9 {
                continue;
            }
            bound += 1.0 / (4.0 * core::f64::consts::PI * d);
        }
        assert!(
            h[0].norm() <= bound * (1.0 + 1e-12),
            "|H| {} > bound {bound}",
            h[0].norm()
        );
    }
}
