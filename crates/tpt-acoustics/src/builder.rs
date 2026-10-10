//! A small builder over the shoebox workflow.
//!
//! [`RoomBuilder`] collects a box-shaped room, a wall material, a source and a
//! receiver, and validates them once. [`Simulation`] then answers the common
//! questions (echogram, impulse response, Sabine T60) without the caller
//! passing geometry around. It wraps the lower-level APIs; nothing here adds
//! physics.

use tpt_acou_core::error::AcousticError;
use tpt_acou_core::medium::Medium;
use tpt_acou_core::vec3::Vec3;
use tpt_acou_ray::ism::ShoeboxRoom;
use tpt_acou_ray::material::BoundaryMaterial;
use tpt_acou_ray::rir::rir_from_echogram;

/// Builder for a [`Simulation`] in a shoebox room.
///
/// Defaults: painted-concrete walls, air at 20 °C. Dimensions, source and
/// receiver are required.
#[derive(Debug, Clone)]
pub struct RoomBuilder {
    dimensions: Option<[f64; 3]>,
    material: BoundaryMaterial,
    source: Option<Vec3>,
    receiver: Option<Vec3>,
    medium: Medium,
}

impl Default for RoomBuilder {
    fn default() -> Self {
        Self::new()
    }
}

impl RoomBuilder {
    /// Starts a builder with the default material and medium.
    #[must_use]
    pub fn new() -> Self {
        Self {
            dimensions: None,
            material: BoundaryMaterial::concrete_painted(),
            source: None,
            receiver: None,
            medium: Medium::air_20c(),
        }
    }

    /// Room size `[length, width, height]` in metres.
    #[must_use]
    pub fn dimensions(mut self, dims: [f64; 3]) -> Self {
        self.dimensions = Some(dims);
        self
    }

    /// Wall material applied to all six faces.
    #[must_use]
    pub fn material(mut self, material: BoundaryMaterial) -> Self {
        self.material = material;
        self
    }

    /// Source position in metres, strictly inside the room.
    #[must_use]
    pub fn source(mut self, position: Vec3) -> Self {
        self.source = Some(position);
        self
    }

    /// Receiver position in metres, strictly inside the room.
    #[must_use]
    pub fn receiver(mut self, position: Vec3) -> Self {
        self.receiver = Some(position);
        self
    }

    /// Propagation medium (default air at 20 °C).
    #[must_use]
    pub fn medium(mut self, medium: Medium) -> Self {
        self.medium = medium;
        self
    }

    /// Validates the configuration and builds the simulation.
    ///
    /// # Errors
    /// Returns [`AcousticError::InvalidArgument`] if the dimensions, source or
    /// receiver are missing, and propagates the room constructor's geometry
    /// errors for non-positive dimensions. Returns
    /// [`AcousticError::InvalidGeometry`] if the source or receiver is outside
    /// the room.
    pub fn build(self) -> Result<Simulation, AcousticError> {
        let dims = self.dimensions.ok_or(AcousticError::InvalidArgument(
            "room dimensions are required",
        ))?;
        let source = self.source.ok_or(AcousticError::InvalidArgument(
            "source position is required",
        ))?;
        let receiver = self.receiver.ok_or(AcousticError::InvalidArgument(
            "receiver position is required",
        ))?;
        let room = ShoeboxRoom::uniform(dims, self.material)?;
        if !room.contains(source) {
            return Err(AcousticError::InvalidGeometry(
                "source must lie strictly inside the room",
            ));
        }
        if !room.contains(receiver) {
            return Err(AcousticError::InvalidGeometry(
                "receiver must lie strictly inside the room",
            ));
        }
        Ok(Simulation {
            room,
            source,
            receiver,
            medium: self.medium,
        })
    }
}

/// A validated source–receiver pair in a shoebox room.
#[derive(Debug, Clone)]
pub struct Simulation {
    room: ShoeboxRoom,
    source: Vec3,
    receiver: Vec3,
    medium: Medium,
}

impl Simulation {
    /// The underlying room.
    #[must_use]
    pub fn room(&self) -> &ShoeboxRoom {
        &self.room
    }

    /// Energy echogram `(arrival_time_s, energy_weight)` up to `max_order`
    /// reflections at `freq_hz`.
    ///
    /// # Errors
    /// Propagates the image-source errors.
    pub fn echogram(&self, freq_hz: f64, max_order: u32) -> Result<Vec<(f64, f64)>, AcousticError> {
        self.room
            .echogram(self.source, self.receiver, freq_hz, max_order, &self.medium)
    }

    /// Sampled impulse response at `sample_rate_hz`, built from the echogram
    /// with the energy-only convention of [`rir_from_echogram`].
    ///
    /// # Errors
    /// Propagates the echogram and RIR errors.
    pub fn impulse_response(
        &self,
        sample_rate_hz: f64,
        max_order: u32,
        seed: u64,
    ) -> Result<Vec<f64>, AcousticError> {
        let echogram = self.echogram(1000.0, max_order)?;
        rir_from_echogram(&echogram, sample_rate_hz, seed)
    }

    /// Sabine reverberation time (s) at `freq_hz`, from the room's mean
    /// absorption.
    #[must_use]
    pub fn sabine_rt60(&self, freq_hz: f64) -> f64 {
        self.room.sabine_rt60(freq_hz)
    }

    /// Distance from source to receiver in metres.
    #[must_use]
    pub fn direct_distance(&self) -> f64 {
        self.source.distance(self.receiver)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn room_builder() -> RoomBuilder {
        RoomBuilder::new()
            .dimensions([8.0, 5.0, 3.0])
            .source(Vec3::new(1.5, 1.2, 1.4))
            .receiver(Vec3::new(6.2, 3.7, 1.6))
    }

    #[test]
    fn builds_with_defaults() {
        let sim = room_builder().build().unwrap();
        // (4.7, 2.5, 0.2) m separation: √(4.7² + 2.5² + 0.2²) = √28.38.
        assert!((sim.direct_distance() - 28.38_f64.sqrt()).abs() < 1e-12);
        assert!(sim.sabine_rt60(1000.0) > 0.0);
    }

    #[test]
    fn direct_arrival_matches_the_geometry() {
        let sim = room_builder().build().unwrap();
        let echo = sim.echogram(1000.0, 0).unwrap();
        assert_eq!(echo.len(), 1);
        let expected = sim.direct_distance() / Medium::air_20c().speed_of_sound;
        assert!((echo[0].0 - expected).abs() < 1e-12);
    }

    #[test]
    fn missing_fields_are_reported() {
        assert!(RoomBuilder::new()
            .source(Vec3::new(1.0, 1.0, 1.0))
            .receiver(Vec3::new(2.0, 2.0, 2.0))
            .build()
            .is_err());
        assert!(RoomBuilder::new()
            .dimensions([8.0, 5.0, 3.0])
            .receiver(Vec3::new(2.0, 2.0, 2.0))
            .build()
            .is_err());
    }

    #[test]
    fn source_and_receiver_must_be_inside() {
        let outside = room_builder().source(Vec3::new(9.0, 1.0, 1.0)).build();
        assert!(matches!(outside, Err(AcousticError::InvalidGeometry(_))));
        let outside = room_builder().receiver(Vec3::new(1.0, 6.0, 1.0)).build();
        assert!(matches!(outside, Err(AcousticError::InvalidGeometry(_))));
    }

    #[test]
    fn impulse_response_has_the_expected_length_and_energy() {
        let sim = room_builder().build().unwrap();
        let rir = sim.impulse_response(48_000.0, 6, 1).unwrap();
        assert!(!rir.is_empty());
        let energy: f64 = rir.iter().map(|h| h * h).sum();
        let echo: Vec<(f64, f64)> = sim.echogram(1000.0, 6).unwrap();
        let weights: f64 = echo.iter().map(|&(_, w)| w).sum();
        // Coincident arrivals can add coherently, so allow a loose bound.
        assert!(energy > 0.0 && energy <= 1.5 * weights + 1e-12);
    }

    #[test]
    fn echogram_contains_no_negative_weights() {
        let sim = room_builder().build().unwrap();
        let echo: Vec<(f64, f64)> = sim.echogram(500.0, 4).unwrap();
        assert!(echo.iter().all(|&(t, w)| t >= 0.0 && w >= 0.0));
    }
}
