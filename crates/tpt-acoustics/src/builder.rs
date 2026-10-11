//! A small builder over the shoebox workflow.
//!
//! [`RoomBuilder`](crate::builder::RoomBuilder) collects a box-shaped room, a
//! wall material, a source and a receiver, and validates them once.
//! [`Simulation`](crate::builder::Simulation) then answers the common
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

    /// Ray-traced T30 (seconds) at 1 kHz from `n_rays` particles, or `None`
    /// if the decay is too short to fit.
    ///
    /// # Errors
    /// Propagates the ray tracer's validation errors.
    pub fn ray_t30(&self, n_rays: u32, seed: u64) -> Result<Option<f64>, AcousticError> {
        use tpt_acou_ray::ism::SPEED_OF_SOUND_REFERENCE;
        use tpt_acou_ray::ray::{rt60_from_curve, schroeder_curve, trace_shoebox, RayTracerConfig};
        let config = RayTracerConfig {
            n_rays,
            time_limit_s: 3.0,
            capture_radius: 0.5,
            min_energy: 1e-9,
            seed,
            speed_of_sound: SPEED_OF_SOUND_REFERENCE,
        };
        let traced = trace_shoebox(&self.room, self.source, self.receiver, 1000.0, &config)?;
        let curve = schroeder_curve(&traced.detected_times, &traced.detected_energies, 0.005);
        Ok(rt60_from_curve(&curve, 0.005))
    }

    /// Speech Transmission Index (male weights) from the hybrid RIR: image
    /// sources before the crossover, ray tail after it. Falls back to the
    /// image-source response alone if the ray tracer captured no energy in the
    /// calibration window.
    ///
    /// # Errors
    /// Propagates the tracer, echogram and STI errors.
    pub fn sti(&self, n_rays: u32, seed: u64) -> Result<f64, AcousticError> {
        use tpt_acou_ray::ism::SPEED_OF_SOUND_REFERENCE;
        use tpt_acou_ray::ray::{trace_shoebox, RayTracerConfig};
        use tpt_acou_ray::rir::hybrid_rir;
        use tpt_acou_ray::sti::{sti_from_rir, STI_WEIGHTS_MALE};
        let fs = 48_000.0;
        let echogram = self.echogram(1000.0, 12)?;
        let config = RayTracerConfig {
            n_rays,
            time_limit_s: 3.0,
            capture_radius: 0.5,
            min_energy: 1e-9,
            seed,
            speed_of_sound: SPEED_OF_SOUND_REFERENCE,
        };
        let traced = trace_shoebox(&self.room, self.source, self.receiver, 1000.0, &config)?;
        let samples = match hybrid_rir(
            &echogram,
            &traced.detected_times,
            &traced.detected_energies,
            fs,
            0.08,
            0.01,
            seed,
        ) {
            Ok(h) => h.samples,
            Err(_) => rir_from_echogram(&echogram, fs, seed)?,
        };
        sti_from_rir(&samples, fs, &STI_WEIGHTS_MALE)
    }

    /// Sound strength G in dB from the image-source echogram (12 orders).
    ///
    /// # Errors
    /// Propagates the echogram errors. Returns `None` inside the `Ok` if the
    /// echogram carries no energy.
    pub fn strength_db(&self) -> Result<Option<f64>, AcousticError> {
        let echogram = self.echogram(1000.0, 12)?;
        Ok(tpt_acou_ray::iso3382::strength_db(&echogram))
    }

    /// Distance from source to receiver in metres.
    #[must_use]
    pub fn direct_distance(&self) -> f64 {
        self.source.distance(self.receiver)
    }

    /// Wave-based impulse response: modal solution below the Schroeder
    /// frequency, image sources above it, signed and phase-consistent. This is
    /// the room's modal wave solve, not the energy-only echogram.
    ///
    /// # Errors
    /// Propagates the broadband solver's sampling and geometry errors.
    pub fn broadband_ir(
        &self,
        sample_rate_hz: f64,
        n_samples: usize,
        max_order: u32,
    ) -> Result<Vec<f64>, AcousticError> {
        tpt_acou_ray::broadband::hybrid_broadband_ir(
            &self.room,
            self.source,
            self.receiver,
            sample_rate_hz,
            n_samples,
            max_order,
            &self.medium,
        )
    }

    /// Three-way impulse response: modal low frequencies and image-source high
    /// frequencies up to `crossover_s`, then a ray-traced energy tail. The
    /// rays are traced at 1 kHz and the tail is calibrated to the broadband
    /// energy near the crossover (see
    /// [`tpt_acou_ray::broadband::three_way_ir`]).
    ///
    /// # Errors
    /// Propagates the tracer, broadband and calibration errors.
    pub fn three_way_impulse_response(
        &self,
        sample_rate_hz: f64,
        n_samples: usize,
        max_order: u32,
        n_rays: u32,
        crossover_s: f64,
        seed: u64,
    ) -> Result<Vec<f64>, AcousticError> {
        use tpt_acou_ray::broadband::three_way_ir;
        use tpt_acou_ray::ism::SPEED_OF_SOUND_REFERENCE;
        use tpt_acou_ray::ray::{trace_shoebox, RayTracerConfig};
        let config = RayTracerConfig {
            n_rays,
            time_limit_s: 3.0,
            capture_radius: 0.5,
            min_energy: 1e-9,
            seed,
            speed_of_sound: SPEED_OF_SOUND_REFERENCE,
        };
        let traced = trace_shoebox(&self.room, self.source, self.receiver, 1000.0, &config)?;
        let out = three_way_ir(
            &self.room,
            self.source,
            self.receiver,
            sample_rate_hz,
            n_samples,
            max_order,
            &self.medium,
            &traced.detected_times,
            &traced.detected_energies,
            crossover_s,
            0.01,
            seed,
        )?;
        Ok(out.samples)
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
    fn ray_t30_sti_and_strength_are_available_from_the_simulation() {
        let sim = room_builder().build().unwrap();
        let t30 = sim.ray_t30(5_000, 3).unwrap();
        if let Some(t) = t30 {
            assert!(t > 0.0 && t.is_finite(), "{t}");
        }
        let sti = sim.sti(5_000, 3).unwrap();
        assert!((0.0..=1.0).contains(&sti), "{sti}");
        let g = sim.strength_db().unwrap().expect("energy present");
        assert!(g.is_finite());
    }

    #[test]
    fn echogram_contains_no_negative_weights() {
        let sim = room_builder().build().unwrap();
        let echo: Vec<(f64, f64)> = sim.echogram(500.0, 4).unwrap();
        assert!(echo.iter().all(|&(t, w)| t >= 0.0 && w >= 0.0));
    }

    #[test]
    fn broadband_and_three_way_responses_have_the_requested_length() {
        let sim = room_builder().build().unwrap();
        let wave = sim.broadband_ir(4_000.0, 2_048, 3).unwrap();
        assert_eq!(wave.len(), 2_048);
        assert!(wave.iter().all(|v| v.is_finite()));
        assert!(wave.iter().any(|v| v.abs() > 0.0));

        let three = sim
            .three_way_impulse_response(4_000.0, 2_048, 3, 20_000, 0.1, 5)
            .unwrap();
        assert_eq!(three.len(), 2_048);
        assert!(three.iter().all(|v| v.is_finite()));
    }
}
