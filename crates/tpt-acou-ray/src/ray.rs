//! Monte Carlo energy ray tracing for late reverberation.
//!
//! Energy particles are shot from the source; each bounce deposits
//! `alpha·E` on the wall (absorbed), redirects diffusely with probability
//! `s(f)` (cosine-weighted Lambertian) or specularly otherwise, and is
//! detected each time its straight segment passes within the receiver's
//! capture radius; detection does not stop the particle. The tracer keeps an
//! explicit energy *partition* `injected = absorbed + remaining`, which is
//! what the energy-conservation property test (and `tpt-acou-verify`) checks.
//!
//! Decay curves use Schroeder backward integration; `T30` is extrapolated
//! from the −5 dB…−35 dB window of the fitted decay slope.

use alloc::vec;
use alloc::vec::Vec;

use tpt_acou_core::error::AcousticError;

use crate::ism::ShoeboxRoom;
use crate::rng::Xoshiro256;
use crate::scene::Scene;
use crate::vec3::Vec3;

/// Tracer configuration. Energies are normalized so the total injected
/// energy equals 1.
#[derive(Debug, Clone)]
pub struct RayTracerConfig {
    /// Number of energy particles.
    pub n_rays: u32,
    /// Time limit per particle in seconds.
    pub time_limit_s: f64,
    /// Receiver capture radius in metres (spherical detector).
    pub capture_radius: f64,
    /// Particles below this energy are terminated (their energy goes to
    /// `remaining`).
    pub min_energy: f64,
    /// RNG seed (deterministic across platforms).
    pub seed: u64,
    /// Speed of sound in m/s.
    pub speed_of_sound: f64,
}

impl Default for RayTracerConfig {
    fn default() -> Self {
        Self {
            n_rays: 10_000,
            time_limit_s: 2.0,
            capture_radius: 0.5,
            min_energy: 1e-7,
            seed: 0x5EED,
            speed_of_sound: crate::ism::SPEED_OF_SOUND_REFERENCE,
        }
    }
}

/// Result of a Monte Carlo run: the energy partition plus detected arrivals.
#[derive(Debug, Clone, Default)]
pub struct RayTraceResult {
    /// Arrival times of captures (s, unsorted). A particle that passes the
    /// receiver several times is recorded once per pass.
    pub detected_times: Vec<f64>,
    /// Energy carried by each capture (same length as times). This is a
    /// detection tally, not part of the energy partition: particles keep
    /// their energy after capture.
    pub detected_energies: Vec<f64>,
    /// Total energy injected (1.0 for a standard run).
    pub injected: f64,
    /// Energy absorbed at boundaries.
    pub absorbed: f64,
    /// Energy still in flight (or below the energy floor / past the time
    /// limit) when tracing ended.
    pub remaining: f64,
    /// Mean number of bounces per particle.
    pub mean_bounces: f64,
}

impl RayTraceResult {
    /// The energy accounting identity `injected = absorbed + remaining` as an
    /// absolute error (0 for exact bookkeeping).
    #[must_use]
    pub fn energy_partition_error(&self) -> f64 {
        (self.injected - (self.absorbed + self.remaining)).abs()
    }
}

/// Closest distance from `p` to segment `a→b`, plus the arc-length position
/// of the closest point along the segment (metres).
#[must_use]
pub fn point_segment_closest(p: Vec3, a: Vec3, b: Vec3) -> (f64, f64) {
    let ab = b - a;
    let len2 = ab.norm_sqr();
    let t = if len2 < 1e-30 {
        0.0
    } else {
        ((p - a).dot(ab) / len2).clamp(0.0, 1.0)
    };
    let closest = a + ab * t;
    (p.distance(closest), t * libm::sqrt(len2))
}

/// True if segment `a→b` passes within `radius` of `p`.
#[must_use]
pub fn segment_near_point(a: Vec3, b: Vec3, p: Vec3, radius: f64) -> bool {
    point_segment_closest(p, a, b).0 <= radius
}

/// Exit intersection of a ray inside the box `[0,lx]×[0,ly]×[0,lz]`:
/// `(distance, inward_unit_normal, wall_index)` with walls numbered as in
/// [`crate::ism`] (`0..6`).
#[must_use]
pub fn box_exit(pos: Vec3, dir: Vec3, lx: f64, ly: f64, lz: f64) -> Option<(f64, Vec3, usize)> {
    let bounds = [lx, ly, lz];
    let p = [pos.x, pos.y, pos.z];
    let d = [dir.x, dir.y, dir.z];
    let mut best_t = f64::INFINITY;
    let mut best_wall = usize::MAX;
    for axis in 0..3 {
        if d[axis].abs() < 1e-15 {
            continue;
        }
        let t = if d[axis] > 0.0 {
            (bounds[axis] - p[axis]) / d[axis]
        } else {
            -p[axis] / d[axis]
        };
        if t > 1e-12 && t < best_t {
            best_t = t;
            best_wall = axis * 2 + usize::from(d[axis] > 0.0);
        }
    }
    let normal = match best_wall {
        0 => Vec3::new(1.0, 0.0, 0.0),
        1 => Vec3::new(-1.0, 0.0, 0.0),
        2 => Vec3::new(0.0, 1.0, 0.0),
        3 => Vec3::new(0.0, -1.0, 0.0),
        4 => Vec3::new(0.0, 0.0, 1.0),
        5 => Vec3::new(0.0, 0.0, -1.0),
        _ => return None,
    };
    Some((best_t, normal, best_wall))
}

/// Cosine-weighted hemisphere sample about `n` (the `normalize(n + u)` trick:
/// adding a uniform sphere vector to the normal is cosine-distributed).
fn cosine_hemisphere(rng: &mut Xoshiro256, n: Vec3) -> Vec3 {
    for _ in 0..8 {
        let u = rng.next_unit_direction();
        let d = (n + u).normalized();
        if d.norm_sqr() > 0.5 {
            return d;
        }
    }
    n
}

impl RayTracerConfig {
    /// Checks the run parameters.
    ///
    /// # Errors
    /// Returns [`AcousticError::InvalidArgument`] for `n_rays == 0`, a
    /// non-positive or non-finite time limit, speed of sound or capture
    /// radius, or a negative or non-finite energy floor.
    pub fn validate(&self) -> Result<(), AcousticError> {
        if self.n_rays == 0 {
            return Err(AcousticError::InvalidArgument("n_rays must be > 0"));
        }
        if !(self.time_limit_s.is_finite() && self.time_limit_s > 0.0) {
            return Err(AcousticError::InvalidArgument(
                "time_limit_s must be finite and > 0",
            ));
        }
        if !(self.speed_of_sound.is_finite() && self.speed_of_sound > 0.0) {
            return Err(AcousticError::InvalidArgument(
                "speed_of_sound must be finite and > 0",
            ));
        }
        if !(self.capture_radius.is_finite() && self.capture_radius >= 0.0) {
            return Err(AcousticError::InvalidArgument(
                "capture_radius must be finite and >= 0",
            ));
        }
        if !(self.min_energy.is_finite() && self.min_energy >= 0.0) {
            return Err(AcousticError::InvalidArgument(
                "min_energy must be finite and >= 0",
            ));
        }
        Ok(())
    }
}

/// Per-ray tally, combined into a [`RayTraceResult`].
#[derive(Debug, Clone, Default)]
struct RayTally {
    /// `(arrival_time_s, energy)` for every receiver pass.
    detections: Vec<(f64, f64)>,
    absorbed: f64,
    remaining: f64,
    bounces: u64,
}

/// Deterministic per-ray seed: splitmix64 of the run seed and the ray index.
/// Each ray has its own stream, so the result does not depend on the order
/// rays are traced in, and the parallel and sequential paths agree exactly.
#[must_use]
pub fn ray_seed(seed: u64, index: u32) -> u64 {
    let mut z = seed ^ (u64::from(index) + 1).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    z ^ (z >> 31)
}

/// Runs `trace` for every ray index and returns the tallies in index order.
#[cfg(feature = "parallel")]
fn run_rays<F>(n_rays: u32, trace: F) -> Vec<RayTally>
where
    F: Fn(u32) -> RayTally + Sync + Send,
{
    use rayon::prelude::*;
    (0..n_rays).into_par_iter().map(trace).collect()
}

/// Runs `trace` for every ray index and returns the tallies in index order.
#[cfg(not(feature = "parallel"))]
fn run_rays<F>(n_rays: u32, trace: F) -> Vec<RayTally>
where
    F: Fn(u32) -> RayTally,
{
    (0..n_rays).map(trace).collect()
}

/// Sums per-ray tallies into the run result.
fn combine(tallies: Vec<RayTally>, n_rays: u32) -> RayTraceResult {
    let mut out = RayTraceResult {
        injected: 1.0,
        ..RayTraceResult::default()
    };
    let mut bounces = 0_u64;
    for t in tallies {
        for (time, energy) in t.detections {
            out.detected_times.push(time);
            out.detected_energies.push(energy);
        }
        out.absorbed += t.absorbed;
        out.remaining += t.remaining;
        bounces += t.bounces;
    }
    out.mean_bounces = bounces as f64 / f64::from(n_rays);
    out
}

/// Traces energy particles through a shoebox room (analytic wall
/// intersection).
///
/// Every ray has its own seeded stream (see [`ray_seed`]), so the result is
/// deterministic for a given seed. With the `parallel` feature the rays are
/// traced on the rayon pool, with identical output.
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] if the configuration fails
/// [`RayTracerConfig::validate`], or if `src` or `rcv` lies outside the room.
pub fn trace_shoebox(
    room: &ShoeboxRoom,
    src: Vec3,
    rcv: Vec3,
    freq_hz: f64,
    config: &RayTracerConfig,
) -> Result<RayTraceResult, AcousticError> {
    config.validate()?;
    if !room.contains(src) {
        return Err(AcousticError::InvalidArgument(
            "source lies outside the room",
        ));
    }
    if !room.contains(rcv) {
        return Err(AcousticError::InvalidArgument(
            "receiver lies outside the room",
        ));
    }
    let [lx, ly, lz] = room.dims();
    let alpha: Vec<f64> = (0..6)
        .map(|w| room.wall(w).absorption_at(freq_hz))
        .collect();
    let scatter: Vec<f64> = (0..6)
        .map(|w| room.wall(w).scattering_at(freq_hz))
        .collect();
    let e0 = 1.0 / f64::from(config.n_rays);

    let trace = |index: u32| -> RayTally {
        let mut rng = Xoshiro256::new(ray_seed(config.seed, index));
        let mut tally = RayTally::default();
        let mut p = src;
        let mut dir = rng.next_unit_direction();
        let mut e = e0;
        let mut t = 0.0_f64;
        loop {
            if t > config.time_limit_s || e < config.min_energy {
                tally.remaining += e;
                break;
            }
            let Some((t_hit, normal, wall)) = box_exit(p, dir, lx, ly, lz) else {
                tally.remaining += e;
                break;
            };
            let from = p;
            let to = from + dir * t_hit;
            if segment_near_point(from, to, rcv, config.capture_radius) {
                let (_, dist_along) = point_segment_closest(rcv, from, to);
                tally
                    .detections
                    .push((t + dist_along / config.speed_of_sound, e));
                // The particle keeps tracing: terminating at first capture
                // would truncate the late decay and bias T30 steep.
            }
            tally.absorbed += e * alpha[wall];
            e *= 1.0 - alpha[wall];
            t += t_hit / config.speed_of_sound;
            tally.bounces += 1;
            dir = if rng.next_f64() < scatter[wall] {
                cosine_hemisphere(&mut rng, normal)
            } else {
                dir.reflect(normal)
            };
            p = to + normal * 1e-9;
        }
        tally
    };
    Ok(combine(run_rays(config.n_rays, trace), config.n_rays))
}

/// Traces energy particles through an arbitrary (closed) triangle scene,
/// intersecting via the BVH and reading coefficients from the hit material.
///
/// For a ray that escapes the scene (open mesh), leftover energy is counted
/// as `remaining`. Seeding and parallelism are as for [`trace_shoebox`].
///
/// # Errors
/// Returns [`AcousticError::InvalidArgument`] if the configuration fails
/// [`RayTracerConfig::validate`].
pub fn trace_scene(
    scene: &Scene,
    src: Vec3,
    rcv: Vec3,
    freq_hz: f64,
    config: &RayTracerConfig,
) -> Result<RayTraceResult, AcousticError> {
    config.validate()?;
    let e0 = 1.0 / f64::from(config.n_rays);

    let trace = |index: u32| -> RayTally {
        let mut rng = Xoshiro256::new(ray_seed(config.seed, index));
        let mut tally = RayTally::default();
        let mut p = src;
        let mut dir = rng.next_unit_direction();
        let mut e = e0;
        let mut t = 0.0_f64;
        loop {
            if t > config.time_limit_s || e < config.min_energy {
                tally.remaining += e;
                break;
            }
            let Some(hit) = scene.hit(p, dir, f64::INFINITY) else {
                tally.remaining += e;
                break;
            };
            let from = p;
            let to = hit.point;
            if segment_near_point(from, to, rcv, config.capture_radius) {
                let (_, dist_along) = point_segment_closest(rcv, from, to);
                tally
                    .detections
                    .push((t + dist_along / config.speed_of_sound, e));
            }
            let material = &scene.materials()[hit.material];
            let a = material.absorption_at(freq_hz);
            let s = material.scattering_at(freq_hz);
            tally.absorbed += e * a;
            e *= 1.0 - a;
            t += hit.t / config.speed_of_sound;
            tally.bounces += 1;
            dir = if rng.next_f64() < s {
                cosine_hemisphere(&mut rng, hit.normal)
            } else {
                dir.reflect(hit.normal)
            };
            p = to + hit.normal * 1e-9;
        }
        tally
    };
    Ok(combine(run_rays(config.n_rays, trace), config.n_rays))
}

/// Backward-integrated (Schroeder) decay curve sampled at `bin_s` spacing:
/// `curve[i] = Σ_{j ≥ i} energy_in_bin[j]`.
#[must_use]
pub fn schroeder_curve(times: &[f64], energies: &[f64], bin_s: f64) -> Vec<f64> {
    if times.is_empty() {
        return vec![];
    }
    let t_max = times.iter().copied().fold(0.0_f64, f64::max);
    let n_bins = libm::ceil(t_max / bin_s) as usize + 1;
    let mut bins = vec![0.0; n_bins];
    for (i, &t) in times.iter().enumerate() {
        let b = ((t / bin_s) as usize).min(n_bins - 1);
        bins[b] += energies[i];
    }
    let mut out = vec![0.0; n_bins];
    let mut acc = 0.0;
    for i in (0..n_bins).rev() {
        acc += bins[i];
        out[i] = acc;
    }
    out
}

/// `T30` reverberation time from a Schroeder curve: least-squares decay-slope
/// fit over the −5 dB…−35 dB window, extrapolated to the −60 dB span.
/// Returns `None` when the curve lacks enough decay data.
#[must_use]
pub fn rt60_from_curve(curve: &[f64], bin_s: f64) -> Option<f64> {
    decay_time_from_curve(curve, bin_s, -5.0, -35.0)
}

/// Decay time from a Schroeder curve, fitted between `upper_db` and
/// `lower_db` (both negative, `upper_db > lower_db`, in dB relative to the
/// curve's first sample) and extrapolated to a 60 dB decay. T30 uses
/// `(-5, -35)`, T20 `(-5, -25)`, EDT `(0, -10)`.
///
/// Returns `None` if the curve does not reach `lower_db` or the window is
/// too short to fit.
#[must_use]
pub fn decay_time_from_curve(
    curve: &[f64],
    bin_s: f64,
    upper_db: f64,
    lower_db: f64,
) -> Option<f64> {
    if curve.len() < 4 || curve[0] <= 0.0 || upper_db <= lower_db {
        return None;
    }
    let peak = curve[0];
    let db = |v: f64| 10.0 * libm::log10((v / peak).max(f64::MIN_POSITIVE));
    let mut i5 = None;
    let mut i35 = None;
    for (i, &v) in curve.iter().enumerate() {
        let l = db(v);
        if i5.is_none() && l <= upper_db {
            i5 = Some(i);
        }
        if l <= lower_db {
            i35 = Some(i);
            break;
        }
    }
    let i5 = i5?;
    let i35 = i35?;
    if i35 <= i5 + 1 {
        return None;
    }
    // Least squares l(t) = m·t + b over the window.
    let (mut sx, mut sy, mut sxx, mut sxy, mut n) = (0.0, 0.0, 0.0, 0.0, 0.0_f64);
    for i in i5..=i35 {
        let t = i as f64 * bin_s;
        let y = db(curve[i]);
        sx += t;
        sy += y;
        sxx += t * t;
        sxy += t * y;
        n += 1.0;
    }
    let denom = n * sxx - sx * sx;
    if denom.abs() < 1e-30 {
        return None;
    }
    let m = (n * sxy - sx * sy) / denom; // dB/s, negative for a decay
    if m >= -1e-9 {
        return None;
    }
    // Extrapolate: the fitted line spans 60 dB in 60/|m| seconds.
    Some(60.0 / -m)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ism::SPEED_OF_SOUND_REFERENCE;
    use crate::material::BoundaryMaterial;

    fn uniform_room(dims: [f64; 3], alpha: f64) -> ShoeboxRoom {
        let mut table = alloc::vec![(125.0, alpha), (4000.0, alpha)];
        table.dedup();
        ShoeboxRoom::uniform(
            dims,
            BoundaryMaterial::new("test", table, alloc::vec![(125.0, 0.0), (4000.0, 0.0)]).unwrap(),
        )
        .unwrap()
    }

    #[test]
    fn box_exit_finds_the_right_wall() {
        // Straight down +z from the center of a 2×2×2 box: exit at +z wall.
        let (t, n, w) = box_exit(
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(0.0, 0.0, 1.0),
            2.0,
            2.0,
            2.0,
        )
        .unwrap();
        assert!((t - 1.0).abs() < 1e-12);
        assert_eq!(w, 5);
        assert!((n.z + 1.0).abs() < 1e-12);
        // Straight -x: wall 0.
        let (_, _, w) = box_exit(
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(-1.0, 0.0, 0.0),
            2.0,
            2.0,
            2.0,
        )
        .unwrap();
        assert_eq!(w, 0);
        // Degenerate direction: no exit.
        assert!(box_exit(
            Vec3::new(1.0, 1.0, 1.0),
            Vec3::new(0.0, 0.0, 0.0),
            2.0,
            2.0,
            2.0
        )
        .is_none());
    }

    #[test]
    fn segment_detection_and_closest_point() {
        let a = Vec3::new(0.0, 0.0, 0.0);
        let b = Vec3::new(10.0, 0.0, 0.0);
        let (d, s) = point_segment_closest(Vec3::new(5.0, 1.0, 0.0), a, b);
        assert!((d - 1.0).abs() < 1e-12);
        assert!((s - 5.0).abs() < 1e-12);
        // Beyond the end: clamped.
        let (d, s) = point_segment_closest(Vec3::new(15.0, 0.0, 0.0), a, b);
        assert!((d - 5.0).abs() < 1e-12 && (s - 10.0).abs() < 1e-12);
        assert!(segment_near_point(a, b, Vec3::new(5.0, 0.4, 0.0), 0.5));
        assert!(!segment_near_point(a, b, Vec3::new(5.0, 0.6, 0.0), 0.5));
    }

    #[test]
    fn energy_partition_is_exact_in_lossless_rooms() {
        let room = uniform_room([8.0, 6.0, 3.0], 0.0);
        let config = RayTracerConfig {
            n_rays: 2_000,
            time_limit_s: 0.5,
            capture_radius: 0.4,
            min_energy: 1e-9,
            seed: 1,
            speed_of_sound: SPEED_OF_SOUND_REFERENCE,
        };
        let res = trace_shoebox(
            &room,
            Vec3::new(2.0, 2.0, 1.5),
            Vec3::new(6.0, 4.0, 1.5),
            1000.0,
            &config,
        )
        .unwrap();
        assert_eq!(res.absorbed, 0.0);
        assert!(
            res.energy_partition_error() < 1e-9,
            "partition error {}",
            res.energy_partition_error()
        );
        // Particles keep tracing after capture, so the receiver saw at least
        // one pass in a room this size.
        assert!(!res.detected_times.is_empty());
    }

    #[test]
    fn energy_partition_holds_with_absorption_and_scattering() {
        let room =
            ShoeboxRoom::uniform([7.0, 5.0, 3.0], BoundaryMaterial::absorber_panel()).unwrap();
        let config = RayTracerConfig {
            n_rays: 2_000,
            time_limit_s: 1.0,
            capture_radius: 0.3,
            min_energy: 1e-8,
            seed: 2,
            speed_of_sound: SPEED_OF_SOUND_REFERENCE,
        };
        let res = trace_shoebox(
            &room,
            Vec3::new(1.0, 1.0, 1.5),
            Vec3::new(6.0, 4.0, 1.5),
            1000.0,
            &config,
        )
        .unwrap();
        assert!(res.absorbed > 0.0);
        assert!(
            res.energy_partition_error() < 1e-9,
            "partition error {}",
            res.energy_partition_error()
        );
    }

    #[test]
    fn rt60_tracks_sabine_within_tolerance() {
        // 10×8×6 m room, alpha = 0.15: Sabine T60 = 0.161·480/(0.15·376) ≈ 1.37 s.
        let room = uniform_room([10.0, 8.0, 6.0], 0.15);
        let sabine = room.sabine_rt60(1000.0);
        let config = RayTracerConfig {
            n_rays: 30_000,
            time_limit_s: 3.0,
            capture_radius: 0.5,
            min_energy: 1e-9,
            seed: 7,
            speed_of_sound: SPEED_OF_SOUND_REFERENCE,
        };
        let res = trace_shoebox(
            &room,
            Vec3::new(2.0, 2.0, 3.0),
            Vec3::new(8.0, 6.0, 3.0),
            1000.0,
            &config,
        )
        .unwrap();
        let curve = schroeder_curve(&res.detected_times, &res.detected_energies, 0.005);
        let t30 = rt60_from_curve(&curve, 0.005).expect("decay long enough for T30");
        let rel = (t30 - sabine).abs() / sabine;
        assert!(
            rel < 0.10,
            "t30 {t30:.3} vs sabine {sabine:.3} (rel {rel:.3})"
        );
    }

    #[test]
    fn rt60_decreases_with_absorption() {
        let mut t60s = [0.0; 2];
        for (i, alpha) in [0.1, 0.4].iter().enumerate() {
            let room = uniform_room([9.0, 7.0, 5.0], *alpha);
            let config = RayTracerConfig {
                n_rays: 20_000,
                time_limit_s: 3.0,
                capture_radius: 0.5,
                min_energy: 1e-9,
                seed: 5,
                speed_of_sound: SPEED_OF_SOUND_REFERENCE,
            };
            let res = trace_shoebox(
                &room,
                Vec3::new(2.0, 2.0, 2.5),
                Vec3::new(7.0, 5.0, 2.5),
                1000.0,
                &config,
            )
            .unwrap();
            let curve = schroeder_curve(&res.detected_times, &res.detected_energies, 0.005);
            t60s[i] = rt60_from_curve(&curve, 0.005).expect("t30");
        }
        assert!(
            t60s[1] < t60s[0],
            "expected decay: {} vs {}",
            t60s[1],
            t60s[0]
        );
    }

    #[test]
    fn scene_tracer_matches_shoebox_partition() {
        let room =
            ShoeboxRoom::uniform([6.0, 5.0, 3.0], BoundaryMaterial::concrete_painted()).unwrap();
        let scene = Scene::from_shoebox(&room).unwrap();
        let config = RayTracerConfig {
            n_rays: 2_000,
            time_limit_s: 1.0,
            capture_radius: 0.3,
            min_energy: 1e-8,
            seed: 3,
            speed_of_sound: SPEED_OF_SOUND_REFERENCE,
        };
        let src = Vec3::new(1.5, 1.5, 1.5);
        let rcv = Vec3::new(4.5, 3.5, 1.5);
        let a = trace_shoebox(&room, src, rcv, 1000.0, &config).unwrap();
        let b = trace_scene(&scene, src, rcv, 1000.0, &config).unwrap();
        // Deterministic seeds make the two paths produce the *same* random
        // sequence; intersection differs only by epsilon nudges.
        assert!(b.energy_partition_error() < 1e-9);
        assert!(b.absorbed > 0.0);
        assert_eq!(
            !a.detected_energies.is_empty(),
            !b.detected_energies.is_empty()
        );
    }

    #[test]
    fn rt60_curve_rejects_short_data() {
        assert!(rt60_from_curve(&[], 0.005).is_none());
        assert!(rt60_from_curve(&[1.0, 0.999, 0.998], 0.005).is_none());
    }
}

#[cfg(test)]
mod determinism_tests {
    use super::*;
    use crate::material::BoundaryMaterial;

    #[test]
    fn same_seed_gives_identical_tallies() {
        let room =
            ShoeboxRoom::uniform([6.0, 5.0, 3.0], BoundaryMaterial::concrete_painted()).unwrap();
        let config = RayTracerConfig {
            n_rays: 800,
            time_limit_s: 0.5,
            capture_radius: 0.4,
            min_energy: 1e-9,
            seed: 123,
            speed_of_sound: crate::ism::SPEED_OF_SOUND_REFERENCE,
        };
        let src = Vec3::new(1.0, 1.0, 1.5);
        let rcv = Vec3::new(4.5, 3.5, 1.5);
        let a = trace_shoebox(&room, src, rcv, 1000.0, &config).unwrap();
        let b = trace_shoebox(&room, src, rcv, 1000.0, &config).unwrap();
        assert_eq!(a.detected_times, b.detected_times);
        assert_eq!(a.detected_energies, b.detected_energies);
        assert_eq!(a.absorbed, b.absorbed);
    }

    #[test]
    fn ray_seeds_differ_per_ray_and_per_run() {
        assert_ne!(ray_seed(1, 0), ray_seed(1, 1));
        assert_ne!(ray_seed(1, 0), ray_seed(2, 0));
    }
}
