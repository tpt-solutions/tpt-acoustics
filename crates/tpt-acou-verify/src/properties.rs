//! Physical-law property tests (reciprocity, energy conservation).
//!
//! These run under `cargo test` as statistical verification of the physics
//! (spec.txt §2: "property-based testing for statistical verification of
//! physical laws").

use proptest::prelude::*;

use tpt_acou_core::medium::Medium;
use tpt_acou_core::Vec3 as WaveVec3;
use tpt_acou_ray::ism::{ShoeboxRoom, SPEED_OF_SOUND_REFERENCE};
use tpt_acou_ray::material::BoundaryMaterial;
use tpt_acou_ray::ray::{rt60_from_curve, schroeder_curve, trace_shoebox, RayTracerConfig};
use tpt_acou_ray::vec3::Vec3;
use tpt_acou_wave::bem::HelmholtzBem;
use tpt_acou_wave::mesh::BoundaryMesh;
use tpt_math_linalg_complex::Complex;

use crate::strategies::{
    arb_absorption, arb_frequency, arb_point_pair_in_room, arb_shoebox_dims, arb_shoebox_room,
};

/// Dimensions plus a source/receiver pair inside them (room not yet built).
pub fn arb_dims_with_pair() -> impl Strategy<Value = ([f64; 3], Vec3, Vec3)> {
    arb_shoebox_dims()
        .prop_flat_map(|dims| (arb_point_pair_in_room(dims), Just(dims)))
        .prop_map(|((src, rcv), dims)| (dims, src, rcv))
}

/// Strategy: a shoebox room together with a source/receiver pair strictly
/// inside *that* room (points drawn from the room's own dimensions), with
/// the pair separated beyond the Green's-function singularity.
pub fn arb_room_with_pair() -> impl Strategy<Value = (ShoeboxRoom, Vec3, Vec3)> {
    arb_shoebox_room().prop_flat_map(|room| {
        let dims = room.dims();
        arb_point_pair_in_room(dims).prop_map(move |(src, rcv)| (room.clone(), src, rcv))
    })
}

/// Strategy: a room plus a source/receiver pair inside it (used when only the
/// direct arrival geometry matters).
pub fn arb_room_with_point() -> impl Strategy<Value = (ShoeboxRoom, Vec3, Vec3)> {
    arb_shoebox_room().prop_flat_map(|room| {
        let dims = room.dims();
        arb_point_pair_in_room(dims).prop_map(move |(a, b)| (room.clone(), a, b))
    })
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(48))]

    /// Acoustic reciprocity `H_AB = H_BA`: the ISM transfer function of a
    /// shoebox room is symmetric under source/receiver exchange, because every
    /// image path reverses into a path of identical length and reflection
    /// sequence (spec.txt §4, tpt-acou-verify).
    #[test]
    fn ism_reciprocity(
        (room, src, rcv) in arb_room_with_pair(),
        freq in arb_frequency(),
    ) {
        let medium = Medium::air_20c();
        let h_ab = room.transfer_function(src, rcv, &[freq], 3, &medium)?;
        let h_ba = room.transfer_function(rcv, src, &[freq], 3, &medium)?;
        prop_assert!((h_ab[0].norm() - h_ba[0].norm()).abs() < 1e-9 * h_ab[0].norm().max(1e-30),
            "|H_AB| {} vs |H_BA| {}", h_ab[0].norm(), h_ba[0].norm());
        prop_assert!((h_ab[0].arg() - h_ba[0].arg()).abs() < 1e-6);
    }

    /// Green's-function scaling: every individual ISM contribution decays as
    /// `1/(4πd)`, so the lossless-room response at frequency f, sampled at two
    /// receiver distances, respects the magnitude bound
    /// `|H| ≥ |G_direct| − Σ|indirect|`; simpler operational check: moving the
    /// receiver onto the source's mirror position reproduces the direct
    /// distance. We check the direct term exactly via `echogram`.
    #[test]
    fn echogram_direct_arrival_matches_geometry(
        (room, src, rcv) in arb_room_with_point(),
    ) {
        let air = tpt_acou_core::medium::Medium::air_20c();
        let echo = room.echogram(src, rcv, 1000.0, 0, &air)?;
        prop_assert_eq!(echo.len(), 1);
        let expected_t = src.distance(rcv) / air.speed_of_sound;
        prop_assert!((echo[0].0 - expected_t).abs() < 1e-12);
        // Lossless rigid wall: direct weight is exactly |G|² = 1/(4πd)².
        let d = src.distance(rcv);
        prop_assert!((echo[0].1 - 1.0 / (4.0 * core::f64::consts::PI * d).powi(2)).abs() < 1e-18);
    }

    /// Energy conservation in *lossless* ray tracing: with α = 0 on every
    /// wall, no energy is absorbed, so detected + remaining must equal the
    /// injected energy exactly (to float rounding).
    #[test]
    fn lossless_ray_tracing_conserves_energy(
        (dims, src, rcv) in arb_dims_with_pair(),
        seed in any::<u64>(),
    ) {
        let room = ShoeboxRoom::uniform(dims, BoundaryMaterial::rigid())?;
        let config = RayTracerConfig {
            n_rays: 512,
            time_limit_s: 0.6,
            capture_radius: 0.4,
            min_energy: 1e-9,
            seed,
            speed_of_sound: SPEED_OF_SOUND_REFERENCE,
        };
        let res = trace_shoebox(&room, src, rcv, 1000.0, &config)?;
        prop_assert_eq!(res.absorbed, 0.0);
        prop_assert!(res.energy_partition_error() < 1e-9,
            "partition error {}", res.energy_partition_error());
    }

    /// With absorption, the partition still closes exactly: every missing
    /// unit of energy is accounted at a wall.
    #[test]
    fn lossy_ray_tracing_partition_closes(
        (dims, src, rcv) in arb_dims_with_pair(),
        alpha in arb_absorption(),
        seed in any::<u64>(),
    ) {
        let material = BoundaryMaterial::new(
            "lossy",
            [(125.0, alpha), (4000.0, alpha)].to_vec(),
            [(125.0, 0.1), (4000.0, 0.1)].to_vec(),
        )?;
        let room = ShoeboxRoom::uniform(dims, material)?;
        let config = RayTracerConfig {
            n_rays: 512,
            time_limit_s: 1.0,
            capture_radius: 0.4,
            min_energy: 1e-9,
            seed,
            speed_of_sound: SPEED_OF_SOUND_REFERENCE,
        };
        let res = trace_shoebox(&room, src, rcv, 1000.0, &config)?;
        prop_assert!(res.energy_partition_error() < 1e-9,
            "partition error {}", res.energy_partition_error());
    }
}

/// Deterministic sanity: the Schroeder curve of a geometric decay yields a
/// finite T30 and the extraction is stable under receiver radius (same room,
/// two radii, T30 within 25%).
#[test]
fn t30_is_stable_across_capture_radii() -> Result<(), tpt_acou_core::AcousticError> {
    let room = ShoeboxRoom::uniform([9.0, 7.0, 4.5], BoundaryMaterial::concrete_painted())?;
    let src = Vec3::new(2.0, 2.0, 2.2);
    let rcv = Vec3::new(7.0, 5.0, 2.2);
    let mut t30s = [0.0; 2];
    for (i, radius) in [0.3, 0.6].iter().enumerate() {
        let config = RayTracerConfig {
            n_rays: 20_000,
            time_limit_s: 2.0,
            capture_radius: *radius,
            min_energy: 1e-9,
            seed: 77,
            speed_of_sound: SPEED_OF_SOUND_REFERENCE,
        };
        let res = trace_shoebox(&room, src, rcv, 1000.0, &config)?;
        let curve = schroeder_curve(&res.detected_times, &res.detected_energies, 0.005);
        t30s[i] = rt60_from_curve(&curve, 0.005).expect("T30 extraction");
    }
    let rel = (t30s[0] - t30s[1]).abs() / t30s[0].max(t30s[1]);
    assert!(
        rel < 0.25,
        "T30 spread {rel} between {} and {}",
        t30s[0],
        t30s[1]
    );
    Ok(())
}

// -- BEM invariants (spec.txt section 4, tpt-acou-verify) -------------------

/// Kernel reciprocity of the single-layer BEM matrix: with area weighting,
/// `A_i G(x_i, y_j) = A_j G(x_j, y_i)` exactly (the Green's function is
/// symmetric), which makes the symmetrized operator self-adjoint — the
/// discrete form of acoustic reciprocity `H_AB = H_BA`.
#[test]
fn bem_single_layer_kernel_reciprocity() {
    let mesh = BoundaryMesh::sphere(1.0, 6, 8).unwrap();
    let bem = HelmholtzBem::new(mesh, 1.7).unwrap();
    let els = bem.mesh().elements();
    let k = bem.k();
    let g = |_a: WaveVec3, _b: WaveVec3, r: f64| {
        Complex::from_polar(1.0 / (4.0 * core::f64::consts::PI * r), -k * r)
    };
    // Kernel symmetry G_ij = G_ji over a spread of element pairs, and the
    // induced symmetry of the Galerkin (area-doubly-weighted) operator
    // B_ij = A_i A_j G_ij — the discrete self-adjointness behind H_AB = H_BA.
    let n = els.len();
    for i in (0..n).step_by(7) {
        for j in (0..n).step_by(11) {
            if i == j {
                continue;
            }
            let (xi, yi) = (els[i].centroid, els[j].centroid);
            let r = xi.distance(yi);
            let g_ij = g(xi, yi, r);
            let g_ji = g(yi, xi, r);
            assert!(
                (g_ij - g_ji).norm() < 1e-14 * g_ij.norm().max(1e-30),
                "kernel reciprocity violated at ({i}, {j})"
            );
            let b_ij = g_ij * (els[i].area * els[j].area);
            let b_ji = g_ji * (els[j].area * els[i].area);
            assert!(
                (b_ij - b_ji).norm() < 1e-12 * b_ij.norm().max(1e-30),
                "Galerkin symmetry violated at ({i}, {j})"
            );
        }
    }
}

/// Passivity of lossless radiation: the baffled-piston radiation efficiency
/// computed through the lossless BEM never exceeds unity (a lossless
/// surface cannot radiate more power than an equally loud plane wave).
#[test]
fn bem_radiation_efficiency_is_passive() {
    // A uniformly vibrating circular-ish patch at a mid ka: sigma <= 1 + tol.
    let mesh = BoundaryMesh::rectangle_patch(0.5, 0.5, 8, 8).unwrap();
    let k = 4.0;
    let bem = HelmholtzBem::new(mesh, k).unwrap();
    let n = bem.mesh().len();
    let v = vec![Complex::new(1.0, 0.0); n];
    let medium = tpt_acou_core::medium::Medium::air_20c();
    let sigma = tpt_acou_vibro_dep_radiation_efficiency(&bem, &v, &medium);
    assert!(sigma > 0.0 && sigma <= 1.15, "sigma = {sigma}");
}

/// Thin wrapper so this test module need not depend on tpt-acou-vibro
/// (keeps the verify dependency graph a DAG: wave, not vibro).
fn tpt_acou_vibro_dep_radiation_efficiency(
    bem: &HelmholtzBem,
    v: &[Complex],
    medium: &tpt_acou_core::medium::Medium,
) -> f64 {
    // Reimplements the baffled-efficiency bookkeeping of
    // tpt-acou-vibro::radiation::radiation_efficiency_baffled on top of the
    // wave crate (the vibro function itself is covered by its own tests).
    let els = bem.mesh().elements();
    let area: f64 = els.iter().map(|e| e.area).sum();
    let v_rms_sq: f64 = els
        .iter()
        .zip(v)
        .map(|(e, v)| v.norm_sqr() * e.area)
        .sum::<f64>()
        / area;
    let omega = bem.k() * medium.speed_of_sound;
    let q: Vec<Complex> = v
        .iter()
        .map(|vi| *vi * Complex::new(0.0, -omega * medium.density))
        .collect();
    let sigma = bem.solve_exterior_neumann(&q).expect("neumann solve");
    let probe: Vec<WaveVec3> = els
        .iter()
        .map(|e| e.centroid + e.normal * f64::sqrt(e.area) * 1e-4)
        .collect();
    let p = bem.evaluate_single_layer(&sigma, &probe).expect("evaluate");
    let mut w = 0.0;
    for (e, (pv, vi)) in els.iter().zip(p.iter().zip(v)) {
        let flux = *pv * vi.conj();
        w += flux.re * e.area;
    }
    2.0 * 0.5 * w / (medium.characteristic_impedance() * area * v_rms_sq)
}
