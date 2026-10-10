//! Reusable proptest strategies for acoustic inputs.

use proptest::prelude::*;

use tpt_acou_ray::ism::ShoeboxRoom;
use tpt_acou_ray::material::BoundaryMaterial;
use tpt_acou_ray::vec3::Vec3;

/// Float in `[min, max]` with a finite, reasonable distribution.
pub fn arb_f64_range(min: f64, max: f64) -> impl Strategy<Value = f64> {
    (min..=max).prop_filter("must be finite", |v: &f64| v.is_finite())
}

/// Shoebox room dimensions in a physically sensible range (2–30 m per axis).
pub fn arb_shoebox_dims() -> impl Strategy<Value = [f64; 3]> {
    (
        arb_f64_range(2.0, 30.0),
        arb_f64_range(2.0, 30.0),
        arb_f64_range(2.0, 30.0),
    )
        .prop_map(|(x, y, z)| [x, y, z])
}

/// A point strictly inside a room of the given dimensions (0.1 m margin).
pub fn arb_point_in_room(dims: [f64; 3]) -> impl Strategy<Value = Vec3> {
    (
        arb_f64_range(0.1, dims[0] - 0.1),
        arb_f64_range(0.1, dims[1] - 0.1),
        arb_f64_range(0.1, dims[2] - 0.1),
    )
        .prop_map(|(x, y, z)| Vec3::new(x, y, z))
}

/// Source/receiver pair strictly inside a room of the given dimensions,
/// separated by at least 0.2 m (avoids the Green's-function singularity).
pub fn arb_point_pair_in_room(dims: [f64; 3]) -> impl Strategy<Value = (Vec3, Vec3)> {
    (arb_point_in_room(dims), arb_point_in_room(dims))
        .prop_filter("points must not coincide", |(a, b)| a.distance(*b) > 0.2)
}

/// A shoebox room mixing preset boundary materials.
pub fn arb_shoebox_room() -> impl Strategy<Value = ShoeboxRoom> {
    let presets = || {
        prop::sample::select(vec![
            BoundaryMaterial::rigid(),
            BoundaryMaterial::concrete_painted(),
            BoundaryMaterial::curtain_heavy(),
            BoundaryMaterial::glass(),
            BoundaryMaterial::absorber_panel(),
        ])
    };
    (arb_shoebox_dims(), presets(), presets(), presets()).prop_map(|(dims, a, b, c)| {
        // Three materials spread over the six walls keeps variety while
        // every wall stays a valid preset.
        ShoeboxRoom::new(dims, [a.clone(), a, b.clone(), b, c.clone(), c]).expect("valid room")
    })
}

/// Absorption coefficient in [0, 0.95] (strictly lossy, no singularity).
pub fn arb_absorption() -> impl Strategy<Value = f64> {
    arb_f64_range(0.0, 0.95)
}

/// Frequency in the geometric-acoustics valid band (100 Hz – 8 kHz).
pub fn arb_frequency() -> impl Strategy<Value = f64> {
    arb_f64_range(100.0, 8000.0)
}

/// Passive specific impedance: complex with positive real part (dissipative)
/// and magnitude in [0.1, 100] × the given reference `rho_c`.
pub fn arb_passive_impedance(rho_c: f64) -> impl Strategy<Value = (f64, f64)> {
    (arb_f64_range(0.05, 100.0), arb_f64_range(-100.0, 100.0)).prop_map(move |(mag, phase_lim)| {
        // Positive-real check: re = mag·cos(phase) > 0 enforced by construction
        // through a phase in (−π/2, π/2).
        let phase = phase_lim * core::f64::consts::FRAC_PI_2 / 100.0;
        let re = mag * phase.cos();
        let im = mag * phase.sin();
        (re.max(1.0e-3) * rho_c, im * rho_c)
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::strategy::ValueTree;

    #[test]
    fn strategies_generate_valid_data() {
        let mut runner = proptest::test_runner::TestRunner::deterministic();
        // Rooms are valid.
        for _ in 0..50 {
            let dims = arb_shoebox_dims().new_tree(&mut runner).unwrap().current();
            assert!(dims.iter().all(|d| (2.0..=30.0).contains(d)));
            // Points are strictly inside those dims.
            let p = arb_point_in_room(dims)
                .new_tree(&mut runner)
                .unwrap()
                .current();
            assert!(p.x > 0.0 && p.x < dims[0]);
            assert!(p.y > 0.0 && p.y < dims[1]);
            assert!(p.z > 0.0 && p.z < dims[2]);
            let (a, b) = arb_point_pair_in_room(dims)
                .new_tree(&mut runner)
                .unwrap()
                .current();
            assert!(a.distance(b) > 0.2);
        }
        // Passive impedances have positive real part.
        for _ in 0..50 {
            let (re, _) = arb_passive_impedance(413.0)
                .new_tree(&mut runner)
                .unwrap()
                .current();
            assert!(re > 0.0);
        }
        // Rooms build.
        let _room = arb_shoebox_room().new_tree(&mut runner).unwrap().current();
    }
}
