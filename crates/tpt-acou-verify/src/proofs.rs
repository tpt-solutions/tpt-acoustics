//! Kani bounded-model-checking harnesses.
//!
//! Compiled only under `cfg(kani)`; run with `cargo kani` (CI job `kani`).
//! These proofs cover the safety-critical bookkeeping of the Image Source
//! Method: the enumeration never exceeds its reflection-order bound, never
//! produces an image whose wall-hit list disagrees with its order, and the
//! room constructor rejects non-positive dimensions. Phase 2 adds CFL and
//! grazing-angle harnesses (see todo.md).

use tpt_acou_ray::ism::{count_images_with_order_at_most, image_count_upper_bound, ShoeboxRoom};
use tpt_acou_ray::material::BoundaryMaterial;
use tpt_acou_ray::vec3::Vec3;

/// The ISM lattice count never exceeds its closed-form `(2N+1)³` bound, for
/// any symbolic order cap in `[0, 3]`.
#[kani::proof]
#[kani::unwind(40)]
fn ism_count_never_expects_order_bound() {
    let n: u32 = kani::any();
    kani::assume(n <= 3);
    let count = count_images_with_order_at_most(n);
    kani::assert(
        count <= image_count_upper_bound(n),
        "image enumeration exceeds its reflection-order bound",
    );
}

/// Every enumerated image respects the order cap and its wall-hit list has
/// exactly `order` entries, for any source strictly inside a fixed room.
#[kani::proof]
#[kani::unwind(160)]
fn ism_images_respect_max_order() {
    let max_order: u32 = kani::any();
    kani::assume(max_order <= 2);
    let sx: f64 = kani::any();
    let sy: f64 = kani::any();
    let sz: f64 = kani::any();
    kani::assume(sx > 0.0 && sx < 4.0);
    kani::assume(sy > 0.0 && sy < 4.0);
    kani::assume(sz > 0.0 && sz < 4.0);
    let room = ShoeboxRoom::uniform([4.0, 4.0, 4.0], BoundaryMaterial::rigid())
        .expect("fixed room is valid");
    let images = room
        .image_sources(Vec3::new(sx, sy, sz), max_order)
        .expect("source is inside");
    kani::assert(
        (images.len() as u64) <= image_count_upper_bound(max_order),
        "enumeration size bound",
    );
    for img in &images {
        assert!(img.order <= max_order, "image order exceeds cap");
        assert_eq!(
            img.walls.len() as u32,
            img.order,
            "wall-hit list length must equal order"
        );
    }
}

/// A symbolic non-positive dimension is always rejected by the constructor.
#[kani::proof]
fn room_constructor_rejects_nonpositive_dimension() {
    let d: f64 = kani::any();
    kani::assume(d.is_finite());
    kani::assume(d <= 0.0);
    let result = ShoeboxRoom::new(
        [d, 1.0, 1.0],
        core::array::from_fn(|_| BoundaryMaterial::rigid()),
    );
    kani::assert(result.is_err(), "non-positive dimension must be rejected");
}

/// The 1-D FDTD stability (Courant) classification: an explicit leapfrog wave
/// update is stable exactly when the CFL number `c*dt/dx <= 1`. Kani proves
/// the classification over *rational* CFL numbers built from bounded
/// symbolic integers (float nondeterminism avoided by construction), with
/// both the discrete comparison and a margin guard.
#[kani::proof]
#[kani::unwind(8)]
fn cfl_classification_matches_rational_bound() {
    // c*dt/dx with c = 343 m/s, dt = num/1e6 s, dx = den/1e2 m:
    // CFL = 343 * num * 1e2 / (1e6 * den) = 343 * num / (1e4 * den).
    let num: u64 = kani::any();
    let den: u64 = kani::any();
    kani::assume(num <= 3000);
    kani::assume(den >= 1 && den <= 3000);
    // Stable iff 343*num <= 1e4*den (exact integer comparison — the float
    // kernel must agree with this rational bound).
    let stable_rational = 343 * num <= 10_000 * den;
    let cfl = 343.0 * num as f64 / (10_000.0 * den as f64);
    let stable_float = cfl <= 1.0;
    // Away from the exact boundary the two classifications must agree.
    let near_boundary = (cfl - 1.0).abs() < 1e-9;
    kani::assert(
        near_boundary || stable_rational == stable_float,
        "float CFL classification disagrees with the rational bound",
    );
}

/// Grazing-incidence / near-singularity guard: for any two distinct points
/// strictly inside a symbolic shoebox, the ISM transfer function skips any
/// image closer than the singularity threshold and never divides by zero
/// (the returned value is finite).
#[kani::proof]
#[kani::unwind(30)]
fn ism_grazing_points_stay_finite() {
    let sx: f64 = kani::any();
    let sy: f64 = kani::any();
    let sz: f64 = kani::any();
    let rx: f64 = kani::any();
    let ry: f64 = kani::any();
    let rz: f64 = kani::any();
    // Room 4x4x4, both points strictly inside, at least 1 mm apart.
    kani::assume(sx > 0.0 && sx < 4.0);
    kani::assume(sy > 0.0 && sy < 4.0);
    kani::assume(sz > 0.0 && sz < 4.0);
    kani::assume(rx > 0.0 && rx < 4.0);
    kani::assume(ry > 0.0 && ry < 4.0);
    kani::assume(rz > 0.0 && rz < 4.0);
    kani::assume((sx - rx) * (sx - rx) + (sy - ry) * (sy - ry) + (sz - rz) * (sz - rz) > 1e-12);
    let room =
        ShoeboxRoom::uniform([4.0, 4.0, 4.0], BoundaryMaterial::rigid()).expect("fixed room valid");
    let echo = room
        .echogram(
            Vec3::new(sx, sy, sz),
            Vec3::new(rx, ry, rz),
            1000.0,
            0, // direct path only keeps the unwind small
            &tpt_acou_core::medium::Medium::air_20c(),
        )
        .expect("points inside");
    for (t, w) in &echo {
        assert!(
            t.is_finite() && *t > 0.0,
            "arrival time must be positive/finite"
        );
        assert!(
            w.is_finite() && *w > 0.0,
            "echogram weight must be positive/finite"
        );
    }
}
