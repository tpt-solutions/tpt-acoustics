#![allow(missing_docs)] // criterion_group! generates an undocumented function
//! Throughput of the ray tracer, the image-source echogram and RIR synthesis.
//!
//! Run with `cargo bench -p tpt-acou-ray`.

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use tpt_acou_core::medium::Medium;
use tpt_acou_ray::ism::{ShoeboxRoom, SPEED_OF_SOUND_REFERENCE};
use tpt_acou_ray::material::BoundaryMaterial;
use tpt_acou_ray::ray::{trace_shoebox, RayTracerConfig};
use tpt_acou_ray::rir::rir_from_echogram;
use tpt_acou_ray::vec3::Vec3;

fn bench_ray_tracer(c: &mut Criterion) {
    let room = ShoeboxRoom::uniform([10.0, 8.0, 6.0], BoundaryMaterial::concrete_painted())
        .expect("valid room");
    let config = RayTracerConfig {
        n_rays: 5_000,
        time_limit_s: 1.0,
        capture_radius: 0.5,
        min_energy: 1e-9,
        seed: 1,
        speed_of_sound: SPEED_OF_SOUND_REFERENCE,
    };
    let src = Vec3::new(2.0, 2.0, 3.0);
    let rcv = Vec3::new(8.0, 6.0, 3.0);
    c.bench_function("ray_shoebox_5k_rays", |b| {
        b.iter(|| {
            trace_shoebox(black_box(&room), src, rcv, 1000.0, black_box(&config))
                .expect("valid run")
        })
    });
}

fn bench_echogram_and_rir(c: &mut Criterion) {
    let room = ShoeboxRoom::uniform([8.0, 5.0, 3.0], BoundaryMaterial::concrete_painted())
        .expect("valid room");
    let air = Medium::air_20c();
    let src = Vec3::new(1.5, 1.2, 1.4);
    let rcv = Vec3::new(6.2, 3.7, 1.6);
    c.bench_function("echogram_order_12", |b| {
        b.iter(|| {
            room.echogram(black_box(src), black_box(rcv), 1000.0, 12, &air)
                .expect("valid echogram")
        })
    });
    let echogram = room.echogram(src, rcv, 1000.0, 12, &air).expect("echogram");
    c.bench_function("rir_from_echogram_48k", |b| {
        b.iter(|| rir_from_echogram(black_box(&echogram), 48_000.0, 1).expect("rir"))
    });
}

criterion_group!(benches, bench_ray_tracer, bench_echogram_and_rir);
criterion_main!(benches);
