//! Phase-1 integration test (spec.txt §6): transmission loss of a simple
//! partition through the full `tpt-fem` + `tpt-acou-vibro` (+ `tpt-acou-wave`)
//! chain.
//!
//! The partition is a simply supported strip (unit width) of 4 mm glass,
//! discretized with Euler–Bernoulli beam elements. `tpt-fem-sparse` assembles
//! K/M, `tpt-fem-modal` extracts modes and harmonic response, the incident
//! (blocked) acoustic pressure drives the structure, and `tpt-acou-wave`'s
//! BEM computes the radiated power of the vibrating strip via
//! `tpt-acou-vibro::radiation`. TL is assembled per 1/3-octave band (BEM for
//! 125–1000 Hz; the diffuse mass law above, where the coarse patch mesh would
//! be beyond the element-per-wavelength limit) and rated with ASTM E413 STC.

use tpt_acou_core::medium::{Medium, Solid};
use tpt_acou_vibro::radiation::radiation_efficiency_baffled;
use tpt_acou_vibro::transmission::{mass_law_diffuse_tl, sound_transmission_class, STC_BAND_HZ};
use tpt_acou_wave::bem::HelmholtzBem;
use tpt_acou_wave::mesh::BoundaryMesh;
use tpt_fem_modal::modal_analysis;
use tpt_fem_sparse::Coo;
use tpt_math_linalg_complex::Complex;

const LENGTH: f64 = 1.5; // strip length (m)
const WIDTH: f64 = 1.0; // strip width (m)
const THICKNESS: f64 = 0.004; // 4 mm
const N_ELEM: usize = 24; // beam elements
const ZETA: f64 = 0.05; // structural damping ratio (robustness choice)

/// Euler–Bernoulli beam element matrices (2 nodes × (w, theta)).
fn beam_element_matrices(ei: f64, rho_a: f64, l: f64) -> ([[f64; 4]; 4], [[f64; 4]; 4]) {
    let ke = [
        [12.0, 6.0 * l, -12.0, 6.0 * l],
        [6.0 * l, 4.0 * l * l, -6.0 * l, 2.0 * l * l],
        [-12.0, -6.0 * l, 12.0, -6.0 * l],
        [6.0 * l, 2.0 * l * l, -6.0 * l, 4.0 * l * l],
    ];
    let me = [
        [156.0, 22.0 * l, 54.0, -13.0 * l],
        [22.0 * l, 4.0 * l * l, 13.0 * l, -3.0 * l * l],
        [54.0, 13.0 * l, 156.0, -22.0 * l],
        [-13.0 * l, -3.0 * l * l, -22.0 * l, 4.0 * l * l],
    ];
    let s = ei / (l * l * l);
    let m = rho_a * l / 420.0;
    (ke.map(|r| r.map(|v| v * s)), me.map(|r| r.map(|v| v * m)))
}

/// Simply supported: w = 0 at both ends (theta free). DOF layout: node i has
/// (w, theta) = (2i, 2i+1). The constrained w DOFs get a large penalty
/// spring to ground so `tpt-fem-modal`'s Lanczos sees a regular system.
fn assemble_strip(ei: f64, rho_a: f64) -> (Coo, Coo) {
    let n_nodes = N_ELEM + 1;
    let mut k = Coo::with_capacity(N_ELEM * 16 + 8);
    let mut m = Coo::with_capacity(N_ELEM * 16 + 8);
    for e in 0..N_ELEM {
        let l = LENGTH / N_ELEM as f64;
        let (ke, me) = beam_element_matrices(ei, rho_a, l);
        for a in 0..4 {
            for b in 0..4 {
                // Local DOF a -> global (node a/2, component a%2).
                let ga = 2 * e + (a / 2) * 2 + (a % 2);
                let gb = 2 * e + (b / 2) * 2 + (b % 2);
                k.push(ga, gb, ke[a][b]);
                m.push(ga, gb, me[a][b]);
            }
        }
    }
    // Stiffness-only penalty (a mass penalty would wreck M's conditioning
    // and the Lanczos positive-definiteness check); constrained modes simply
    // move to very high frequency, far above the retained band.
    let penalty_k = ei / (LENGTH * LENGTH * LENGTH) * 1.0e8;
    let w0 = 0_usize;
    let w_last = 2 * (n_nodes - 1);
    k.push(w0, w0, penalty_k);
    k.push(w_last, w_last, penalty_k);
    (k, m)
}

#[test]
fn partition_transmission_loss_and_stc() {
    let solid = Solid::new(2500.0, 70.0e9, 0.22).unwrap(); // glass
    let ei = solid.plate_bending_stiffness(THICKNESS).unwrap() * WIDTH; // strip bending stiffness
    let rho_a = solid.density * THICKNESS * WIDTH;

    let (k, m) = assemble_strip(ei, rho_a);
    let n_dofs = 2 * (N_ELEM + 1);

    // --- tpt-fem-modal: eigenanalysis -------------------------------------
    let modal = modal_analysis(&k, &m, 0.0, 24, 48, ZETA).unwrap();
    let w1 = modal.omega[0];
    // Analytic simply supported strip: omega1 = (pi/L)^2 sqrt(EI/(rho A)).
    let f1_exact = (core::f64::consts::PI * core::f64::consts::PI / (LENGTH * LENGTH))
        * libm::sqrt(ei / rho_a)
        / core::f64::consts::TAU;
    let rel = (w1 / core::f64::consts::TAU - f1_exact).abs() / f1_exact;
    assert!(
        rel < 0.05,
        "strip f1 {:.3} Hz vs exact {:.3} Hz",
        w1 / core::f64::consts::TAU,
        f1_exact
    );

    let air = Medium::air_20c();
    let p_inc = 1.0_f64; // 1 Pa RMS incident plane wave (normal incidence)

    // --- per-band radiated power and TL -----------------------------------
    let mut spectrum: Vec<(f64, f64)> = Vec::new();
    for &band in STC_BAND_HZ.iter() {
        let f = band;
        if f <= 1000.0 {
            // Blocked (rigid-wall doubled) pressure drives the strip; the
            // force on w-DOFs is the blocked-pressure surface integral.
            // Consistent nodal loads: blocked pressure (2 p_inc) integrated
            // over the strip width with linear shape functions per element.
            let mut force = vec![0.0_f64; n_dofs];
            let dxs = LENGTH / N_ELEM as f64;
            let q = 2.0 * p_inc * WIDTH; // line load per unit length
            for e in 0..N_ELEM {
                force[2 * e] += q * dxs / 2.0;
                force[2 * (e + 1)] += q * dxs / 2.0;
            }
            let resp = modal.frequency_response(&m, &force, &[core::f64::consts::TAU * f]);
            let disp = &resp[0].displacement;
            // Transverse velocities of the w-DOFs.
            let v: Vec<Complex> = (0..=N_ELEM)
                .map(|e| {
                    let (re, im) = disp[2 * e];
                    Complex::new(re, im) * Complex::new(0.0, core::f64::consts::TAU * f)
                })
                .collect();

            // Map strip velocities onto the partition surface mesh.
            let nx = 24_usize;
            let ny = 8_usize;
            let mesh = BoundaryMesh::rectangle_patch(LENGTH, WIDTH, nx, ny).unwrap();
            let bem =
                HelmholtzBem::new(mesh, core::f64::consts::TAU * f / air.speed_of_sound).unwrap();
            let els = bem.mesh().elements();
            let mut v_n = Vec::with_capacity(els.len());
            for el in els {
                let x_frac = (el.centroid.x + LENGTH / 2.0) / LENGTH;
                let idx = (x_frac * N_ELEM as f64).round().clamp(0.0, N_ELEM as f64) as usize;
                v_n.push(v[idx]);
            }
            let eff = radiation_efficiency_baffled(&bem, &v_n, &air).unwrap();
            let area = LENGTH * WIDTH;
            let v2: f64 = els
                .iter()
                .zip(&v_n)
                .map(|(e, v)| v.norm_sqr() * e.area)
                .sum::<f64>()
                / area;
            let rho_c = air.characteristic_impedance();
            let w_rad = eff * rho_c * area * v2;
            // Incident plane-wave power on the partition (RMS pressure).
            let w_inc = p_inc * p_inc / (2.0 * rho_c) * area;
            let tl = 10.0 * libm::log10(w_inc / w_rad.max(1e-30));
            spectrum.push((f, tl));
        } else {
            // Above the BEM mesh's element-per-wavelength limit: exact
            // diffuse mass law for the same surface mass (documented hybrid).
            let tl = mass_law_diffuse_tl(solid.density * THICKNESS, f, &air);
            spectrum.push((f, tl));
        }
    }

    // --- sanity of the spectrum -------------------------------------------
    // At 500 Hz the resonant (FEM) TL should sit within 8 dB of the diffuse
    // mass law for the same surface mass.
    let tl500 = spectrum.iter().find(|(f, _)| *f == 500.0).unwrap().1;
    let ml500 = mass_law_diffuse_tl(solid.density * THICKNESS, 500.0, &air);
    assert!(
        (tl500 - ml500).abs() < 8.0,
        "TL(500) {tl500:.1} dB vs diffuse mass law {ml500:.1} dB"
    );
    // TL grows with frequency overall (mass-law dominated at HF).
    let tl125 = spectrum.iter().find(|(f, _)| *f == 125.0).unwrap().1;
    assert!(
        tl500 > tl125 - 2.0,
        "TL(125) {tl125:.1} should not exceed TL(500) {tl500:.1} much"
    );

    // --- ASTM E413 rating --------------------------------------------------
    let stc = sound_transmission_class(&spectrum).unwrap();
    let stc_masslaw = sound_transmission_class(
        &STC_BAND_HZ
            .iter()
            .map(|&f| (f, mass_law_diffuse_tl(solid.density * THICKNESS, f, &air)))
            .collect::<Vec<_>>(),
    )
    .unwrap();
    assert!(
        (stc - stc_masslaw).abs() <= 8,
        "STC {stc} vs mass-law STC {stc_masslaw}"
    );
}
