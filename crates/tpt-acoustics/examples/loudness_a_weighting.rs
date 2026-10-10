//! A-weighting against the IEC 61672 reference points, and the phon/sone map.
//!
//! The nominal A-weighting corrections are −19.1 dB at 100 Hz, 0 dB at 1 kHz
//! and −2.5 dB at 10 kHz. A 40-phon level is 1 sone by definition.
//!
//! Run with `cargo run --example loudness_a_weighting`.

use tpt_acoustics::psycho::loudness::phon_to_sone;
use tpt_acoustics::psycho::weighting::a_weighting_db;

fn main() {
    let reference = [(100.0, -19.1), (1000.0, 0.0), (10_000.0, -2.5)];
    for (f, nominal) in reference {
        let got = a_weighting_db(f);
        println!("A({f:>6.0} Hz) = {got:>6.2} dB  (IEC nominal {nominal:.1} dB)");
        assert!((got - nominal).abs() < 0.1, "A({f}) = {got}");
    }

    let sone = phon_to_sone(40.0);
    println!("40 phon = {sone} sone");
    assert!((sone - 1.0).abs() < 1e-12);
}
