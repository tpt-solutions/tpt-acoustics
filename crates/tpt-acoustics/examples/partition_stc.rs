//! Mass-law transmission loss and the ASTM E413 STC of a single leaf.
//!
//! A 5 kg/m² leaf in air. The diffuse-field mass law gains about 6 dB per
//! doubling of frequency in its mass-controlled range, which this example
//! checks between 1 and 2 kHz. The STC is then rated from the 1/3-octave
//! spectrum.
//!
//! Run with `cargo run --release --example partition_stc`.

use tpt_acoustics::core::medium::Medium;
use tpt_acoustics::core::AcousticError;
use tpt_acoustics::vibro::transmission::{
    mass_law_diffuse_tl, sound_transmission_class, STC_BAND_HZ,
};

fn main() -> Result<(), AcousticError> {
    let air = Medium::air_20c();
    let mass = 5.0;

    let spectrum: Vec<(f64, f64)> = STC_BAND_HZ
        .iter()
        .map(|&f| (f, mass_law_diffuse_tl(mass, f, &air)))
        .collect();

    println!("{:>8}  {:>8}", "f [Hz]", "TL [dB]");
    for &(f, tl) in &spectrum {
        println!("{f:>8.0}  {tl:>8.1}");
    }

    let tl_1k = mass_law_diffuse_tl(mass, 1000.0, &air);
    let tl_2k = mass_law_diffuse_tl(mass, 2000.0, &air);
    let slope = tl_2k - tl_1k;
    println!("slope 1 to 2 kHz    : {slope:.2} dB per octave");
    assert!((4.5..=6.5).contains(&slope), "mass-law slope {slope:.2} dB");

    // A single 5 kg/m² leaf rates in the low 20s: STC tracks the mass-law TL
    // near 500 Hz (about 18 dB here), so a value far above that would signal
    // a bug in the rating search.
    let stc = sound_transmission_class(&spectrum)?;
    println!("STC                 : {stc}");
    assert!(
        (18..=26).contains(&stc),
        "STC {stc} outside the expected range"
    );
    Ok(())
}
