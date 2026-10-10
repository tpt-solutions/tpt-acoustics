//! Octave-band absorption data for common finishes.
//!
//! Each entry gives the random-incidence absorption coefficient in six octave
//! bands, 125 Hz to 4 kHz, with its source. The values are typical
//! laboratory averages from published vendor charts, measured in reverberation
//! rooms (ASTM C423 style). They vary with installation, mounting and product,
//! so treat them as typical values and use the manufacturer's data for a
//! specific product.
//!
//! The database is deliberately short. Entries are added only when a
//! specific, checkable source has been identified. Each value set in this file
//! comes from one named chart, and the chart's name is in `source`.

use crate::material::BoundaryMaterial;

/// The six octave-band centres of the database, Hz.
pub const OCTAVE_BANDS_HZ: [f64; 6] = [125.0, 250.0, 500.0, 1000.0, 2000.0, 4000.0];

/// One finish: name, absorption per band, and where the values come from.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OctaveAbsorption {
    /// Lookup key, lower case with underscores.
    pub name: &'static str,
    /// Absorption coefficient per band, matching [`OCTAVE_BANDS_HZ`].
    pub alpha: [f64; 6],
    /// The published chart the values were read from.
    pub source: &'static str,
}

/// The finishes in the database.
pub const FINISHES: [OctaveAbsorption; 3] = [
    OctaveAbsorption {
        name: "concrete_sealed_or_painted",
        alpha: [0.01, 0.01, 0.02, 0.02, 0.02, 0.02],
        source: "Sound Absorption Coefficient Chart (125 Hz–4 kHz), commercial-acoustics.com",
    },
    OctaveAbsorption {
        name: "carpet_on_concrete",
        alpha: [0.02, 0.06, 0.14, 0.37, 0.60, 0.65],
        source: "Sound Absorption Coefficients Guide, commercial-acoustics.com",
    },
    OctaveAbsorption {
        name: "suspended_acoustic_ceiling_tile",
        alpha: [0.25, 0.45, 0.78, 0.92, 0.90, 0.85],
        source: "Sound Absorption Coefficient Chart (125 Hz–4 kHz), commercial-acoustics.com",
    },
];

/// Looks up a finish by name.
#[must_use]
pub fn lookup(name: &str) -> Option<&'static OctaveAbsorption> {
    FINISHES.iter().find(|f| f.name == name)
}

impl OctaveAbsorption {
    /// The finish as a wall material: the six band values, with zero
    /// scattering.
    #[must_use]
    pub fn to_material(&self) -> BoundaryMaterial {
        let absorption: alloc::vec::Vec<(f64, f64)> = OCTAVE_BANDS_HZ
            .iter()
            .zip(self.alpha)
            .map(|(&f, a)| (f, a))
            .collect();
        let scattering: alloc::vec::Vec<(f64, f64)> =
            OCTAVE_BANDS_HZ.iter().map(|&f| (f, 0.0)).collect();
        // The table is valid by construction (values in [0, 1], bands ascending).
        BoundaryMaterial::new(self.name, absorption, scattering)
            .expect("database entries are valid by construction")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn every_entry_is_a_valid_absorption_profile() {
        for f in &FINISHES {
            assert!(
                f.alpha.iter().all(|a| (0.0..=1.0).contains(a)),
                "{}",
                f.name
            );
            assert!(!f.source.is_empty());
        }
    }

    #[test]
    fn lookup_finds_entries_and_rejects_unknown_names() {
        assert!(lookup("carpet_on_concrete").is_some());
        assert!(lookup("unobtainium").is_none());
    }

    #[test]
    fn material_reproduces_the_band_values() {
        let f = lookup("suspended_acoustic_ceiling_tile").unwrap();
        let m = f.to_material();
        for (&band, &a) in OCTAVE_BANDS_HZ.iter().zip(&f.alpha) {
            assert!((m.absorption_at(band) - a).abs() < 1e-12, "band {band}");
        }
    }
}
