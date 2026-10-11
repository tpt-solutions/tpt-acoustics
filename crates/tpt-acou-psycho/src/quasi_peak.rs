//! Quasi-peak detector model fitted to the tone-burst limits of ITU-R BS.468-4.
//!
//! **What this is.** A two-channel peak-hold model. Each channel charges with
//! its own time constant while the input is above its output, and discharges
//! with its own time constant otherwise. The two channels are summed with
//! weights that add to 1, so a steady input gives a steady reading of 1. The
//! parameters were found by a numerical search, so that the readings meet the
//! limits of BS.468-4 Tables 2 and 3 (single and repeated 5 kHz tone-bursts).
//!
//! **What this is not.** BS.468-4 gives no circuit and no time constants. It
//! defines the detector by tone-burst limits (§2.1 and §2.2). The internal
//! structure here is ours. Only the burst limits are checked, so this is not
//! a BS.468-4 instrument and makes no claim about the other requirements of
//! §2 (overload, reversibility, overswing, calibration, input impedance).
//!
//! The limits come from the freely downloadable ITU-R BS.468-4 recommendation.
//! The tests in this module use those limits and drive the sampled detector,
//! so they check the implementation, not only the fitted formulas.

use alloc::vec::Vec;

use tpt_acou_core::error::AcousticError;

/// Parameters of the two-channel peak-hold model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct QuasiPeakParams {
    /// Weight of channel 1; channel 2 has weight `1 − weight`.
    pub weight: f64,
    /// Charge time constant of channel 1 (s).
    pub charge_s1: f64,
    /// Charge time constant of channel 2 (s).
    pub charge_s2: f64,
    /// Discharge time constant of channel 1 (s).
    pub discharge_s1: f64,
    /// Discharge time constant of channel 2 (s).
    pub discharge_s2: f64,
}

impl QuasiPeakParams {
    /// Parameters found by the search described in the module docs.
    #[must_use]
    pub const fn bs468_tone_burst_fit() -> Self {
        Self {
            weight: 0.389,
            charge_s1: 1.755e-3,
            charge_s2: 38.26e-3,
            discharge_s1: 190.1e-3,
            discharge_s2: 880.8e-3,
        }
    }
}

/// A sampled quasi-peak detector.
#[derive(Debug, Clone, PartialEq)]
pub struct QuasiPeakDetector {
    params: QuasiPeakParams,
    sample_rate_hz: f64,
    y1: f64,
    y2: f64,
}

impl QuasiPeakDetector {
    /// A detector at rest, sampling its input at `sample_rate_hz`.
    ///
    /// # Errors
    /// Returns [`AcousticError::InvalidArgument`] if the sample rate is not
    /// positive and finite, a time constant is not positive and finite, or the
    /// weight is outside `[0, 1]`.
    pub fn new(params: QuasiPeakParams, sample_rate_hz: f64) -> Result<Self, AcousticError> {
        if !(sample_rate_hz.is_finite() && sample_rate_hz > 0.0) {
            return Err(AcousticError::InvalidArgument(
                "sample rate must be finite and > 0",
            ));
        }
        let taus = [
            params.charge_s1,
            params.charge_s2,
            params.discharge_s1,
            params.discharge_s2,
        ];
        if taus.iter().any(|t| !(t.is_finite() && *t > 0.0)) {
            return Err(AcousticError::InvalidArgument(
                "time constants must be finite and > 0",
            ));
        }
        if !(params.weight.is_finite() && (0.0..=1.0).contains(&params.weight)) {
            return Err(AcousticError::InvalidArgument("weight must lie in [0, 1]"));
        }
        Ok(Self {
            params,
            sample_rate_hz,
            y1: 0.0,
            y2: 0.0,
        })
    }

    /// Resets both channels to zero.
    pub fn reset(&mut self) {
        self.y1 = 0.0;
        self.y2 = 0.0;
    }

    /// Advances one sample of a non-negative envelope and returns the reading.
    ///
    /// Each channel moves toward `x` with the charge constant while `x` is above
    /// its state, and with the discharge constant otherwise. Each step is the
    /// exact exponential update for a constant input over one sample.
    pub fn step(&mut self, x: f64) -> f64 {
        let dt = 1.0 / self.sample_rate_hz;
        let p = &self.params;
        self.y1 = advance(self.y1, x, dt, p.charge_s1, p.discharge_s1);
        self.y2 = advance(self.y2, x, dt, p.charge_s2, p.discharge_s2);
        p.weight * self.y1 + (1.0 - p.weight) * self.y2
    }

    /// Runs the detector over an envelope and returns the reading at each
    /// sample. The state carries over from previous calls.
    pub fn run(&mut self, envelope: &[f64]) -> Vec<f64> {
        envelope.iter().map(|&x| self.step(x)).collect()
    }
}

/// One exact exponential step toward `x`, with charge or discharge constant.
fn advance(y: f64, x: f64, dt: f64, charge_s: f64, discharge_s: f64) -> f64 {
    let tau = if x > y { charge_s } else { discharge_s };
    let k = 1.0 - libm::exp(-dt / tau);
    y + (x - y) * k
}

#[cfg(test)]
mod tests {
    use super::*;

    const FS: f64 = 100_000.0;

    /// BS.468-4 Table 2: burst duration (ms), lower and upper limits (% of full
    /// scale) with the steady signal at 80 %.
    const TABLE_2: [(f64, f64, f64); 8] = [
        (1.0, 13.5, 21.4),
        (2.0, 22.4, 31.6),
        (5.0, 34.0, 46.0),
        (10.0, 41.0, 55.0),
        (20.0, 44.0, 60.0),
        (50.0, 50.0, 68.0),
        (100.0, 58.0, 78.0),
        (200.0, 68.0, 92.0),
    ];

    /// BS.468-4 Table 3: bursts per second, lower and upper limits (% of the
    /// steady reading) for 5 ms bursts.
    const TABLE_3: [(f64, f64, f64); 3] =
        [(2.0, 43.0, 53.0), (10.0, 72.0, 82.0), (100.0, 94.0, 100.0)];

    fn detector() -> QuasiPeakDetector {
        QuasiPeakDetector::new(QuasiPeakParams::bs468_tone_burst_fit(), FS).unwrap()
    }

    #[test]
    fn steady_input_reads_one() {
        let mut d = detector();
        let out = d.run(&alloc::vec![1.0; 2 * FS as usize]);
        assert!((out.last().unwrap() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn single_bursts_meet_table_2() {
        for (ms, lo, hi) in TABLE_2 {
            let mut d = detector();
            let burst = (ms * 1e-3 * FS).round() as usize;
            let mut env = alloc::vec![1.0; burst];
            env.resize(burst + (0.3 * FS) as usize, 0.0);
            let reading = d.run(&env).into_iter().fold(0.0_f64, f64::max);
            // The steady signal reads 80 % of full scale, so scale by 80.
            let pct = 80.0 * reading;
            assert!(
                lo <= pct && pct <= hi,
                "{ms} ms: {pct:.1} % outside [{lo}, {hi}]"
            );
        }
    }

    #[test]
    fn repeated_bursts_meet_table_3() {
        for (rate, lo, hi) in TABLE_3 {
            let mut d = detector();
            let total = (6.0 * FS) as usize;
            let period = (FS / rate).round() as usize;
            let burst = (5e-3 * FS).round() as usize;
            let mut env = alloc::vec![0.0; total];
            let mut start = 0;
            while start < total {
                for v in env.iter_mut().skip(start).take(burst) {
                    *v = 1.0;
                }
                start += period;
            }
            let out = d.run(&env);
            // Steady state: the last second of the run, relative to a steady input (1.0).
            let steady = out[total - FS as usize..]
                .iter()
                .copied()
                .fold(0.0_f64, f64::max);
            let pct = 100.0 * steady;
            assert!(
                lo <= pct && pct <= hi,
                "{rate}/s: {pct:.1} % outside [{lo}, {hi}]"
            );
        }
    }

    #[test]
    fn rejects_invalid_configuration() {
        let p = QuasiPeakParams::bs468_tone_burst_fit();
        assert!(QuasiPeakDetector::new(p, 0.0).is_err());
        assert!(QuasiPeakDetector::new(QuasiPeakParams { weight: 1.5, ..p }, FS).is_err());
        assert!(QuasiPeakDetector::new(
            QuasiPeakParams {
                charge_s1: 0.0,
                ..p
            },
            FS
        )
        .is_err());
    }

    #[test]
    fn output_stays_within_the_input_range() {
        let mut d = detector();
        let env: Vec<f64> = (0..10_000)
            .map(|i| if (i / 500) % 2 == 0 { 0.7 } else { 0.0 })
            .collect();
        for y in d.run(&env) {
            assert!((-1e-12..=0.7 + 1e-12).contains(&y), "{y}");
        }
    }
}
