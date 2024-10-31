use core::f32;
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MzToleramce {
    Absolute((f64, f64)),
    Ppm((f64, f64)),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum RtTolerance {
    Absolute((f32, f32)),
    Pct((f32, f32)),
    None,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MobilityTolerance {
    Absolute((f32, f32)),
    Pct((f32, f32)),
    None,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum QuadTolerance {
    Absolute((f32, f32, u8)),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DefaultTolerance {
    pub ms: MzToleramce,
    pub rt: RtTolerance,
    pub mobility: MobilityTolerance,
    pub quad: QuadTolerance,
    pub num_ms1_isotopes: usize,
    pub num_ms2_isotopes: usize,
}

impl Default for DefaultTolerance {
    fn default() -> Self {
        DefaultTolerance {
            ms: MzToleramce::Ppm((20.0, 20.0)),
            rt: RtTolerance::Absolute((5.0, 5.0)),
            mobility: MobilityTolerance::Pct((3.0, 3.0)),
            quad: QuadTolerance::Absolute((0.1, 0.1, 1)),
            num_ms1_isotopes: 3,
            num_ms2_isotopes: 1,
        }
    }
}

fn mass_to_isotope_mzs(mass: f64, charge: u8, num_isotopes: usize) -> Vec<f64> {
    let mut out = Vec::with_capacity(num_isotopes);
    const PROTON_MASS: f64 = 1.007276;
    for i in 1..(num_isotopes + 1) {
        let mass = mass + i as f64 * PROTON_MASS;
        let mz = mass / charge as f64;
        out.push(mz);
    }
    out
}

fn mz_to_isotope_mzs(mz: f64, charge: u8, num_isotopes: usize) -> Vec<f64> {
    let mut out = Vec::with_capacity(num_isotopes);
    const PROTON_MASS: f64 = 1.007276;
    let proton_mass_frac = PROTON_MASS / charge as f64;
    for i in 0..num_isotopes {
        let mz = mz + (i as f64 * proton_mass_frac);
        out.push(mz);
    }
    out
}

pub trait Tolerance {
    fn mz_range(&self, mz: f64) -> (f64, f64);
    fn rt_range(&self, rt: f32) -> Option<(f32, f32)>;
    fn mobility_range(&self, mobility: f32) -> Option<(f32, f32)>;
    fn quad_range(&self, precursor_mz: f64, precursor_charge: u8) -> (f32, f32);
    fn num_ms1_isotopes(&self) -> usize;
    fn isotope_mzs_mass(&self, monoisotopic_mass: f64, charge: u8) -> Vec<f64> {
        mass_to_isotope_mzs(monoisotopic_mass, charge, self.num_ms1_isotopes())
    }
    fn isotope_mzs_mz(&self, mz: f64, charge: u8) -> Vec<f64> {
        mz_to_isotope_mzs(mz, charge, self.num_ms1_isotopes())
    }
    fn num_ms2_isotopes(&self) -> usize;
}

impl Tolerance for DefaultTolerance {
    fn mz_range(&self, mz: f64) -> (f64, f64) {
        match self.ms {
            MzToleramce::Absolute((low, high)) => (mz - low, mz + high),
            MzToleramce::Ppm((low, high)) => {
                let low = mz * low / 1e6;
                let high = mz * high / 1e6;
                (mz - low, mz + high)
            }
        }
    }

    fn rt_range(&self, rt: f32) -> Option<(f32, f32)> {
        match self.rt {
            RtTolerance::Absolute((low, high)) => Some((rt - low, rt + high)),
            RtTolerance::Pct((low, high)) => {
                let low = rt * low / 100.0;
                let high = rt * high / 100.0;
                Some((rt - low, rt + high))
            }
            RtTolerance::None => None,
        }
    }

    fn mobility_range(&self, mobility: f32) -> Option<(f32, f32)> {
        match self.mobility {
            MobilityTolerance::Absolute((low, high)) => Some((mobility - low, mobility + high)),
            MobilityTolerance::Pct((low, high)) => {
                let low = mobility * (low / 100.0);
                let high = mobility * (high / 100.0);
                Some((mobility - low, mobility + high))
            }
            MobilityTolerance::None => None,
        }
    }

    fn quad_range(&self, precursor_mz: f64, precursor_charge: u8) -> (f32, f32) {
        // Should this be a recoverable error?
        if precursor_charge == 0 {
            panic!("Precursor charge is 0, inputs: self: {:?}, precursor_mz: {:?}, precursor_charge: {:?}", self, precursor_mz, precursor_charge);
        };
        match self.quad {
            QuadTolerance::Absolute((low, high, num_isotopes)) => {
                let max_offset = (1.0 / precursor_charge as f32) * num_isotopes as f32;
                let f32_mz = precursor_mz as f32;
                let mz_low = f32_mz - low - max_offset;
                let mz_high = f32_mz + high + max_offset;
                assert!(mz_low <= mz_high);
                assert!(
                    mz_low > 0.0,
                    "Precursor mz is 0 or less, inputs: self: {:?}, precursor_mz: {:?}, precursor_charge: {:?}",
                    self,
                    precursor_mz,
                    precursor_charge
                );
                (mz_low, mz_high)
            }
        }
    }

    fn num_ms1_isotopes(&self) -> usize {
        self.num_ms1_isotopes
    }

    fn num_ms2_isotopes(&self) -> usize {
        self.num_ms2_isotopes
    }
}

/// A trait that can be implemented by types that can convert
/// elution groups into queries.
///
/// The elution group here is generic but most regularly it will be
/// an `ElutionGroup` struct.
pub trait ToleranceAdapter<QF, T> {
    fn query_from_elution_group(&self, tol: &dyn Tolerance, elution_group: &T) -> QF;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_isotope_mzs_neutral() {
        let test_vals = vec![
            (100.0, 1, vec![101.0, 102.0, 103.0]),
            (100.0, 2, vec![100.5, 101.0, 101.5]),
            (100.0, 3, vec![100.33333, 100.666666, 101.0]),
        ];

        for (monoisotopic_mass, charge, expected) in test_vals {
            let out = mass_to_isotope_mzs(monoisotopic_mass, charge, 3);
            assert_eq!(out.len(), expected.len());
            let abs_diff: Vec<f64> = out
                .iter()
                .zip(expected.iter())
                .map(|(a, b)| (a - b).abs())
                .collect();
            for ad in abs_diff.iter() {
                // Very tight tolerances here ...
                assert!(*ad < 0.01);
            }
        }
    }

    #[test]
    fn test_isotope_mzs_mz() {
        let test_vals = vec![
            (100.0, 1, vec![10.0, 101.0, 102.0]),
            (100.0, 2, vec![100.0, 100.5, 101.5]),
            (100.0, 3, vec![100.0, 100.3333, 101.66666]),
        ];

        for (monoisotopic_mass, charge, expected) in test_vals {
            let out = mz_to_isotope_mzs(monoisotopic_mass, charge, 3);
            assert_eq!(out.len(), expected.len());
            let abs_diff: Vec<f64> = out
                .iter()
                .zip(expected.iter())
                .map(|(a, b)| (a - b).abs())
                .collect();

            for ad in abs_diff.iter() {
                // Very tight tolerances here ...
                assert!(*ad < 0.01);
            }
        }
    }
}
