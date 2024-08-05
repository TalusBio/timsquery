use core::f32;

use serde::{Deserialize, Serialize};

use crate::models::elution_group::ElutionGroup;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MzToleramce {
    Absolute((f64, f64)),
    Ppm((f64, f64)),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum RtTolerance {
    Absolute((f32, f32)),
    Pct((f32, f32)),
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum MobilityTolerance {
    Absolute((f32, f32)),
    Pct((f32, f32)),
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
}

impl Default for DefaultTolerance {
    fn default() -> Self {
        DefaultTolerance {
            ms: MzToleramce::Ppm((20.0, 20.0)),
            rt: RtTolerance::Absolute((5.0, 5.0)),
            mobility: MobilityTolerance::Pct((3.0, 3.0)),
            quad: QuadTolerance::Absolute((0.1, 0.1, 1)),
        }
    }
}

pub trait Tolerance {
    fn mz_range(&self, mz: f64) -> (f64, f64);
    fn rt_range(&self, rt: f32) -> (f32, f32);
    fn mobility_range(&self, mobility: f32) -> (f32, f32);
    fn quad_range(&self, precursor_mz: f64, precursor_charge: u8) -> (f32, f32);
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

    fn rt_range(&self, rt: f32) -> (f32, f32) {
        match self.rt {
            RtTolerance::Absolute((low, high)) => (rt - low, rt + high),
            RtTolerance::Pct((low, high)) => {
                let low = rt * low / 100.0;
                let high = rt * high / 100.0;
                (rt - low, rt + high)
            }
        }
    }

    fn mobility_range(&self, mobility: f32) -> (f32, f32) {
        match self.mobility {
            MobilityTolerance::Absolute((low, high)) => (mobility - low, mobility + high),
            MobilityTolerance::Pct((low, high)) => {
                let low = mobility * low / 100.0;
                let high = mobility * high / 100.0;
                (mobility - low, mobility + high)
            }
        }
    }

    fn quad_range(&self, precursor_mz: f64, precursor_charge: u8) -> (f32, f32) {
        match self.quad {
            QuadTolerance::Absolute((low, high, num_isotopes)) => {
                let max_offset = (1.0 / precursor_charge as f32) * num_isotopes as f32;
                let f32_mz = precursor_mz as f32;
                (f32_mz - low - max_offset, f32_mz + high + max_offset)
            }
        }
    }
}

pub trait ToleranceAdapter<QF> {
    fn query_from_elution_group(&self, tol: &dyn Tolerance, elution_group: &ElutionGroup) -> QF;
}
