use std::fmt::Display;
use std::hash::Hash;
use std::sync::Arc;

use timsrust::{AcquisitionType, Frame, MSLevel, QuadrupoleSettings};

use crate::sort_by_indices_multi;
use crate::utils::sorting::argsort_by;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct SingleQuadrupoleSettingIndex {
    pub major_index: usize,
    pub sub_index: usize,
}

impl Display for SingleQuadrupoleSettingIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self)
    }
}

#[derive(Debug, Clone, Copy)]
pub struct SingleQuadrupoleSettingRanges {
    pub scan_start: usize,
    pub scan_end: usize,
    pub isolation_mz: f64,
    pub isolation_width: f64,
    pub isolation_high: f64,
    pub isolation_low: f64,
    pub collision_energy: f64,
}

impl Display for SingleQuadrupoleSettingRanges {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "SingleQuadrupoleSettingRanges {{ scan_start: {}, scan_end: {}, isolation_mz: {}, isolation_width: {}, isolation_high: {}, isolation_low: {}, collision_energy: {} }}",
            self.scan_start, self.scan_end, self.isolation_mz, self.isolation_width, self.isolation_high, self.isolation_low, self.collision_energy
        )
    }
}

#[derive(Debug, Clone, Copy)]
pub struct SingleQuadrupoleSetting {
    pub index: SingleQuadrupoleSettingIndex,
    pub ranges: SingleQuadrupoleSettingRanges,
}

impl Display for SingleQuadrupoleSetting {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "SingleQuadrupoleSetting {{\n index: {},\n ranges: {}\n }}",
            self.index, self.ranges
        )
    }
}

impl Hash for SingleQuadrupoleSetting {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.index.hash(state);
    }
}

impl PartialEq for SingleQuadrupoleSetting {
    fn eq(&self, other: &Self) -> bool {
        let out = self.index == other.index;
        if out && cfg!(debug_assertions) {
            debug_assert!(
                self.ranges.scan_start == other.ranges.scan_start,
                "{:?} != {:?} in scan_start",
                self.ranges.scan_start,
                other.ranges.scan_start
            );
            debug_assert!(self.ranges.scan_end == other.ranges.scan_end);
            debug_assert!(self.ranges.isolation_mz == other.ranges.isolation_mz);
            debug_assert!(self.ranges.isolation_width == other.ranges.isolation_width);
            debug_assert!(self.ranges.isolation_high == other.ranges.isolation_high);
            debug_assert!(self.ranges.isolation_low == other.ranges.isolation_low);
            debug_assert!(self.ranges.collision_energy == other.ranges.collision_energy);
        };
        out
    }
}

impl Eq for SingleQuadrupoleSetting {}

pub fn expand_quad_settings(quad_settings: &QuadrupoleSettings) -> Vec<SingleQuadrupoleSetting> {
    let mut out = Vec::with_capacity(quad_settings.scan_ends.len());
    for i in 0..quad_settings.scan_ends.len() {
        let isolation_width = quad_settings.isolation_width[i];
        let half_width = isolation_width / 2.0;
        let isolation_center = quad_settings.isolation_mz[i];
        let collision_energy = quad_settings.collision_energy[i];
        let scan_start = quad_settings.scan_starts[i];
        let scan_end = quad_settings.scan_ends[i];

        let index = SingleQuadrupoleSettingIndex {
            major_index: quad_settings.index,
            sub_index: i,
        };

        let ranges = SingleQuadrupoleSettingRanges {
            scan_start,
            scan_end,
            isolation_mz: quad_settings.isolation_mz[i],
            isolation_width,
            isolation_high: isolation_center + half_width,
            isolation_low: isolation_center - half_width,
            collision_energy,
        };

        out.push(SingleQuadrupoleSetting {
            index: index,
            ranges: ranges,
        });
    }
    out
}
