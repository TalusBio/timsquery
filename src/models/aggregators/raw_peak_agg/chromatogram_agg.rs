use serde::Serialize;
use std::collections::BTreeMap;

use super::super::streaming_aggregator::{RunningStatsCalculator, StreamingAggregatorError};
use crate::models::frames::raw_peak::RawPeak;
use crate::traits::aggregator::Aggregator;

#[derive(Debug, Clone, Serialize)]
pub struct ExtractedIonChromatomobilogram {
    pub rt_tree: BTreeMap<u32, u64>,
    pub scan_tree: BTreeMap<usize, u64>,
    pub tof_tree: BTreeMap<u32, u64>,
    pub id: u64,
}

impl ExtractedIonChromatomobilogram {
    pub fn new(id: u64) -> Self {
        Self {
            rt_tree: BTreeMap::new(),
            scan_tree: BTreeMap::new(),
            tof_tree: BTreeMap::new(),
            id,
        }
    }
}

impl Aggregator<RawPeak> for ExtractedIonChromatomobilogram {
    type Output = ChromatomobilogramVectorArrayTuples;

    fn add(&mut self, peak: &RawPeak) {
        let u64_intensity = peak.intensity as u64;

        // In theory I could use a power of 2 to have a better preservation of
        // the precision.
        // TODO make this macro ... right now it feels very repetitive.
        let rt_miliseconds = (peak.retention_time * 1000.0) as u32;

        self.rt_tree
            .entry(rt_miliseconds)
            .and_modify(|curr| *curr += u64_intensity)
            .or_insert(u64_intensity);

        self.scan_tree
            .entry(peak.scan_index)
            .and_modify(|curr| *curr += u64_intensity)
            .or_insert(u64_intensity);

        self.tof_tree
            .entry(peak.tof_index)
            .and_modify(|curr| *curr += u64_intensity)
            .or_insert(u64_intensity);
    }

    fn fold(&mut self, other: Self) {
        panic!("Not Implemented;")
    }

    fn finalize(self) -> ChromatomobilogramVectorArrayTuples {
        ChromatomobilogramVectorArrayTuples {
            scan_indices: self.scan_tree.into_iter().collect(),
            tof_indices: self.tof_tree.into_iter().collect(),
            retention_times: self
                .rt_tree
                .into_iter()
                .map(|(k, v)| ((k as f32) / 1000.0, v))
                .collect(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct ChromatomobilogramStats {
    // TODO OPTIMIZE THIS ... as needed.
    // In theory we can optimize this to make a single aggregator struct
    // that shares the weight (intensity), since all will have the same weight
    // and retention times.
    pub scan_tree: BTreeMap<u32, RunningStatsCalculator>,
    pub tof_tree: BTreeMap<u32, RunningStatsCalculator>,
    pub id: u64,
}

impl ChromatomobilogramStats {
    pub fn new(id: u64) -> Self {
        Self {
            scan_tree: BTreeMap::new(),
            tof_tree: BTreeMap::new(),
            id,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct ChromatomobilogramStatsArrays {
    pub retention_time_miliseconds: Vec<u32>,
    pub tof_index_means: Vec<f64>,
    pub tof_index_sds: Vec<f64>,
    pub scan_index_means: Vec<f64>,
    pub scan_index_sds: Vec<f64>,
}

impl Aggregator<RawPeak> for ChromatomobilogramStats {
    type Output = ChromatomobilogramStatsArrays;

    fn add(&mut self, peak: &RawPeak) {
        let u64_intensity = peak.intensity as u64;
        let rt_miliseconds = (peak.retention_time * 1000.0) as u32;

        // TODO make this macro as well ...
        self.scan_tree
            .entry(rt_miliseconds)
            .and_modify(|curr| curr.add(peak.scan_index as f64, u64_intensity))
            .or_insert(RunningStatsCalculator::new(
                u64_intensity,
                peak.scan_index as f64,
            ));
        self.tof_tree
            .entry(rt_miliseconds)
            .and_modify(|curr| curr.add(peak.tof_index as f64, u64_intensity))
            .or_insert(RunningStatsCalculator::new(
                u64_intensity,
                peak.tof_index as f64,
            ));
    }

    fn fold(&mut self, _other: Self) {
        panic!("Not Implemented;")
    }

    fn finalize(self) -> ChromatomobilogramStatsArrays {
        ChromatomobilogramStatsArrays {
            retention_time_miliseconds: self.scan_tree.keys().map(|x| *x).collect::<Vec<u32>>(),
            tof_index_means: self
                .tof_tree
                .values()
                .map(|x| x.mean().unwrap())
                .collect::<Vec<f64>>(),
            tof_index_sds: self
                .tof_tree
                .values()
                .map(|x| x.standard_deviation().unwrap())
                .collect::<Vec<f64>>(),
            scan_index_means: self
                .scan_tree
                .values()
                .map(|x| x.mean().unwrap())
                .collect::<Vec<f64>>(),
            scan_index_sds: self
                .scan_tree
                .values()
                .map(|x| x.standard_deviation().unwrap())
                .collect::<Vec<f64>>(),
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct ChromatomobilogramVectorArrayTuples {
    pub scan_indices: Vec<(usize, u64)>,
    pub tof_indices: Vec<(u32, u64)>,
    pub retention_times: Vec<(f32, u64)>,
}
