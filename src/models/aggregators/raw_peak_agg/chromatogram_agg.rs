use super::super::streaming_aggregator::RunningStatsCalculator;
use crate::models::frames::raw_peak::RawPeak;
use crate::sort_by_indices_multi;
use crate::traits::aggregator::Aggregator;
use crate::utils::sorting::argsort_by;

use serde::Serialize;
use std::collections::BTreeMap;
use std::collections::HashMap;

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

// type MappingCollection<T1, T2> = BTreeMap<T1, T2>;
pub type MappingCollection<T1, T2> = HashMap<T1, T2>;

#[derive(Debug, Clone)]
pub struct ScanTofStatsCalculatorPair {
    pub scan: RunningStatsCalculator,
    pub tof: RunningStatsCalculator,
}

impl ScanTofStatsCalculatorPair {
    pub fn new(intensity: u64, scan_index: usize, tof_index: u32) -> Self {
        let scan_index = scan_index as f64;
        let tof_index = tof_index as f64;

        Self {
            scan: RunningStatsCalculator::new(intensity, scan_index),
            tof: RunningStatsCalculator::new(intensity, tof_index),
        }
    }
    pub fn add(&mut self, intensity: u64, scan_index: usize, tof_index: u32) {
        self.scan.add(scan_index as f64, intensity);
        self.tof.add(tof_index as f64, intensity);
    }

    pub fn weight(&self) -> u64 {
        self.scan.weight()
    }
}

#[derive(Debug, Clone)]
pub struct ChromatomobilogramStats {
    // TODO OPTIMIZE THIS ... as needed.
    // In theory we can optimize this to make a single aggregator struct
    // that shares the weight (intensity), since all will have the same weight
    // and retention times.
    pub scan_tof_mapping: MappingCollection<u32, ScanTofStatsCalculatorPair>,
    pub id: u64,
}

impl ChromatomobilogramStats {
    pub fn new(id: u64) -> Self {
        Self {
            scan_tof_mapping: MappingCollection::new(),
            id,
        }
    }
}

#[derive(Debug, Clone, Serialize, Default)]
pub struct ChromatomobilogramStatsArrays {
    pub retention_time_miliseconds: Vec<u32>,
    pub tof_index_means: Vec<f64>,
    pub tof_index_sds: Vec<f64>,
    pub scan_index_means: Vec<f64>,
    pub scan_index_sds: Vec<f64>,
    pub intensities: Vec<u64>,
}

impl ChromatomobilogramStatsArrays {
    // TODO use default instead of new everywhere ..
    pub fn new() -> Self {
        Self {
            retention_time_miliseconds: Vec::new(),
            tof_index_means: Vec::new(),
            tof_index_sds: Vec::new(),
            scan_index_means: Vec::new(),
            scan_index_sds: Vec::new(),
            intensities: Vec::new(),
        }
    }

    pub fn fold(&mut self, other: Self) {
        self.retention_time_miliseconds
            .extend(other.retention_time_miliseconds);
        self.tof_index_means.extend(other.tof_index_means);
        self.tof_index_sds.extend(other.tof_index_sds);
        self.scan_index_means.extend(other.scan_index_means);
        self.scan_index_sds.extend(other.scan_index_sds);
        self.intensities.extend(other.intensities);
    }

    pub fn sort_by_rt(&mut self) {
        let mut indices = argsort_by(&self.retention_time_miliseconds, |x| *x);
        sort_by_indices_multi!(
            &mut indices,
            &mut self.retention_time_miliseconds,
            &mut self.tof_index_means,
            &mut self.tof_index_sds,
            &mut self.scan_index_means,
            &mut self.scan_index_sds,
            &mut self.intensities
        );
    }
}

impl Aggregator<RawPeak> for ChromatomobilogramStats {
    type Output = ChromatomobilogramStatsArrays;

    fn add(&mut self, peak: &RawPeak) {
        let u64_intensity = peak.intensity as u64;
        let rt_miliseconds = (peak.retention_time * 1000.0) as u32;

        self.scan_tof_mapping
            .entry(rt_miliseconds)
            .and_modify(|curr| {
                curr.add(u64_intensity, peak.scan_index, peak.tof_index);
            })
            .or_insert(ScanTofStatsCalculatorPair::new(
                u64_intensity,
                peak.scan_index,
                peak.tof_index,
            ));
    }

    fn finalize(self) -> ChromatomobilogramStatsArrays {
        type VecTuple = (Vec<f64>, Vec<f64>);
        type OutVecTubples = ((VecTuple, VecTuple), Vec<u64>);
        let ((scan_data, tof_data), intensities): OutVecTubples = self
            .scan_tof_mapping
            .values()
            .map(|pair| {
                (
                    (
                        (
                            pair.scan.mean().unwrap(),
                            pair.scan.standard_deviation().unwrap(),
                        ),
                        (
                            pair.tof.mean().unwrap(),
                            pair.tof.standard_deviation().unwrap(),
                        ),
                    ),
                    pair.tof.weight(),
                )
            })
            .unzip();

        let (scan_means, scan_sds) = scan_data;
        let (tof_means, tof_sds) = tof_data;

        ChromatomobilogramStatsArrays {
            retention_time_miliseconds: self.scan_tof_mapping.keys().cloned().collect::<Vec<u32>>(),
            tof_index_means: tof_means,
            tof_index_sds: tof_sds,
            scan_index_means: scan_means,
            scan_index_sds: scan_sds,
            intensities,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct ChromatomobilogramVectorArrayTuples {
    pub scan_indices: Vec<(usize, u64)>,
    pub tof_indices: Vec<(u32, u64)>,
    pub retention_times: Vec<(f32, u64)>,
}
