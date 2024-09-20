use crate::models::frames::raw_peak::RawPeak;
use crate::traits::aggregator::Aggregator;
use serde::Serialize;

#[derive(Debug, Clone, Copy)]
pub struct RawPeakIntensityAggregator {
    pub id: u64,
    pub intensity: u64,
}

impl RawPeakIntensityAggregator {
    pub fn new(id: u64) -> Self {
        Self { id, intensity: 0 }
    }
}

impl Aggregator<RawPeak> for RawPeakIntensityAggregator {
    type Output = u64;

    fn add(&mut self, peak: &RawPeak) {
        self.intensity += peak.intensity as u64;
    }

    fn fold(&mut self, other: Self) {
        self.intensity += other.intensity;
    }

    fn finalize(self) -> u64 {
        self.intensity
    }
}

#[derive(Debug, Clone)]
pub struct RawPeakVectorAggregator {
    pub id: u64,
    pub peaks: RawPeakVectorArrays,
}

impl RawPeakVectorAggregator {
    pub fn new(id: u64) -> Self {
        Self {
            id,
            peaks: RawPeakVectorArrays::new(),
        }
    }
}

impl RawPeakVectorArrays {
    pub fn new() -> Self {
        Self {
            scans: Vec::new(),
            tofs: Vec::new(),
            intensities: Vec::new(),
            retention_times: Vec::new(),
        }
    }
}

#[derive(Debug, Clone, Serialize, Default)]
pub struct RawPeakVectorArrays {
    pub scans: Vec<usize>,
    pub tofs: Vec<u32>,
    pub intensities: Vec<u32>,
    pub retention_times: Vec<f32>,
}

impl Aggregator<RawPeak> for RawPeakVectorAggregator {
    type Output = RawPeakVectorArrays;

    fn add(&mut self, peak: &RawPeak) {
        self.peaks.scans.push(peak.scan_index);
        self.peaks.tofs.push(peak.tof_index);
        self.peaks.intensities.push(peak.intensity);
        self.peaks.retention_times.push(peak.retention_time);
    }

    fn fold(&mut self, other: Self) {
        self.peaks.scans.extend(other.peaks.scans);
        self.peaks.tofs.extend(other.peaks.tofs);
        self.peaks.intensities.extend(other.peaks.intensities);
        self.peaks
            .retention_times
            .extend(other.peaks.retention_times);
    }

    fn finalize(self) -> RawPeakVectorArrays {
        self.peaks
    }
}
