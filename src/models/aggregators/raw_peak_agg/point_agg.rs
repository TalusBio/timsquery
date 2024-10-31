use crate::models::frames::raw_peak::RawPeak;
use crate::traits::aggregator::{Aggregator, NoContext, ProvidesContext};
use serde::Serialize;

#[derive(Debug, Clone, Copy)]
pub struct RawPeakIntensityAggregator<T> {
    pub id: u64,
    pub intensity: u64,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: Send + Sync> RawPeakIntensityAggregator<T> {
    pub fn new(id: u64) -> Self {
        Self {
            id,
            intensity: 0,
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: Send + Sync + Clone> Aggregator for RawPeakIntensityAggregator<T> {
    type Item = RawPeak;
    type Context = NoContext;
    type Output = u64;

    fn add(&mut self, peak: impl Into<RawPeak>) {
        let peak = peak.into();
        self.intensity += peak.intensity as u64;
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

impl Aggregator for RawPeakVectorAggregator {
    type Item = RawPeak;
    type Context = NoContext;
    type Output = RawPeakVectorArrays;

    fn add(&mut self, peak: impl Into<RawPeak>) {
        let peak = peak.into();
        self.peaks.scans.push(peak.scan_index);
        self.peaks.tofs.push(peak.tof_index);
        self.peaks.intensities.push(peak.intensity);
        self.peaks.retention_times.push(peak.retention_time);
    }

    fn finalize(self) -> RawPeakVectorArrays {
        self.peaks
    }
}
