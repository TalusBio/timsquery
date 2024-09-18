use crate::models::frames::raw_peak::RawPeak;
use crate::traits::aggregator::Aggregator;

pub struct RawPeakIntensityAggregator {
    pub intensity: u64,
}

impl Aggregator<RawPeak> for RawPeakIntensityAggregator {
    type Output = u64;

    fn add(&mut self, peak: &RawPeak) {
        self.intensity += peak.intensity as u64;
    }

    fn fold(&mut self, other: Self) {
        self.intensity += other.intensity;
    }

    fn finalize(&self) -> u64 {
        self.intensity
    }
}
