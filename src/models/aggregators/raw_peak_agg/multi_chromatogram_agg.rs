use super::chromatogram_agg::{
    ChromatomobilogramStatsArrays, MappingCollection, ScanTofStatsCalculatorPair,
};
use crate::models::frames::raw_peak::RawPeak;
use crate::traits::aggregator::Aggregator;
use serde::Serialize;

#[derive(Debug, Clone)]
pub struct MultiCMGStats {
    pub scan_tof_mapping: MappingCollection<(usize, u32), ScanTofStatsCalculatorPair>,
    pub id: u64,
}

impl MultiCMGStats {
    pub fn new(id: u64) -> Self {
        Self {
            scan_tof_mapping: MappingCollection::new(),
            id,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct MultiCMGStatsArrays {
    pub transition_stats: MappingCollection<usize, ChromatomobilogramStatsArrays>,
    id: u64,
}

impl Aggregator<(RawPeak, usize)> for MultiCMGStats {
    type Output = MultiCMGStatsArrays;

    fn add(&mut self, peak: &(RawPeak, usize)) {
        let (peak, transition) = peak;
        let u64_intensity = peak.intensity as u64;
        let rt_miliseconds = (peak.retention_time * 1000.0) as u32;

        self.scan_tof_mapping
            .entry((*transition, rt_miliseconds))
            .and_modify(|curr| {
                curr.add(u64_intensity, peak.scan_index, peak.tof_index);
            })
            .or_insert(ScanTofStatsCalculatorPair::new(
                u64_intensity,
                peak.scan_index,
                peak.tof_index,
            ));
    }

    fn finalize(self) -> MultiCMGStatsArrays {
        let mut transition_stats = MappingCollection::new();
        for ((transition, rt_ms), scan_tof_mapping) in self.scan_tof_mapping.into_iter() {
            transition_stats
                .entry(transition)
                .and_modify(|curr: &mut ChromatomobilogramStatsArrays| {
                    curr.retention_time_miliseconds.push(rt_ms);
                    curr.scan_index_means
                        .push(scan_tof_mapping.scan.mean().unwrap());
                    curr.scan_index_sds
                        .push(scan_tof_mapping.scan.standard_deviation().unwrap());
                    curr.tof_index_means
                        .push(scan_tof_mapping.tof.mean().unwrap());
                    curr.tof_index_sds
                        .push(scan_tof_mapping.tof.standard_deviation().unwrap());
                })
                .or_insert_with(|| {
                    let mut out = ChromatomobilogramStatsArrays::new();
                    out.retention_time_miliseconds.push(rt_ms);
                    out.scan_index_means
                        .push(scan_tof_mapping.scan.mean().unwrap());
                    out.scan_index_sds
                        .push(scan_tof_mapping.scan.standard_deviation().unwrap());
                    out.tof_index_means
                        .push(scan_tof_mapping.tof.mean().unwrap());
                    out.tof_index_sds
                        .push(scan_tof_mapping.tof.standard_deviation().unwrap());
                    out
                });
        }

        MultiCMGStatsArrays {
            transition_stats,
            id: self.id,
        }
    }
}
