use super::super::chromatogram_agg::ScanTofStatsCalculatorPair;
use super::base::{
    ParitionedCMGAggregator,
    PartitionedCMGArrayStats,
    PartitionedCMGScoredStatsArrays,
};
use crate::models::frames::raw_peak::RawPeak;
use crate::models::queries::MsLevelContext;
use crate::traits::aggregator::Aggregator;
use serde::Serialize;
use std::hash::Hash;

use timsrust::converters::{
    Scan2ImConverter,
    Tof2MzConverter,
};

#[derive(Debug, Clone)]
pub struct MultiCMGStatsAgg<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub converters: (Tof2MzConverter, Scan2ImConverter),
    pub ms1_stats: ParitionedCMGAggregator<usize>,
    pub ms2_stats: ParitionedCMGAggregator<FH>,
    pub id: u64,
    pub context: Option<MsLevelContext<usize, FH>>,
    pub buffer: Option<ScanTofStatsCalculatorPair>,
}

#[derive(Debug, Clone)]
pub struct MultiCMGStatsFactory<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub converters: (Tof2MzConverter, Scan2ImConverter),
    pub _phantom: std::marker::PhantomData<FH>,
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync> MultiCMGStatsFactory<FH> {
    pub fn build(&self, id: u64) -> MultiCMGStatsAgg<FH> {
        MultiCMGStatsAgg {
            converters: (self.converters.0, self.converters.1),
            ms1_stats: Default::default(),
            ms2_stats: Default::default(),
            id,
            context: None,
            buffer: None,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct MultiCMGArrayStats<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub ms1_stats: PartitionedCMGArrayStats<usize>,
    pub ms2_stats: PartitionedCMGArrayStats<FH>,
    pub id: u64,
}

#[derive(Debug, Clone, Serialize)]
pub struct NaturalFinalizedMultiCMGStatsArrays<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub ms1_stats: PartitionedCMGScoredStatsArrays<usize>,
    pub ms2_stats: PartitionedCMGScoredStatsArrays<FH>,
    pub id: u64,
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug> Aggregator
    for MultiCMGStatsAgg<FH>
{
    type Context = MsLevelContext<usize, FH>;
    type Item = RawPeak;
    type Output = NaturalFinalizedMultiCMGStatsArrays<FH>;

    fn add(&mut self, peak: impl Into<Self::Item>) {
        let peak = peak.into();
        let u64_intensity = peak.intensity as u64;
        let rt_miliseconds = (peak.retention_time * 1000.0) as u32;

        match self.get_context() {
            MsLevelContext::MS1(i) => {
                self.ms1_stats.uniq_rts.insert(rt_miliseconds);
                self.ms1_stats
                    .scan_tof_mapping
                    .entry((i, rt_miliseconds))
                    .and_modify(|curr| {
                        curr.add(u64_intensity, peak.scan_index, peak.tof_index);
                    })
                    .or_insert(ScanTofStatsCalculatorPair::new(
                        u64_intensity,
                        peak.scan_index,
                        peak.tof_index,
                    ));
            }
            MsLevelContext::MS2(i) => {
                self.ms2_stats.uniq_rts.insert(rt_miliseconds);
                self.ms2_stats
                    .scan_tof_mapping
                    .entry((i, rt_miliseconds))
                    .and_modify(|curr| {
                        curr.add(u64_intensity, peak.scan_index, peak.tof_index);
                    })
                    .or_insert(ScanTofStatsCalculatorPair::new(
                        u64_intensity,
                        peak.scan_index,
                        peak.tof_index,
                    ));
            }
        }
    }

    fn set_context(&mut self, context: Self::Context) {
        match &context {
            MsLevelContext::MS1(i) => {
                self.ms1_stats.uniq_ids.insert(*i);
            }
            MsLevelContext::MS2(i) => {
                self.ms2_stats.uniq_ids.insert(i.clone());
            }
        }
        self.context = Some(context);
    }

    fn supports_context(&self) -> bool {
        true
    }

    fn get_context(&self) -> Self::Context {
        match &self.context {
            Some(context) => context.clone(),
            None => panic!("No context set"),
        }
    }

    fn finalize(self) -> NaturalFinalizedMultiCMGStatsArrays<FH> {
        let mz_converter = &self.converters.0;
        let mobility_converter = &self.converters.1;

        let ms1_stats = PartitionedCMGScoredStatsArrays::new(
            PartitionedCMGArrayStats::from(self.ms1_stats.finalize()),
            mz_converter,
            mobility_converter,
        );
        let ms2_stats = PartitionedCMGScoredStatsArrays::new(
            PartitionedCMGArrayStats::from(self.ms2_stats.finalize()),
            mz_converter,
            mobility_converter,
        );

        NaturalFinalizedMultiCMGStatsArrays {
            ms2_stats,
            ms1_stats,
            id: self.id,
        }
    }
}
