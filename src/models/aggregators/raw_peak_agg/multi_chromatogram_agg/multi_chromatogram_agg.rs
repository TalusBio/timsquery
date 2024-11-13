use super::super::chromatogram_agg::ScanTofStatsCalculatorPair;
use super::aggregator::ParitionedCMGAggregator;
use super::arrays::{
    PartitionedCMGArrayStats,
    PartitionedCMGArrays,
};
use super::scored_arrays::PartitionedCMGScoredStatsArrays;
use crate::errors::{
    DataProcessingError,
    Result,
};
use crate::models::elution_group::ElutionGroup;
use crate::models::frames::raw_peak::RawPeak;
use crate::models::queries::MsLevelContext;
use crate::traits::aggregator::Aggregator;
use crate::utils::dtw::dtw_max;
use serde::Serialize;
use std::hash::Hash;

use timsrust::converters::{
    Scan2ImConverter,
    Tof2MzConverter,
};

#[derive(Debug, Clone)]
pub struct MultiCMGStatsAgg<FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug> {
    pub converters: (Tof2MzConverter, Scan2ImConverter),
    pub ms1_stats: ParitionedCMGAggregator<usize>,
    pub ms2_stats: ParitionedCMGAggregator<FH>,
    pub id: u64,
    pub context: Option<MsLevelContext<usize, FH>>,
    pub buffer: Option<ScanTofStatsCalculatorPair>,
    // TODO: Make this a reference instead of a clone.... maybe.
    pub elution_group_ref: ElutionGroup<FH>,
}

#[derive(Debug, Clone)]
pub struct MultiCMGStatsFactory<FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug> {
    pub converters: (Tof2MzConverter, Scan2ImConverter),
    pub _phantom: std::marker::PhantomData<FH>,
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug> MultiCMGStatsFactory<FH> {
    pub fn build_with_elution_group(
        &self,
        elution_group: &ElutionGroup<FH>,
    ) -> MultiCMGStatsAgg<FH> {
        // TODO: RN this is a super hacky ... IDEALLY we would either keep a reference OR
        // preserve only the expected intensities.
        let elution_group_ref = elution_group.clone();

        let fragment_keys = elution_group
            .fragment_mzs
            .keys()
            .cloned()
            .collect::<Vec<FH>>();

        let precursor_keys: Vec<usize> = elution_group
            .precursor_mzs
            .iter()
            .enumerate()
            .map(|x| x.0)
            .collect();

        let exp_fragment_intensities: Option<Vec<f32>> =
            if elution_group.expected_fragment_intensity.is_some() {
                let tmp_expect_inten = elution_group.expected_fragment_intensity.as_ref().unwrap();
                let out = fragment_keys
                    .iter()
                    .map(|x| *(tmp_expect_inten.get(x).unwrap_or(&0.0)))
                    .collect();
                Some(out)
            } else {
                None
            };

        MultiCMGStatsAgg {
            converters: (self.converters.0, self.converters.1),
            ms1_stats: ParitionedCMGAggregator::new(
                precursor_keys,
                elution_group.expected_precursor_intensity.to_owned(),
            ),
            ms2_stats: ParitionedCMGAggregator::new(fragment_keys, exp_fragment_intensities),
            id: elution_group.id,
            context: None,
            buffer: None,
            elution_group_ref,
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

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug>
    NaturalFinalizedMultiCMGStatsArrays<FH>
{
    pub fn finalized_score(&self) -> Result<Vec<f64>> {
        let ms1_cosine = match self.ms1_stats.scores.cosine_similarity.as_ref() {
            Some(x) => x,
            None => return Err(DataProcessingError::ExpectedNonEmptyData.into()),
        };
        let ms1_rt = self.ms1_stats.retention_time_miliseconds.as_slice();
        let ms2_cosine = match self.ms2_stats.scores.cosine_similarity.as_ref() {
            Some(x) => x,
            None => return Err(DataProcessingError::ExpectedNonEmptyData.into()),
        };
        let ms2_rt = self.ms2_stats.retention_time_miliseconds.as_slice();
        let ms2_lzb = self.ms2_stats.scores.lazyerscore_vs_baseline.as_slice();

        // Should I project MS1 -> MS2 or MS2 -> MS1?
        let proj_mz2_cosine = dtw_max(ms1_rt, ms2_rt, ms2_cosine)?;
        let proj_mz2_lzb = dtw_max(ms1_rt, ms2_rt, ms2_lzb)?;

        assert!(proj_mz2_cosine.len() == proj_mz2_lzb.len());
        assert!(proj_mz2_cosine.len() == ms1_rt.len());

        let mut out = Vec::with_capacity(proj_mz2_cosine.len());
        for (ms1_cos, (ms2_cos, ms2_lzb)) in ms1_cosine
            .iter()
            .zip(proj_mz2_cosine.iter().zip(proj_mz2_lzb.iter()))
        {
            let score = ms1_cos * ms2_cos * ms2_lzb;
            out.push(score);
        }

        Ok(out)
    }
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
            MsLevelContext::MS1(_i) => {
                self.ms1_stats.add(
                    rt_miliseconds,
                    peak.scan_index,
                    peak.tof_index,
                    u64_intensity,
                );
            }
            MsLevelContext::MS2(_i) => {
                self.ms2_stats.add(
                    rt_miliseconds,
                    peak.scan_index,
                    peak.tof_index,
                    u64_intensity,
                );
            }
        }
    }

    fn set_context(&mut self, context: Self::Context) {
        match &context {
            MsLevelContext::MS1(i) => {
                self.ms1_stats.set_context(*i).unwrap();
            }
            MsLevelContext::MS2(i) => {
                self.ms2_stats.set_context(i.clone()).unwrap();
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

        // TODO clean this up ... 4 conversions in one line is feeling a bit weird.
        let ms1_stats = PartitionedCMGScoredStatsArrays::new(
            PartitionedCMGArrayStats::from(PartitionedCMGArrays::from(self.ms1_stats)),
            mz_converter,
            mobility_converter,
        );
        let ms2_stats = PartitionedCMGScoredStatsArrays::new(
            PartitionedCMGArrayStats::from(PartitionedCMGArrays::from(self.ms2_stats)),
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
