use super::{
    super::chromatogram_agg::{
        ChromatomobilogramStatsArrays,
        MappingCollection,
        ScanTofStatsCalculatorPair,
    },
    base::{
        ParitionedCMGAggregator,
        PartitionedCMGArrays,
    },
};
use crate::{
    models::{
        aggregators::{
            rolling_calculators::{
                calculate_lazy_hyperscore,
                calculate_value_vs_baseline,
            },
            streaming_aggregator::RunningStatsCalculator,
        },
        frames::raw_peak::RawPeak,
        queries::MsLevelContext,
    },
    traits::aggregator::Aggregator,
    utils::math::{
        lnfact,
        lnfact_float,
    },
};
use serde::Serialize;
use std::{
    collections::{
        BTreeMap,
        HashSet,
    },
    hash::Hash,
};
use tracing::{
    debug,
    warn,
};

use timsrust::converters::{
    ConvertableDomain,
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
pub struct PartitionedCMGArrayStats<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub retention_time_miliseconds: Vec<u32>,
    pub weighted_scan_index_mean: Vec<f64>,
    pub summed_intensity: Vec<u64>,

    // This should be the same as the sum of log intensities.
    pub log_intensity_products: Vec<f64>,
    pub npeaks: Vec<usize>,
    pub scan_index_means: MappingCollection<FH, Vec<f64>>,
    pub tof_index_means: MappingCollection<FH, Vec<f64>>,
    // TODO consider if I want to add the standard deviations ... RN they dont
    // seem to be that useful.
    pub intensities: MappingCollection<FH, Vec<u64>>,
}

#[derive(Debug, Clone, Serialize)]
pub struct FinalizedMultiCMGStatsArrays<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub ms1_stats: PartitionedCMGArrayStats<usize>,
    pub ms2_stats: PartitionedCMGArrayStats<FH>,
    pub id: u64,
}

// This name is starting to get really long ...
#[derive(Debug, Clone, Serialize)]
pub struct _NaturalFinalizedMultiCMGStatsArrays<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub retention_time_miliseconds: Vec<u32>,
    pub average_mobility: Vec<f64>,
    pub summed_intensity: Vec<u64>,
    pub npeaks: Vec<usize>,
    pub lazy_hyperscore: Vec<f64>,
    pub lazy_hyperscore_vs_baseline: Vec<f64>,
    pub norm_hyperscore_vs_baseline: Vec<f64>,
    pub lazyerscore: Vec<f64>,
    pub lazyerscore_vs_baseline: Vec<f64>,
    pub norm_lazyerscore_vs_baseline: Vec<f64>,
    pub transition_mobilities: MappingCollection<FH, Vec<f64>>,
    pub transition_mzs: MappingCollection<FH, Vec<f64>>,
    pub transition_intensities: MappingCollection<FH, Vec<u64>>,
    pub apex_primary_score_index: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct NaturalFinalizedMultiCMGStatsArrays<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub ms1_stats: _NaturalFinalizedMultiCMGStatsArrays<usize>,
    pub ms2_stats: _NaturalFinalizedMultiCMGStatsArrays<FH>,
    pub id: u64,
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync> _NaturalFinalizedMultiCMGStatsArrays<FH> {
    pub fn new(
        other: PartitionedCMGArrayStats<FH>,
        mz_converter: &Tof2MzConverter,
        mobility_converter: &Scan2ImConverter,
    ) -> Self {
        let lazy_hyperscore = calculate_lazy_hyperscore(&other.npeaks, &other.summed_intensity);
        let basline_window_len = 1 + (other.retention_time_miliseconds.len() / 10);
        let lazy_hyperscore_vs_baseline =
            calculate_value_vs_baseline(&lazy_hyperscore, basline_window_len);
        let lazyerscore: Vec<f64> = other
            .log_intensity_products
            .iter()
            .map(|x| lnfact_float(*x))
            .collect();
        let lazyerscore_vs_baseline = calculate_value_vs_baseline(&lazyerscore, basline_window_len);

        // Set 0 the NANs
        // Q: Can I do this in-place? Will the comnpiler do it for me?
        let lazy_hyperscore_vs_baseline: Vec<f64> = lazy_hyperscore_vs_baseline
            .into_iter()
            .map(|x| if x.is_nan() { 0.0 } else { x })
            .collect();
        let lazyerscore_vs_baseline: Vec<f64> = lazyerscore_vs_baseline
            .into_iter()
            .map(|x| if x.is_nan() { 0.0 } else { x })
            .collect();

        // Calculate the standard deviation of the lazyscores v baseline
        let mut sd_calculator_hyperscore = RunningStatsCalculator::new(1, 0.0);
        let mut sd_calculator_lazyscore = RunningStatsCalculator::new(1, 0.0);

        (0..lazyerscore_vs_baseline.len()).for_each(|i| {
            sd_calculator_hyperscore.add(lazy_hyperscore_vs_baseline[i], 1);
            sd_calculator_lazyscore.add(lazyerscore_vs_baseline[i], 1);
        });
        let sd_hyperscore = sd_calculator_hyperscore.standard_deviation().unwrap();
        let sd_lazyscore = sd_calculator_lazyscore.standard_deviation().unwrap();

        // Calculate the normalized scores
        let norm_hyperscore_vs_baseline: Vec<f64> = lazy_hyperscore_vs_baseline
            .iter()
            .map(|x| x / sd_hyperscore)
            .collect();
        let norm_lazyerscore_vs_baseline: Vec<f64> = lazyerscore_vs_baseline
            .iter()
            .map(|x| x / sd_lazyscore)
            .collect();

        let mut apex_primary_score_index = 0;
        let mut max_primary_score = 0.0f64;
        let primary_scores = &norm_lazyerscore_vs_baseline;
        for (i, val) in primary_scores.iter().enumerate() {
            if max_primary_score.is_nan() || *val > max_primary_score {
                max_primary_score = *val;
                apex_primary_score_index = i;
            }
        }

        assert!(
            lazy_hyperscore_vs_baseline.len() == other.retention_time_miliseconds.len(),
            "Failed sanity check"
        );

        _NaturalFinalizedMultiCMGStatsArrays {
            retention_time_miliseconds: other.retention_time_miliseconds,
            average_mobility: other
                .weighted_scan_index_mean
                .into_iter()
                .map(|x| {
                    let out = mobility_converter.convert(x);
                    if !(0.5..=2.0).contains(&out) {
                        debug!("Bad mobility value: {:?}, input was {:?}", out, x);
                    }
                    out
                })
                .collect(),
            summed_intensity: other.summed_intensity,
            npeaks: other.npeaks,
            transition_mobilities: other
                .scan_index_means
                .into_iter()
                .map(|(k, v)| {
                    let new_v = v
                        .into_iter()
                        .map(|x| mobility_converter.convert(x))
                        .collect();
                    (k, new_v)
                })
                .collect(),
            transition_mzs: other
                .tof_index_means
                .into_iter()
                .map(|(k, v)| {
                    let new_v = v.into_iter().map(|x| mz_converter.convert(x)).collect();
                    (k, new_v)
                })
                .collect(),
            transition_intensities: other.intensities,
            lazy_hyperscore,
            lazy_hyperscore_vs_baseline,
            apex_primary_score_index,
            lazyerscore,
            lazyerscore_vs_baseline,
            norm_hyperscore_vs_baseline,
            norm_lazyerscore_vs_baseline,
        }
    }
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync> From<PartitionedCMGArrays<FH>>
    for PartitionedCMGArrayStats<FH>
{
    fn from(other: PartitionedCMGArrays<FH>) -> Self {
        // TODO ... maybe refactor this ... RN its king of ugly.

        let mut out = PartitionedCMGArrayStats {
            retention_time_miliseconds: Vec::new(),
            scan_index_means: MappingCollection::new(),
            tof_index_means: MappingCollection::new(),
            weighted_scan_index_mean: Vec::new(),
            intensities: MappingCollection::new(),
            summed_intensity: Vec::new(),
            log_intensity_products: Vec::new(),
            npeaks: Vec::new(),
        };

        let unique_rts = other
            .transition_stats
            .iter()
            .flat_map(|(_k, v)| v.retention_time_miliseconds.clone())
            .collect::<HashSet<u32>>();
        let mut unique_rts = unique_rts.into_iter().collect::<Vec<u32>>();
        unique_rts.sort();

        // Q: Is this the most efficient way to do this?
        // I think having the btrees as the unit of integration might not be the best idea.
        // ... If I want to preserve the sparsity, I can use a hashmap and the sort it.
        let mut summed_intensity_tree: BTreeMap<u32, u64> =
            BTreeMap::from_iter(unique_rts.iter().map(|x| (*x, 0)));
        let mut intensity_logsums_tree: BTreeMap<u32, f64> =
            BTreeMap::from_iter(unique_rts.iter().map(|x| (*x, 0.0)));
        let mut npeaks_tree: BTreeMap<u32, usize> =
            BTreeMap::from_iter(unique_rts.iter().map(|x| (*x, 0)));
        let mut weighted_tof_index_mean_tree: BTreeMap<u32, RunningStatsCalculator> =
            BTreeMap::new();
        let mut weighted_scan_index_mean_tree: BTreeMap<u32, RunningStatsCalculator> =
            BTreeMap::new();

        for (k, v) in other.transition_stats.into_iter() {
            type ScanTofIntensityTuple = (f64, f64, u64);
            type ScanTofIntensityVecs = (Vec<f64>, (Vec<f64>, Vec<u64>));
            let mut tmp_tree: BTreeMap<u32, ScanTofIntensityTuple> = BTreeMap::new();

            for i in 0..v.retention_time_miliseconds.len() {
                let rt = v.retention_time_miliseconds[i];
                let scan = v.scan_index_means[i];
                let tof = v.tof_index_means[i];
                let inten = v.intensities[i];

                if inten > 100 {
                    npeaks_tree.entry(rt).and_modify(|curr| *curr += 1);
                }
                tmp_tree.entry(rt).or_insert((scan, tof, inten));
                summed_intensity_tree
                    .entry(rt)
                    .and_modify(|curr| *curr += inten);

                intensity_logsums_tree.entry(rt).and_modify(|curr| {
                    if inten > 10 {
                        *curr += (inten as f64).ln()
                    }
                });
                weighted_tof_index_mean_tree
                    .entry(rt)
                    .and_modify(|curr| curr.add(tof, inten))
                    .or_insert(RunningStatsCalculator::new(inten, tof));
                weighted_scan_index_mean_tree
                    .entry(rt)
                    .and_modify(|curr| curr.add(scan, inten))
                    .or_insert(RunningStatsCalculator::new(inten, scan));
            }

            // Now we fill with nans the missing values.
            // Q: Do I really need to do this here?
            for rt in unique_rts.iter() {
                tmp_tree.entry(*rt).or_insert((f64::NAN, f64::NAN, 0));
            }

            let (out_scans, (out_tofs, out_intens)): ScanTofIntensityVecs = tmp_tree
                .into_iter()
                .map(|(_, (scan, tof, inten))| (scan, (tof, inten)))
                .unzip();
            out.scan_index_means.insert(k.clone(), out_scans);
            out.tof_index_means.insert(k.clone(), out_tofs);
            out.intensities.insert(k.clone(), out_intens);
        }

        out.retention_time_miliseconds = unique_rts;
        out.summed_intensity = summed_intensity_tree.into_values().collect();
        out.npeaks = npeaks_tree.into_values().collect();
        out.weighted_scan_index_mean = weighted_scan_index_mean_tree
            .into_values()
            .map(|x| {
                let out = x.mean().expect("At least one value should be present");
                if !(0.0..=1000.0).contains(&out) {
                    warn!("Bad mobility value: {:?}, input was {:?}", out, x);
                }

                out
            })
            .collect();

        // Note: The log of products is the same as the sum of logs.
        out.log_intensity_products = intensity_logsums_tree.into_values().collect();
        out
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

        let ms1_stats = _NaturalFinalizedMultiCMGStatsArrays::new(
            PartitionedCMGArrayStats::from(self.ms1_stats.finalize()),
            mz_converter,
            mobility_converter,
        );
        let ms2_stats = _NaturalFinalizedMultiCMGStatsArrays::new(
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
