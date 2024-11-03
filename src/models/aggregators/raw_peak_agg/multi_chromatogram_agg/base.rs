use super::super::chromatogram_agg::{
    ChromatomobilogramStatsArrays,
    MappingCollection,
    ScanTofStatsCalculatorPair,
};
use crate::models::aggregators::rolling_calculators::{
    calculate_lazy_hyperscore,
    calculate_value_vs_baseline,
};
use crate::models::aggregators::streaming_aggregator::RunningStatsCalculator;
use crate::utils::math::lnfact_float;
use serde::Serialize;
use std::collections::{
    BTreeMap,
    HashSet,
};
use std::hash::Hash;
use timsrust::converters::{
    ConvertableDomain,
    Scan2ImConverter,
    Tof2MzConverter,
};
use tracing::{
    debug,
    warn,
};

#[derive(Debug, Clone)]
pub struct ParitionedCMGAggregator<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub scan_tof_mapping: MappingCollection<(FH, u32), ScanTofStatsCalculatorPair>,
    pub uniq_rts: HashSet<u32>,
    pub uniq_ids: HashSet<FH>,
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

// This name is starting to get really long ...
#[derive(Debug, Clone, Serialize)]
pub struct PartitionedCMGScoredStatsArrays<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
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

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync> Default for ParitionedCMGAggregator<FH> {
    fn default() -> Self {
        Self {
            scan_tof_mapping: MappingCollection::new(),
            uniq_rts: HashSet::new(),
            uniq_ids: HashSet::new(),
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct PartitionedCMGArrays<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub transition_stats: MappingCollection<FH, ChromatomobilogramStatsArrays>,
}

// TODO reimplement ... dont really like how non-idiomatic this finalize is.
// From-Into might be a better way.
impl<FH: Clone + Eq + Serialize + Hash + Send + Sync> ParitionedCMGAggregator<FH> {
    pub fn finalize(self) -> PartitionedCMGArrays<FH> {
        let mut transition_stats = MappingCollection::new();

        for id_key in self.uniq_ids.iter() {
            let mut id_cmgs = ChromatomobilogramStatsArrays::new();
            for rt_key in self.uniq_rts.iter() {
                let scan_tof_mapping = self.scan_tof_mapping.get(&(id_key.clone(), *rt_key));
                if let Some(scan_tof_mapping) = scan_tof_mapping {
                    id_cmgs.retention_time_miliseconds.push(*rt_key);
                    id_cmgs
                        .scan_index_means
                        .push(scan_tof_mapping.scan.mean().unwrap());
                    id_cmgs
                        .scan_index_sds
                        .push(scan_tof_mapping.scan.standard_deviation().unwrap());
                    id_cmgs
                        .tof_index_means
                        .push(scan_tof_mapping.tof.mean().unwrap());
                    id_cmgs
                        .tof_index_sds
                        .push(scan_tof_mapping.tof.standard_deviation().unwrap());
                    id_cmgs.intensities.push(scan_tof_mapping.tof.weight());
                }
            }
            id_cmgs.sort_by_rt();
            transition_stats.insert(id_key.clone(), id_cmgs);
        }

        PartitionedCMGArrays { transition_stats }
    }
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync> PartitionedCMGScoredStatsArrays<FH> {
    // TODO: Refactor this function ... its getting pretty big.
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

        PartitionedCMGScoredStatsArrays {
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
