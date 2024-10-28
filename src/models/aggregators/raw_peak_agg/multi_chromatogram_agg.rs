use super::chromatogram_agg::{
    ChromatomobilogramStatsArrays, MappingCollection, ScanTofStatsCalculatorPair,
};
use crate::models::aggregators::rolling_calculators::rolling_median;
use crate::models::aggregators::streaming_aggregator::RunningStatsCalculator;
use crate::models::frames::raw_peak::RawPeak;
use crate::traits::aggregator::Aggregator;
use crate::utils::math::{lnfact, lnfact_float};
use serde::Serialize;
use std::collections::{BTreeMap, HashSet};
use std::hash::Hash;
use tracing::{debug, warn};

use timsrust::converters::{ConvertableDomain, Scan2ImConverter, Tof2MzConverter};

#[derive(Debug, Clone)]
pub struct MultiCMGStats<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub scan_tof_mapping: MappingCollection<(FH, u32), ScanTofStatsCalculatorPair>,
    pub converters: (Tof2MzConverter, Scan2ImConverter),
    pub uniq_rts: HashSet<u32>,
    pub uniq_ids: HashSet<FH>,
    pub id: u64,
}

#[derive(Debug, Clone)]
pub struct MultiCMGStatsFactory<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub converters: (Tof2MzConverter, Scan2ImConverter),
    pub _phantom: std::marker::PhantomData<FH>,
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync> MultiCMGStatsFactory<FH> {
    pub fn build(&self, id: u64) -> MultiCMGStats<FH> {
        MultiCMGStats {
            scan_tof_mapping: MappingCollection::new(),
            converters: (self.converters.0, self.converters.1),
            uniq_rts: HashSet::new(),
            uniq_ids: HashSet::new(),
            id,
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct MultiCMGStatsArrays<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub transition_stats: MappingCollection<FH, ChromatomobilogramStatsArrays>,
    id: u64,
}

#[derive(Debug, Clone, Serialize)]
pub struct FinalizedMultiCMGStatsArrays<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
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
    pub id: u64,
}

// This name is starting to get really long ...
#[derive(Debug, Clone, Serialize)]
pub struct NaturalFinalizedMultiCMGStatsArrays<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
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
    pub id: u64,
}

// Reference python code

// log1p_intensitties = np.log1p(arrays["summed_intensity"])
//     lazy_hyperscore = LN_FACTORIALS[arrays["npeaks"]] + (2 * log1p_intensitties)
//
//     five_pct_len = int(len(arrays["retention_time_miliseconds"]) / 100) * 5
//     lazy_hyperscore_roll_median = (
//         pd.Series(lazy_hyperscore).rolling(window=five_pct_len, center=True).median()
// ).array
//     hyperscore_off_baseline = lazy_hyperscore - lazy_hyperscore_roll_median
//     hyperscore_off_baseline = np.nan_to_num(hyperscore_off_baseline, 0)
//         iqr_hyperscore_off_baseline = np.percentile(
//             hyperscore_off_baseline, 75
//     ) - np.percentile(hyperscore_off_baseline, 25)
//
//     scaled_hyperscore_off_baseline = hyperscore_off_baseline / iqr_hyperscore_off_baseline
//     )
// )

fn calculate_lazy_hyperscore(npeaks: &[usize], summed_intensity: &[u64]) -> Vec<f64> {
    let mut scores = vec![0.0; npeaks.len()];
    for i in 0..npeaks.len() {
        let npeaks_i = npeaks[i];
        let summed_intensity_i = summed_intensity[i];
        let log1p_intensities_i = (summed_intensity_i as f64 + 1.0).ln();
        scores[i] = lnfact(npeaks_i as u16) + (2.0 * log1p_intensities_i);
    }
    scores
}

fn calculate_value_vs_baseline(vals: &[f64], baseline_window_size: usize) -> Vec<f64> {
    let baseline = rolling_median(vals, baseline_window_size, f64::NAN);
    vals.iter()
        .zip(baseline.iter())
        .map(|(x, y)| x - y)
        .collect()
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync> NaturalFinalizedMultiCMGStatsArrays<FH> {
    pub fn new(
        other: FinalizedMultiCMGStatsArrays<FH>,
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

        NaturalFinalizedMultiCMGStatsArrays {
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
            id: other.id,
        }
    }
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync> From<MultiCMGStatsArrays<FH>>
    for FinalizedMultiCMGStatsArrays<FH>
{
    fn from(other: MultiCMGStatsArrays<FH>) -> Self {
        // TODO ... maybe refactor this ... RN its king of ugly.

        let mut out = FinalizedMultiCMGStatsArrays {
            retention_time_miliseconds: Vec::new(),
            scan_index_means: MappingCollection::new(),
            tof_index_means: MappingCollection::new(),
            weighted_scan_index_mean: Vec::new(),
            intensities: MappingCollection::new(),
            summed_intensity: Vec::new(),
            log_intensity_products: Vec::new(),
            npeaks: Vec::new(),
            id: other.id,
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

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync> Aggregator for MultiCMGStats<FH> {
    type Item = (RawPeak, FH);
    type Output = NaturalFinalizedMultiCMGStatsArrays<FH>;

    fn add(&mut self, peak: &(RawPeak, FH)) {
        let (peak, transition) = peak;
        let u64_intensity = peak.intensity as u64;
        let rt_miliseconds = (peak.retention_time * 1000.0) as u32;

        self.uniq_rts.insert(rt_miliseconds);
        self.uniq_ids.insert(transition.clone());

        self.scan_tof_mapping
            .entry((transition.clone(), rt_miliseconds))
            .and_modify(|curr| {
                curr.add(u64_intensity, peak.scan_index, peak.tof_index);
            })
            .or_insert(ScanTofStatsCalculatorPair::new(
                u64_intensity,
                peak.scan_index,
                peak.tof_index,
            ));
    }

    fn finalize(self) -> NaturalFinalizedMultiCMGStatsArrays<FH> {
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

        let tmp = MultiCMGStatsArrays {
            transition_stats,
            id: self.id,
        };
        NaturalFinalizedMultiCMGStatsArrays::new(
            FinalizedMultiCMGStatsArrays::from(tmp),
            &self.converters.0,
            &self.converters.1,
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_value_vs_baseline() {
        let vals = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0];
        let baseline_window_size = 3;
        let _baseline = rolling_median(&vals, baseline_window_size, f64::NAN);
        let out = calculate_value_vs_baseline(&vals, baseline_window_size);
        let expect_val = vec![f64::NAN, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, f64::NAN];
        let all_close = out
            .iter()
            .zip(expect_val.iter())
            .filter(|(a, b)| ((!a.is_nan()) && (!b.is_nan())))
            .all(|(a, b)| (a - b).abs() < 1e-6);

        let all_match_nan = out
            .iter()
            .zip(expect_val.iter())
            .filter(|(a, b)| ((a.is_nan()) || (b.is_nan())))
            .all(|(a, b)| a.is_nan() && b.is_nan());

        assert!(all_close, "Expected {:?}, got {:?}", expect_val, out);
        assert!(all_match_nan, "Expected {:?}, got {:?}", expect_val, out);
    }
}
