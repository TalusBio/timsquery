use super::chromatogram_agg::{
    ChromatomobilogramStatsArrays, MappingCollection, ScanTofStatsCalculatorPair,
};
use crate::models::aggregators::rolling_calculators::rolling_median;
use crate::models::aggregators::streaming_aggregator::RunningStatsCalculator;
use crate::models::frames::raw_peak::RawPeak;
use crate::traits::aggregator::Aggregator;
use crate::utils::math::lnfact;
use serde::Serialize;
use std::collections::{BTreeMap, HashSet};
use std::hash::Hash;

use timsrust::converters::{ConvertableDomain, Scan2ImConverter, Tof2MzConverter};

#[derive(Debug, Clone)]
pub struct MultiCMGStats<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub scan_tof_mapping: MappingCollection<(FH, u32), ScanTofStatsCalculatorPair>,
    pub converters: (Tof2MzConverter, Scan2ImConverter),
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
    pub transition_mobilities: MappingCollection<FH, Vec<f64>>,
    pub transition_mzs: MappingCollection<FH, Vec<f64>>,
    pub transition_intensities: MappingCollection<FH, Vec<u64>>,
    pub apex_hyperscore_index: usize,
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
        let lazy_hyperscore_vs_baseline = calculate_value_vs_baseline(
            &lazy_hyperscore,
            1 + (other.retention_time_miliseconds.len() / 10),
        );
        let mut apex_hyperscore_index = 0;
        let mut max_hyperscore = 0.0f64;
        for (i, val) in lazy_hyperscore_vs_baseline.iter().enumerate() {
            if max_hyperscore.is_nan() || *val > max_hyperscore {
                max_hyperscore = *val;
                apex_hyperscore_index = i;
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
                .map(|x| mobility_converter.convert(x))
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
            apex_hyperscore_index,
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
        let mut summed_intensity_tree: BTreeMap<u32, u64> =
            BTreeMap::from_iter(unique_rts.iter().map(|x| (*x, 0)));
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

                tmp_tree.entry(rt).or_insert((scan, tof, inten));
                summed_intensity_tree
                    .entry(rt)
                    .and_modify(|curr| *curr += inten);
                npeaks_tree.entry(rt).and_modify(|curr| *curr += 1);
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
            .map(|x| x.mean().expect("At least one value should be present"))
            .collect();
        out
    }
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync> Aggregator<(RawPeak, FH)>
    for MultiCMGStats<FH>
{
    type Output = NaturalFinalizedMultiCMGStatsArrays<FH>;

    fn add(&mut self, peak: &(RawPeak, FH)) {
        let (peak, transition) = peak;
        let u64_intensity = peak.intensity as u64;
        let rt_miliseconds = (peak.retention_time * 1000.0) as u32;

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
                    curr.intensities.push(scan_tof_mapping.tof.weight());
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
                    out.intensities.push(scan_tof_mapping.tof.weight());
                    out
                });
        }

        for (_k, v) in transition_stats.iter_mut() {
            v.sort_by_rt();
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
