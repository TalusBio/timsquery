use super::arrays::PartitionedCMGArrayStats;
use crate::models::aggregators::rolling_calculators::{
    calculate_lazy_hyperscore,
    calculate_value_vs_baseline,
};
use crate::models::aggregators::streaming_aggregator::RunningStatsCalculator;
use crate::utils::math::lnfact_float;
use serde::Serialize;
use std::collections::HashMap;
use std::f64;
use std::hash::Hash;
use timsrust::converters::{
    ConvertableDomain,
    Scan2ImConverter,
    Tof2MzConverter,
};
use tracing::debug;

// This name is starting to get really long ...
#[derive(Debug, Clone, Serialize)]
pub struct PartitionedCMGScoredStatsArrays<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub retention_time_miliseconds: Vec<u32>,
    pub average_mobility: Vec<f64>,
    pub summed_intensity: Vec<u64>,
    pub scores: ChromatogramScores,
    pub transition_mobilities: HashMap<FH, Vec<f64>>,
    pub transition_mzs: HashMap<FH, Vec<f64>>,
    pub transition_intensities: HashMap<FH, Vec<u64>>,
    pub apex_primary_score_index: usize,
}

#[derive(Debug, Clone, Serialize)]
pub struct ChromatogramScores {
    pub lazy_hyperscore: Vec<f64>,
    pub lazyerscore: Vec<f64>,
    pub lazy_hyperscore_vs_baseline: Vec<f64>,
    pub lazyerscore_vs_baseline: Vec<f64>,
    pub norm_hyperscore_vs_baseline: Vec<f64>,
    pub norm_lazyerscore_vs_baseline: Vec<f64>,
    pub npeaks: Vec<u8>,
}

impl ChromatogramScores {
    pub fn new<T: Clone + Eq + Serialize + Hash + Send + Sync>(
        arrays: &PartitionedCMGArrayStats<T>,
    ) -> Self {
        let lazy_hyperscore = calculate_lazy_hyperscore(&arrays.npeaks, &arrays.summed_intensity);
        let basline_window_len = 1 + (arrays.retention_time_miliseconds.len() / 20);
        let lazy_hyperscore_vs_baseline =
            calculate_value_vs_baseline(&lazy_hyperscore, basline_window_len);
        let lazyerscore: Vec<f64> = arrays
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

        // Make sure all arrays are the same length
        assert!(
            lazy_hyperscore_vs_baseline.len() == arrays.retention_time_miliseconds.len(),
            "Failed sanity check"
        );
        assert!(
            lazyerscore_vs_baseline.len() == arrays.retention_time_miliseconds.len(),
            "Failed sanity check"
        );

        ChromatogramScores {
            lazy_hyperscore,
            lazyerscore,
            lazy_hyperscore_vs_baseline,
            lazyerscore_vs_baseline,
            norm_hyperscore_vs_baseline,
            norm_lazyerscore_vs_baseline,
            npeaks: arrays.npeaks.to_owned(),
        }
    }

    fn get_apex_score_index(&self) -> usize {
        let mut apex_primary_score_index = 0;
        let mut max_primary_score = 0.0f64;
        let primary_scores = &self.norm_lazyerscore_vs_baseline;
        for (i, val) in primary_scores.iter().enumerate() {
            if max_primary_score.is_nan() || *val > max_primary_score {
                max_primary_score = *val;
                apex_primary_score_index = i;
            }
        }
        apex_primary_score_index
    }
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync> PartitionedCMGScoredStatsArrays<FH> {
    pub fn new(
        other: PartitionedCMGArrayStats<FH>,
        mz_converter: &Tof2MzConverter,
        mobility_converter: &Scan2ImConverter,
    ) -> Self {
        let scores = ChromatogramScores::new(&other);
        let apex_primary_score_index = scores.get_apex_score_index();
        let average_mobility = other
            .weighted_scan_index_mean
            .into_iter()
            .map(|x| {
                let out = mobility_converter.convert(x);
                if !(0.5..=2.0).contains(&out) {
                    debug!("Bad mobility value: {:?}, input was {:?}", out, x);
                }
                out
            })
            .collect();

        let transition_mobilities = other
            .scan_index_means
            .into_iter()
            .map(|(k, v)| {
                let new_v = v
                    .into_iter()
                    .map(|x| mobility_converter.convert(x))
                    .collect();
                (k, new_v)
            })
            .collect();

        let transition_mzs = other
            .tof_index_means
            .into_iter()
            .map(|(k, v)| {
                let new_v = v.into_iter().map(|x| mz_converter.convert(x)).collect();
                (k, new_v)
            })
            .collect();

        PartitionedCMGScoredStatsArrays {
            retention_time_miliseconds: other.retention_time_miliseconds,
            average_mobility,
            summed_intensity: other.summed_intensity,
            transition_mobilities,
            transition_mzs,
            transition_intensities: other.intensities,
            apex_primary_score_index,
            scores,
        }
    }
}
