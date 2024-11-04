use super::super::chromatogram_agg::ChromatomobilogramStatsArrays;
use super::aggregator::ParitionedCMGAggregator;
use crate::models::aggregators::streaming_aggregator::RunningStatsCalculator;
use serde::Serialize;
use std::collections::{
    BTreeMap,
    HashMap,
    HashSet,
};
use std::f64;
use std::hash::Hash;
use tracing::warn;

#[derive(Debug, Clone, Serialize)]
pub struct PartitionedCMGArrayStats<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub retention_time_miliseconds: Vec<u32>,
    pub weighted_scan_index_mean: Vec<f64>,
    pub summed_intensity: Vec<u64>,

    // This should be the same as the sum of log intensities.
    pub log_intensity_products: Vec<f64>,
    pub npeaks: Vec<usize>,
    pub scan_index_means: HashMap<FH, Vec<f64>>,
    pub tof_index_means: HashMap<FH, Vec<f64>>,
    // TODO consider if I want to add the standard deviations ... RN they dont
    // seem to be that useful.
    pub intensities: HashMap<FH, Vec<u64>>,
}

#[derive(Debug, Clone, Serialize)]
pub struct PartitionedCMGArrays<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub transition_stats: HashMap<FH, ChromatomobilogramStatsArrays>,
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug>
    From<ParitionedCMGAggregator<FH>> for PartitionedCMGArrays<FH>
{
    fn from(mut item: ParitionedCMGAggregator<FH>) -> Self {
        item.flush_buffer();
        let mut transition_stats = HashMap::new();
        let mut uniq_rts: Vec<u32> = item
            .scan_tof_calc
            .iter()
            .flatten()
            .flatten()
            .map(|x| *x.0)
            .collect();
        uniq_rts.sort_unstable();
        uniq_rts.dedup();

        for (id_ind, id_key) in item.keys.iter().enumerate() {
            let mut id_cmgs = ChromatomobilogramStatsArrays::new();
            let local_id_mapping = &item.scan_tof_calc[id_ind];
            if local_id_mapping.is_none() {
                continue;
            }
            let local_id_mapping = local_id_mapping.as_ref().unwrap();

            for rt_key in uniq_rts.iter() {
                let scan_tof_mapping = local_id_mapping.get(&rt_key);
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

        Self { transition_stats }
    }
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync> From<PartitionedCMGArrays<FH>>
    for PartitionedCMGArrayStats<FH>
{
    fn from(other: PartitionedCMGArrays<FH>) -> Self {
        // TODO ... maybe refactor this ... RN its king of ugly.

        let mut out = PartitionedCMGArrayStats {
            retention_time_miliseconds: Vec::new(),
            scan_index_means: HashMap::new(),
            tof_index_means: HashMap::new(),
            weighted_scan_index_mean: Vec::new(),
            intensities: HashMap::new(),
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
