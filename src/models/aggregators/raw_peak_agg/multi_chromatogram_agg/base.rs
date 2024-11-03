use super::super::chromatogram_agg::{
    ChromatomobilogramStatsArrays,
    MappingCollection,
    ScanTofStatsCalculatorPair,
};
use crate::{
    models::{
        aggregators::{
            rolling_calculators::rolling_median,
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

#[derive(Debug, Clone)]
pub struct ParitionedCMGAggregator<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub scan_tof_mapping: MappingCollection<(FH, u32), ScanTofStatsCalculatorPair>,
    pub uniq_rts: HashSet<u32>,
    pub uniq_ids: HashSet<FH>,
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
