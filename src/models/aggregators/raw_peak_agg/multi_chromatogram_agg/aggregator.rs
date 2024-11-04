use super::super::chromatogram_agg::{
    ChromatomobilogramStatsArrays,
    ScanTofStatsCalculatorPair,
};
use crate::errors::{
    DataProcessingError,
    Result,
};
use crate::utils::correlation::rolling_cosine_similarity;
use nohash_hasher::BuildNoHashHasher;
use serde::Serialize;
use std::collections::{
    HashMap,
    HashSet,
};
use std::f64;
use std::hash::Hash;
use std::sync::Arc;

type SparseRTCollection = HashMap<u32, ScanTofStatsCalculatorPair, BuildNoHashHasher<u32>>;

// pub uniq_rts: Arc<[u32]>,
// Move upstream
#[derive(Debug, Clone)]
pub struct ParitionedCMGAggregator<
    FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug,
> {
    pub scan_tof_calc: Vec<Option<SparseRTCollection>>,
    pub keys: Vec<FH>,
    pub context_key_num: usize,
    pub expected_intensities: Option<Vec<f32>>,
    pub context_buffer: SparseRTCollection,
}

impl<FH: Clone + Eq + Serialize + Hash + Send + Sync + std::fmt::Debug>
    ParitionedCMGAggregator<FH>
{
    pub fn new(keys: Vec<FH>) -> Self {
        let mut scan_tof_calc = Vec::with_capacity(keys.len());
        for _ in 0..keys.len() {
            scan_tof_calc.push(None);
        }
        Self {
            scan_tof_calc,
            keys,
            context_key_num: 0,
            expected_intensities: None,
            context_buffer: SparseRTCollection::with_hasher(BuildNoHashHasher::default()),
        }
    }

    pub fn flush_buffer(&mut self) {
        if !self.context_buffer.is_empty() {
            // Swap the buffer with the current one.
            let new_buffer = SparseRTCollection::with_hasher(BuildNoHashHasher::default());
            self.scan_tof_calc[self.context_key_num] = Some(new_buffer);
            std::mem::swap(
                &mut self.context_buffer,
                self.scan_tof_calc[self.context_key_num].as_mut().unwrap(),
            );
        }
    }

    pub fn set_context(&mut self, context: FH) -> Result<()> {
        self.flush_buffer();

        // Find the position in the keys.
        let pos = self.keys.iter().position(|x| x == &context);
        if let Some(pos) = pos {
            self.context_key_num = pos;
            Ok(())
        } else {
            let msg = format!(
                "Context Not Found, wante any of {:?}, got {:?}",
                self.keys, context
            );
            Err(crate::TimsqueryError::Other(msg))
        }
    }

    pub fn add(&mut self, rt_ms: u32, scan_index: usize, tof_index: u32, intensity: u64) {
        self.context_buffer
            .entry(rt_ms)
            .and_modify(|curr| {
                curr.add(intensity, scan_index, tof_index);
            })
            .or_insert(ScanTofStatsCalculatorPair::new(
                intensity, scan_index, tof_index,
            ));
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct PartitionedCMGArrays<FH: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub transition_stats: HashMap<FH, ChromatomobilogramStatsArrays>,
}

fn avg_correlation_vs_hashmap<'a, T>(
    ref_arr: &'a [f64],
    hm: &'a HashMap<T, Vec<f64>>,
    window_size: usize,
) -> Result<Vec<f64>>
where
    T: Hash + Eq + Clone,
{
    let mut num_added = 0;
    let mut added = vec![0.0; ref_arr.len()];
    for target in hm.values() {
        let res = rolling_cosine_similarity(ref_arr, target, window_size);
        if res.is_err() {
            continue;
        }
        let res = res.unwrap();
        for (i, val) in res.iter().enumerate() {
            if val.is_nan() {
                continue;
            }
            added[i] += val;
        }
        num_added += 1;
    }

    if num_added == 0 {
        return Err(DataProcessingError::ExpectedNonEmptyData.into());
    }

    let added_f64 = num_added as f64;
    for x in added.iter_mut() {
        *x /= added_f64;
    }
    Ok(added)
}
