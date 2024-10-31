use crate::traits::aggregator::ProvidesContext;
use crate::utils::tolerance_ranges::IncludedRange;
use std::collections::HashMap;
use std::hash::Hash;
use timsrust::converters::{ConvertableDomain, Scan2ImConverter, Tof2MzConverter};

#[derive(Debug, Clone)]
pub struct PrecursorIndexQuery {
    pub frame_index_range: IncludedRange<usize>,
    pub rt_range_seconds: IncludedRange<f32>,
    pub mz_index_ranges: Vec<IncludedRange<u32>>,
    pub mobility_index_range: IncludedRange<usize>,
    pub isolation_mz_range: IncludedRange<f32>,
}

#[derive(Debug, Clone)]
pub enum MsLevelContext<T1, T2> {
    MS1(T1),
    MS2(T2),
}

#[derive(Debug, Clone)]
pub struct FragmentGroupIndexQuery<FH: Clone + Eq + Hash + Send + Sync + Copy> {
    pub mz_index_ranges: HashMap<FH, IncludedRange<u32>>,
    pub precursor_query: PrecursorIndexQuery,
}

impl<FH: Clone + Eq + Hash + Send + Sync + Copy> ProvidesContext for FragmentGroupIndexQuery<FH> {
    type Context = MsLevelContext<FH, usize>;
}

#[derive(Debug, Clone)]
pub struct NaturalPrecursorQuery {
    pub rt_range: IncludedRange<f32>,
    pub mobility_range: IncludedRange<f32>,
    pub mz_ranges: Vec<IncludedRange<f64>>,
    pub isolation_mz_range: IncludedRange<f32>,
}

impl NaturalPrecursorQuery {
    pub fn as_precursor_query(
        &self,
        mz_converter: &Tof2MzConverter,
        im_converter: &Scan2ImConverter,
    ) -> PrecursorIndexQuery {
        PrecursorIndexQuery {
            frame_index_range: (
                im_converter.invert(self.mobility_range.start()).round() as usize,
                im_converter.invert(self.mobility_range.end()).round() as usize,
            )
                .into(),
            rt_range_seconds: self.rt_range.clone(),
            mz_index_ranges: self
                .mz_ranges
                .iter()
                .map(|mz_range| {
                    (
                        mz_converter.invert(mz_range.start()).round() as u32,
                        mz_converter.invert(mz_range.end()).round() as u32,
                    )
                        .into()
                })
                .collect(),
            mobility_index_range: (
                im_converter.invert(self.mobility_range.start()).round() as usize,
                im_converter.invert(self.mobility_range.end()).round() as usize,
            )
                .into(),
            isolation_mz_range: self.isolation_mz_range.clone(),
        }
    }
}
