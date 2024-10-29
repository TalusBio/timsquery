use std::collections::HashMap;
use std::hash::Hash;

#[derive(Debug, Clone)]
pub struct PrecursorIndexQuery {
    pub frame_index_range: (usize, usize),
    pub rt_range_seconds: (f32, f32),
    pub mz_index_range: (u32, u32),
    pub mobility_index_range: (usize, usize),
    pub isolation_mz_range: (f32, f32),
}

#[derive(Debug, Clone)]
pub struct FragmentGroupIndexQuery<FH: Clone + Eq + Hash + Send + Sync + Copy> {
    pub mz_index_ranges: HashMap<FH, (u32, u32)>,
    pub precursor_query: PrecursorIndexQuery,
}

pub struct NaturalPrecursorQuery {
    pub rt_range: (f32, f32),
    pub mobility_range: (f32, f32),
    pub mz_range: (f64, f64),
    pub isolation_mz_range: (f32, f32),
}

pub struct NaturalFragmentQuery<'a> {
    pub mz_range: (f64, f64),
    pub precursor_query: &'a NaturalPrecursorQuery,
}
