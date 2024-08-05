pub struct PrecursorIndexQuery {
    pub frame_index_range: (usize, usize),
    pub mz_index_range: (u32, u32),
    pub mobility_index_range: (usize, usize),
    pub isolation_mz_range: (f32, f32),
}

pub struct FragmentGroupIndexQuery {
    pub mz_index_ranges: Vec<(u32, u32)>,
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
