use crate::models::adapters::FragmentIndexAdapter;
use crate::models::frames::raw_frames::frame_elems_matching;
use crate::models::frames::raw_peak::RawPeak;
use crate::models::queries::{FragmentGroupIndexQuery, NaturalPrecursorQuery, PrecursorIndexQuery};
use crate::traits::aggregator::Aggregator;
use crate::traits::queriable_data::QueriableData;
use crate::utils::tolerance_ranges::IncludedRange;
use crate::ElutionGroup;
use crate::ToleranceAdapter;
use rayon::iter::ParallelIterator;
use serde::Serialize;
use std::fmt::Debug;
use std::hash::Hash;
use timsrust::converters::ConvertableDomain;
use timsrust::readers::{FrameReader, FrameReaderError, MetadataReader};
use timsrust::TimsRustError;
use timsrust::{Frame, Metadata};
use tracing::trace;

pub struct RawFileIndex {
    file_reader: FrameReader,
    adapter: FragmentIndexAdapter,
}

impl RawFileIndex {
    pub fn from_path(path: &str) -> Result<Self, TimsRustError> {
        let file_reader = FrameReader::new(path)?;

        let sql_path = std::path::Path::new(path).join("analysis.tdf");
        let meta_converters = MetadataReader::new(&sql_path)?;
        let adapter = meta_converters.into();
        Ok(Self {
            file_reader,
            adapter,
        })
    }

    fn apply_on_query<
        'c,
        'b: 'c,
        'a: 'b,
        FH: Clone + Eq + Serialize + Hash + Debug + Send + Sync + Copy,
    >(
        &'a self,
        fqs: &'b FragmentGroupIndexQuery<FH>,
        fun: &'c mut dyn for<'r> FnMut(RawPeak, FH),
    ) {
        trace!("RawFileIndex::apply_on_query");
        trace!("FragmentGroupIndexQuery: {:?}", fqs);
        let frames: Vec<Frame> = self
            .read_frames_in_range(&fqs.precursor_query.frame_index_range.into())
            .filter(|x| x.is_ok())
            .map(|x| x.unwrap())
            .collect();

        let pq = &fqs.precursor_query;
        let iso_mz_range = pq.isolation_mz_range;
        for frame in frames {
            for (fh, tof_range) in fqs.mz_index_ranges.iter() {
                let scan_range = pq.mobility_index_range;
                let peaks = frame_elems_matching(
                    &frame,
                    *tof_range,
                    scan_range,
                    Some((iso_mz_range.start() as f64, iso_mz_range.end() as f64).into()),
                );
                for peak in peaks {
                    fun(peak, *fh);
                }
            }
        }
    }

    fn read_frames_in_range<'b, 'a: 'b>(
        &'a self,
        range: &'b (usize, usize),
    ) -> impl ParallelIterator<Item = Result<Frame, FrameReaderError>> + 'b {
        let lambda_use = |x: &Frame| x.index >= range.0 && x.index < range.1;
        self.file_reader.parallel_filter(lambda_use)
    }

    fn queries_from_elution_elements_impl<
        FH: Clone + Eq + Serialize + Hash + Send + Sync + Copy,
    >(
        &self,
        tol: &dyn crate::traits::tolerance::Tolerance,
        elution_elements: &crate::models::elution_group::ElutionGroup<FH>,
    ) -> FragmentGroupIndexQuery<FH> {
        self.adapter.query_from_elution_group(tol, elution_elements)
    }
}

impl<FH: Hash + Copy + Clone + Serialize + Eq + Debug + Send + Sync>
    QueriableData<FragmentGroupIndexQuery<FH>, (RawPeak, FH)> for RawFileIndex
{
    fn query(&self, fragment_query: &FragmentGroupIndexQuery<FH>) -> Vec<(RawPeak, FH)> {
        let mut out = Vec::new();
        self.apply_on_query(fragment_query, &mut |x, y| out.push((x, y)));
        out
    }

    fn add_query<A, O, AG>(&self, fragment_query: &FragmentGroupIndexQuery<FH>, aggregator: &mut AG)
    where
        A: From<(RawPeak, FH)> + Send + Sync + Clone + Copy,
        AG: Aggregator<Item = A, Output = O>,
    {
        self.apply_on_query(fragment_query, &mut |x, y| aggregator.add((x, y)));
    }

    fn add_query_multi_group<A, O, AG>(
        &self,
        fragment_queries: &[FragmentGroupIndexQuery<FH>],
        aggregator: &mut [AG],
    ) where
        A: From<(RawPeak, FH)> + Send + Sync + Clone + Copy,
        AG: Aggregator<Item = A, Output = O>,
    {
        fragment_queries
            .iter()
            .zip(aggregator.iter_mut())
            .for_each(|(query, agg)| self.add_query(query, agg));
    }
}

impl<FH: Hash + Serialize + Eq + Clone + Send + Sync + Copy>
    ToleranceAdapter<FragmentGroupIndexQuery<FH>, ElutionGroup<FH>> for RawFileIndex
{
    fn query_from_elution_group(
        &self,
        tol: &dyn crate::traits::tolerance::Tolerance,
        elution_group: &crate::models::elution_group::ElutionGroup<FH>,
    ) -> FragmentGroupIndexQuery<FH> {
        self.queries_from_elution_elements_impl(tol, elution_group)
    }
}
