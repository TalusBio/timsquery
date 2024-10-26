use std::sync::Arc;

use crate::models::frames::raw_frames::frame_elems_matching;
use crate::models::frames::raw_peak::RawPeak;
use crate::models::queries::{FragmentGroupIndexQuery, NaturalPrecursorQuery, PrecursorIndexQuery};
use crate::traits::aggregator::Aggregator;
use crate::traits::indexed_data::QueriableData;
use crate::ElutionGroup;
use crate::ToleranceAdapter;
use log::trace;
use rayon::iter::ParallelIterator;
use serde::Serialize;
use std::fmt::Debug;
use std::hash::Hash;
use timsrust::converters::ConvertableDomain;
use timsrust::readers::{FrameReader, FrameReaderError, MetadataReader};
use timsrust::TimsRustError;
use timsrust::{Frame, Metadata, QuadrupoleSettings};

pub struct RawFileIndex {
    file_reader: FrameReader,
    meta_converters: Metadata,
}

impl RawFileIndex {
    pub fn from_path(path: &str) -> Result<Self, TimsRustError> {
        let file_reader = FrameReader::new(path)?;

        let sql_path = std::path::Path::new(path).join("analysis.tdf");
        let meta_converters = MetadataReader::new(&sql_path)?;
        Ok(Self {
            file_reader,
            meta_converters,
        })
    }

    pub fn convert_precursor_query(&self, query: &NaturalPrecursorQuery) -> PrecursorIndexQuery {
        PrecursorIndexQuery {
            frame_index_range: (
                self.meta_converters.rt_converter.invert(query.rt_range.0) as usize,
                self.meta_converters.rt_converter.invert(query.rt_range.1) as usize,
            ),
            mz_index_range: (
                self.meta_converters.mz_converter.invert(query.mz_range.0) as u32,
                self.meta_converters.mz_converter.invert(query.mz_range.1) as u32,
            ),
            mobility_index_range: (
                self.meta_converters
                    .im_converter
                    .invert(query.mobility_range.0) as usize,
                self.meta_converters
                    .im_converter
                    .invert(query.mobility_range.1) as usize,
            ),
            isolation_mz_range: query.isolation_mz_range,
        }
    }

    fn apply_on_query<
        'c,
        'b: 'c,
        'a: 'b,
        FH: Clone + Eq + Serialize + Hash + Debug + Send + Sync + Copy,
    >(
        &'a self,
        fqs: &'b FragmentGroupIndexQuery<FH>,
        fun: &'c mut dyn for<'r> FnMut(RawPeak, Arc<QuadrupoleSettings>),
    ) {
        trace!("RawFileIndex::apply_on_query");
        trace!("FragmentGroupIndexQuery: {:?}", fqs);
        let frames: Vec<Frame> = self
            .read_frames_in_range(&fqs.precursor_query.frame_index_range)
            .filter(|x| x.is_ok())
            .map(|x| x.unwrap())
            .collect();

        let pq = &fqs.precursor_query;
        let iso_mz_range = pq.isolation_mz_range;
        for frame in frames {
            for (_, tof_range) in fqs.mz_index_ranges.iter() {
                let scan_range = pq.mobility_index_range;
                let peaks = frame_elems_matching(
                    &frame,
                    *tof_range,
                    scan_range,
                    Some((iso_mz_range.0 as f64, iso_mz_range.1 as f64)),
                );
                for peak in peaks {
                    fun(peak, frame.quadrupole_settings.clone());
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

    fn query_from_elution_group_impl<FH: Clone + Eq + Serialize + Hash + Send + Sync>(
        &self,
        tol: &dyn crate::traits::tolerance::Tolerance,
        elution_group: &crate::models::elution_group::ElutionGroup<FH>,
    ) -> PrecursorIndexQuery {
        let rt_range = tol.rt_range(elution_group.rt_seconds);
        let mz_range = tol.mz_range(elution_group.precursor_mz);
        let quad_range = tol.quad_range(elution_group.precursor_mz, elution_group.precursor_charge);

        let mobility_range = tol.mobility_range(elution_group.mobility);
        let mobility_range = match mobility_range {
            Some(mobility_range) => mobility_range,
            None => (
                self.meta_converters.lower_im as f32,
                self.meta_converters.upper_im as f32,
            ),
        };
        let mut min_scan_index =
            self.meta_converters.im_converter.invert(mobility_range.0) as usize;
        let mut max_scan_index =
            self.meta_converters.im_converter.invert(mobility_range.1) as usize;

        if min_scan_index > max_scan_index {
            std::mem::swap(&mut min_scan_index, &mut max_scan_index);
        }

        let rt_range = match rt_range {
            Some(rt_range) => rt_range,
            None => (
                self.meta_converters.lower_rt as f32,
                self.meta_converters.upper_rt as f32,
            ),
        };

        PrecursorIndexQuery {
            frame_index_range: (
                self.meta_converters.rt_converter.invert(rt_range.0) as usize,
                self.meta_converters.rt_converter.invert(rt_range.1) as usize,
            ),
            mz_index_range: (
                self.meta_converters.mz_converter.invert(mz_range.0) as u32,
                self.meta_converters.mz_converter.invert(mz_range.1) as u32,
            ),
            mobility_index_range: (min_scan_index, max_scan_index),
            isolation_mz_range: quad_range,
        }
    }

    fn queries_from_elution_elements_impl<
        FH: Clone + Eq + Serialize + Hash + Send + Sync + Copy,
    >(
        &self,
        tol: &dyn crate::traits::tolerance::Tolerance,
        elution_elements: &crate::models::elution_group::ElutionGroup<FH>,
    ) -> FragmentGroupIndexQuery<FH> {
        let precursor_query = self.query_from_elution_group_impl(tol, elution_elements);
        // TODO: change this unwrap and use explicitly the lack of fragment mzs.
        // Does that mean its onlyt a precursor query?
        // Why is it an option?
        //
        // let fragment_mzs = elution_elements.fragment_mzs.as_ref().unwrap();

        let fqs = elution_elements
            .fragment_mzs
            .iter()
            .map(|(k, v)| {
                let mz_range = tol.mz_range(*v);
                (
                    k.clone(),
                    (
                        self.meta_converters.mz_converter.invert(mz_range.0) as u32,
                        self.meta_converters.mz_converter.invert(mz_range.1) as u32,
                    ),
                )
            })
            .collect();

        FragmentGroupIndexQuery {
            mz_index_ranges: fqs,
            precursor_query,
        }
    }
}

impl<FH: Hash + Copy + Clone + Serialize + Eq + Debug + Send + Sync>
    QueriableData<FragmentGroupIndexQuery<FH>, RawPeak> for RawFileIndex
{
    fn query(&self, fragment_query: &FragmentGroupIndexQuery<FH>) -> Vec<RawPeak> {
        let mut out = Vec::new();
        self.apply_on_query(fragment_query, &mut |peak, _| out.push(peak));
        out
    }

    fn add_query<O, AG: Aggregator<Item = RawPeak, Output = O>>(
        &self,
        fragment_query: &FragmentGroupIndexQuery<FH>,
        aggregator: &mut AG,
    ) {
        self.apply_on_query(fragment_query, &mut |peak, _| aggregator.add(&peak));
    }

    fn add_query_multi_group<O, AG: crate::Aggregator<Item = RawPeak, Output = O>>(
        &self,
        fragment_queries: &[FragmentGroupIndexQuery<FH>],
        aggregator: &mut [AG],
    ) {
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
