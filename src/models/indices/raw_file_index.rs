use core::panic;
use std::process::Output;
use std::sync::Arc;

use rayon::iter::ParallelIterator;
use timsrust::converters::ConvertableDomain;
use timsrust::readers::{
    FrameReader, FrameReaderError, MetadataReader, QuadrupoleSettingsReader,
};
use timsrust::TimsRustError;
use timsrust::{Frame, Metadata, QuadrupoleSettings};
// use timsrust::io::;

use crate::models::frames::raw_frames::frame_elems_matching;
use crate::models::frames::raw_peak::RawPeak;
use crate::models::queries::{
    FragmentGroupIndexQuery, NaturalFragmentQuery, NaturalPrecursorQuery, PrecursorIndexQuery,
};
use crate::traits::indexed_data::IndexedData;
use crate::ToleranceAdapter;

pub struct RawFileIndex {
    file_reader: FrameReader,
    meta_converters: Metadata,
    quad_settings_reader: Vec<QuadrupoleSettings>,
}

impl RawFileIndex {
    pub fn from_path(path: &str) -> Result<Self, TimsRustError> {
        let file_reader = FrameReader::new(&path)?;

        let sql_path = std::path::Path::new(path).join("analysis.tdf");
        let meta_converters = MetadataReader::new(&sql_path)?;
        let quad_settings_reader = QuadrupoleSettingsReader::new(&sql_path)?;
        Ok(Self {
            file_reader,
            meta_converters,
            quad_settings_reader,
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

    fn apply_on_query<'c, 'b: 'c, 'a: 'b>(
        &'a self,
        fqs: &'b FragmentGroupIndexQuery,
        fun: &'c mut dyn for<'r> FnMut(RawPeak, Arc<QuadrupoleSettings>),
    ) {
        let frames: Vec<Frame> = self
            .read_frames_in_range(&fqs.precursor_query.frame_index_range)
            .filter(|x| x.is_ok())
            .map(|x| x.unwrap())
            .collect();

        let pq = &fqs.precursor_query;
        let iso_mz_range = pq.isolation_mz_range;
        for frame in frames {
            for tof_range in fqs.mz_index_ranges.iter() {
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

    fn query_from_elution_group_impl(
        &self,
        tol: &dyn crate::traits::tolerance::Tolerance,
        elution_group: &crate::models::elution_group::ElutionGroup,
    ) -> PrecursorIndexQuery {
        let rt_range = tol.rt_range(elution_group.rt_seconds);
        let mobility_range = tol.mobility_range(elution_group.mobility);
        let mz_range = tol.mz_range(elution_group.precursor_mz);
        let quad_range = tol.quad_range(elution_group.precursor_mz, elution_group.precursor_charge);

        PrecursorIndexQuery {
            frame_index_range: (
                self.meta_converters.rt_converter.invert(rt_range.0) as usize,
                self.meta_converters.rt_converter.invert(rt_range.1) as usize,
            ),
            mz_index_range: (
                self.meta_converters.mz_converter.invert(mz_range.0) as u32,
                self.meta_converters.mz_converter.invert(mz_range.1) as u32,
            ),
            mobility_index_range: (
                self.meta_converters.im_converter.invert(mobility_range.0) as usize,
                self.meta_converters.im_converter.invert(mobility_range.1) as usize,
            ),
            isolation_mz_range: quad_range,
        }
    }

    fn queries_from_elution_elements_impl<'d>(
        &self,
        tol: &dyn crate::traits::tolerance::Tolerance,
        elution_elements: &crate::models::elution_group::ElutionGroup,
    ) -> FragmentGroupIndexQuery {
        let precursor_query = self.query_from_elution_group_impl(tol, elution_elements);
        let fragment_mzs = elution_elements.fragment_mzs.as_ref().unwrap();
        let mut fqs = Vec::with_capacity(fragment_mzs.len());
        for fragment in fragment_mzs {
            let mz_range = tol.mz_range(*fragment);
            fqs.push((
                self.meta_converters.mz_converter.invert(mz_range.0) as u32,
                self.meta_converters.mz_converter.invert(mz_range.1) as u32,
            ));
        }

        FragmentGroupIndexQuery {
            mz_index_ranges: fqs,
            precursor_query,
        }
    }
}

impl IndexedData<FragmentGroupIndexQuery, RawPeak> for RawFileIndex {
    fn query(&self, fragment_query: &FragmentGroupIndexQuery) -> Option<Vec<RawPeak>> {
        let mut out = Vec::new();
        self.apply_on_query(&fragment_query, &mut |peak, _| out.push(peak));
        Some(out)
    }

    fn add_query<O, AG: crate::Aggregator<RawPeak, Output = O>>(
        &self,
        fragment_query: &FragmentGroupIndexQuery,
        aggregator: &mut AG,
    ) {
        self.apply_on_query(&fragment_query, &mut |peak, _| aggregator.add(&peak));
    }

    fn add_query_multi_group<O, AG: crate::Aggregator<RawPeak, Output = O>>(
        &self,
        fragment_queries: &[FragmentGroupIndexQuery],
        aggregator: &mut [AG],
    ) {
        fragment_queries
            .iter()
            .zip(aggregator.iter_mut())
            .for_each(|(query, agg)| self.add_query(query, agg));
    }
}

impl ToleranceAdapter<FragmentGroupIndexQuery> for RawFileIndex {
    fn query_from_elution_group(
        &self,
        tol: &dyn crate::traits::tolerance::Tolerance,
        elution_group: &crate::models::elution_group::ElutionGroup,
    ) -> FragmentGroupIndexQuery {
        self.queries_from_elution_elements_impl(tol, elution_group)
    }
}

pub struct RawPeakIntensityAggregator {
    pub intensity: u64,
}

impl crate::traits::aggregator::Aggregator<RawPeak> for RawPeakIntensityAggregator {
    type Output = u64;

    fn add(&mut self, peak: &RawPeak) {
        self.intensity += peak.intensity as u64;
    }

    fn fold(&mut self, other: Self) {
        self.intensity += other.intensity;
    }

    fn finalize(&self) -> u64 {
        self.intensity
    }
}
