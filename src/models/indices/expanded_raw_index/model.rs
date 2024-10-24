use crate::models::adapters::FragmentIndexAdapter;
use crate::models::elution_group::ElutionGroup;
use crate::models::frames::expanded_frame::ExpandedFrame;
use crate::models::frames::expanded_frame::{
    expand_and_arrange_frames, expand_and_split_frame, par_expand_and_centroid_frames,
    ExpandedFrameSlice, SortedState,
};
use crate::models::frames::peak_in_quad::PeakInQuad;
use crate::models::frames::raw_peak::RawPeak;
use crate::models::frames::single_quad_settings::{
    get_matching_quad_settings, SingleQuadrupoleSetting, SingleQuadrupoleSettingIndex,
};
use crate::models::queries::FragmentGroupIndexQuery;
use crate::traits::indexed_data::IndexedData;
use crate::ToleranceAdapter;
use log::{debug, info};
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashMap;
use std::hash::Hash;
use std::time::Instant;
use timsrust::converters::{
    ConvertableDomain, Frame2RtConverter, Scan2ImConverter, Tof2MzConverter,
};
use timsrust::readers::{FrameReader, FrameReaderError, MetadataReader};
use timsrust::{QuadrupoleSettings, TimsRustError};

type QuadSettingsIndex = usize;

#[derive(Debug)]
pub struct ExpandedRawFrameIndex {
    bundled_ms1_frames: ExpandedSliceBundle,
    bundled_frames: HashMap<SingleQuadrupoleSettingIndex, ExpandedSliceBundle>,
    flat_quad_settings: Vec<SingleQuadrupoleSetting>,
    rt_converter: Frame2RtConverter,
    pub mz_converter: Tof2MzConverter,
    pub im_converter: Scan2ImConverter,
    adapter: FragmentIndexAdapter,
}

#[derive(Debug, Clone)]
pub struct ExpandedSliceBundle {
    slices: Vec<ExpandedFrameSlice<SortedState>>,
    rts: Vec<f64>,
    frame_indices: Vec<usize>,
}

impl ExpandedSliceBundle {
    pub fn new(mut slices: Vec<ExpandedFrameSlice<SortedState>>) -> Self {
        slices.sort_unstable_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap());
        let rts = slices.iter().map(|x| x.rt).collect();
        let frame_indices = slices.iter().map(|x| x.frame_index).collect();
        Self {
            slices,
            rts,
            frame_indices,
        }
    }

    pub fn query_peaks<F>(
        &self,
        tof_range: (u32, u32),
        scan_range: Option<(usize, usize)>,
        frame_index_range: (usize, usize),
        f: &mut F,
    ) where
        F: FnMut(PeakInQuad),
    {
        // Binary search the rt if needed.
        let frame_indices = self.frame_indices.as_slice();
        let low = frame_indices.partition_point(|x| x < &frame_index_range.0);
        let high = frame_indices.partition_point(|x| x <= &frame_index_range.1);

        for i in low..high {
            let slice = &self.slices[i];
            slice.query_peaks(tof_range, scan_range, f);
        }
    }
}

impl ExpandedRawFrameIndex {
    pub fn query_peaks<F>(
        &self,
        tof_range: (u32, u32),
        precursor_mz_range: (f64, f64),
        scan_range: Option<(usize, usize)>,
        frame_index_range: (usize, usize),
        f: &mut F,
    ) where
        F: FnMut(PeakInQuad),
    {
        let matching_quads: Vec<SingleQuadrupoleSettingIndex> =
            get_matching_quad_settings(&self.flat_quad_settings, precursor_mz_range, scan_range)
                .collect();
        self.query_precursor_peaks(&matching_quads, tof_range, scan_range, frame_index_range, f);
    }

    fn query_precursor_peaks<F>(
        &self,
        matching_quads: &[SingleQuadrupoleSettingIndex],
        tof_range: (u32, u32),
        scan_range: Option<(usize, usize)>,
        frame_index_range: (usize, usize),
        f: &mut F,
    ) where
        F: FnMut(PeakInQuad),
    {
        for quad in matching_quads {
            let tqi = self
                .bundled_frames
                .get(quad)
                .expect("Only existing quads should be queried.");

            tqi.query_peaks(tof_range, scan_range, frame_index_range, &mut *f)
        }
    }

    pub fn from_path_centroided(path: &str) -> Result<Self, TimsRustError> {
        let file_reader = FrameReader::new(path)?;

        let sql_path = std::path::Path::new(path).join("analysis.tdf");
        let meta_converters = MetadataReader::new(&sql_path)?;

        let st = Instant::now();
        let all_frames = file_reader.get_all().into_iter().flatten().collect();
        let read_elap = st.elapsed();
        info!("Reading all frames took {:#?}", read_elap);
        let st = Instant::now();
        // TODO: Expose this parameter as a config option.
        let centroided_split_frames = par_expand_and_centroid_frames(
            all_frames,
            1.5,
            15.0,
            &meta_converters.im_converter,
            &meta_converters.mz_converter,
        );
        let centroided_elap = st.elapsed();
        info!("Centroiding took {:#?}", centroided_elap);

        let mut out_ms2_frames = HashMap::new();
        let mut out_ms1_frames: Option<ExpandedSliceBundle> = None;

        let mut flat_quad_settings = Vec::new();
        centroided_split_frames
            .into_iter()
            .for_each(|(q, frameslices)| match q {
                None => {
                    out_ms1_frames = Some(ExpandedSliceBundle::new(frameslices));
                }
                Some(q) => {
                    flat_quad_settings.push(q);
                    out_ms2_frames.insert(q.index, ExpandedSliceBundle::new(frameslices));
                }
            });

        let adapter = FragmentIndexAdapter::from(meta_converters.clone());

        let out = Self {
            bundled_ms1_frames: out_ms1_frames.expect("At least one ms1 frame should be present"),
            bundled_frames: out_ms2_frames,
            flat_quad_settings,
            rt_converter: meta_converters.rt_converter,
            mz_converter: meta_converters.mz_converter,
            im_converter: meta_converters.im_converter,
            adapter,
        };

        Ok(out)
    }

    pub fn from_path(path: &str) -> Result<Self, TimsRustError> {
        // NOTE: I am just copy-pasting the centroided version. If I keep both I will
        // abstract it and make dispatch in a config ...

        let file_reader = FrameReader::new(path)?;

        let sql_path = std::path::Path::new(path).join("analysis.tdf");
        let meta_converters = MetadataReader::new(&sql_path)?;

        let st = Instant::now();
        let all_frames = file_reader.get_all().into_iter().flatten().collect();
        let read_elap = st.elapsed();
        info!("Reading all frames took {:#?}", read_elap);
        let st = Instant::now();
        let centroided_split_frames = expand_and_arrange_frames(all_frames);
        let centroided_elap = st.elapsed();
        info!("Splitting took {:#?}", centroided_elap);

        let mut out_ms2_frames = HashMap::new();
        let mut out_ms1_frames: Option<ExpandedSliceBundle> = None;

        let mut flat_quad_settings = Vec::new();
        centroided_split_frames
            .into_iter()
            .for_each(|(q, frameslices)| match q {
                None => {
                    out_ms1_frames = Some(ExpandedSliceBundle::new(frameslices));
                }
                Some(q) => {
                    flat_quad_settings.push(q);
                    out_ms2_frames.insert(q.index, ExpandedSliceBundle::new(frameslices));
                }
            });

        let adapter = FragmentIndexAdapter::from(meta_converters.clone());

        let out = Self {
            bundled_ms1_frames: out_ms1_frames.expect("At least one ms1 frame should be present"),
            bundled_frames: out_ms2_frames,
            flat_quad_settings,
            rt_converter: meta_converters.rt_converter,
            mz_converter: meta_converters.mz_converter,
            im_converter: meta_converters.im_converter,
            adapter,
        };

        Ok(out)
    }
}

impl<FH: Eq + Hash + Copy + Serialize + Send + Sync>
    IndexedData<FragmentGroupIndexQuery<FH>, (RawPeak, FH)> for ExpandedRawFrameIndex
{
    fn query(&self, fragment_query: &FragmentGroupIndexQuery<FH>) -> Vec<(RawPeak, FH)> {
        todo!();
        // let precursor_mz_range = (
        //     fragment_query.precursor_query.isolation_mz_range.0 as f64,
        //     fragment_query.precursor_query.isolation_mz_range.0 as f64,
        // );
        // let scan_range = Some(fragment_query.precursor_query.mobility_index_range);
        // let frame_index_range = Some(FrameRTTolerance::FrameIndex(
        //     fragment_query.precursor_query.frame_index_range,
        // ));

        // fragment_query
        //     .mz_index_ranges
        //     .iter()
        //     .flat_map(|(fh, tof_range)| {
        //         let mut local_vec: Vec<(RawPeak, FH)> = vec![];
        //         self.query_peaks(
        //             *tof_range,
        //             precursor_mz_range,
        //             scan_range,
        //             frame_index_range,
        //             &mut |x| local_vec.push((RawPeak::from(x), *fh)),
        //         );

        //         local_vec
        //     })
        //     .collect()
    }

    fn add_query<O, AG: crate::Aggregator<(RawPeak, FH), Output = O>>(
        &self,
        fragment_query: &FragmentGroupIndexQuery<FH>,
        aggregator: &mut AG,
    ) {
        todo!();
        // let precursor_mz_range = (
        //     fragment_query.precursor_query.isolation_mz_range.0 as f64,
        //     fragment_query.precursor_query.isolation_mz_range.0 as f64,
        // );
        // let scan_range = Some(fragment_query.precursor_query.mobility_index_range);
        // let frame_index_range = Some(FrameRTTolerance::FrameIndex(
        //     fragment_query.precursor_query.frame_index_range,
        // ));

        // fragment_query
        //     .mz_index_ranges
        //     .iter()
        //     .for_each(|(fh, tof_range)| {
        //         self.query_peaks(
        //             *tof_range,
        //             precursor_mz_range,
        //             scan_range,
        //             frame_index_range,
        //             &mut |peak| aggregator.add(&(RawPeak::from(peak), *fh)),
        //         );
        //     })
    }

    fn add_query_multi_group<O, AG: crate::Aggregator<(RawPeak, FH), Output = O>>(
        &self,
        fragment_queries: &[FragmentGroupIndexQuery<FH>],
        aggregator: &mut [AG],
    ) {
        fragment_queries
            .par_iter()
            .zip(aggregator.par_iter_mut())
            .for_each(|(fragment_query, agg)| {
                let precursor_mz_range = (
                    fragment_query.precursor_query.isolation_mz_range.0 as f64,
                    fragment_query.precursor_query.isolation_mz_range.1 as f64,
                );
                assert!(precursor_mz_range.0 <= precursor_mz_range.1);
                assert!(precursor_mz_range.0 > 0.0);
                let scan_range = Some(fragment_query.precursor_query.mobility_index_range);
                let frame_index_range = fragment_query.precursor_query.frame_index_range;

                let local_quad_vec: Vec<SingleQuadrupoleSettingIndex> = get_matching_quad_settings(
                    &self.flat_quad_settings,
                    precursor_mz_range,
                    scan_range,
                )
                .collect();

                for (fh, tof_range) in fragment_query.mz_index_ranges.clone().into_iter() {
                    self.query_precursor_peaks(
                        &local_quad_vec,
                        tof_range,
                        scan_range,
                        frame_index_range,
                        &mut |peak| agg.add(&(RawPeak::from(peak), fh)),
                    );
                }
            });
    }
}

impl<FH: Eq + Hash + Copy + Serialize + Send + Sync>
    IndexedData<FragmentGroupIndexQuery<FH>, RawPeak> for ExpandedRawFrameIndex
{
    fn query(&self, fragment_query: &FragmentGroupIndexQuery<FH>) -> Vec<RawPeak> {
        todo!();
        // let precursor_mz_range = (
        //     fragment_query.precursor_query.isolation_mz_range.0 as f64,
        //     fragment_query.precursor_query.isolation_mz_range.0 as f64,
        // );
        // let scan_range = Some(fragment_query.precursor_query.mobility_index_range);
        // let frame_index_range = Some(FrameRTTolerance::FrameIndex(
        //     fragment_query.precursor_query.frame_index_range,
        // ));

        // fragment_query
        //     .mz_index_ranges
        //     .iter()
        //     .flat_map(|(fh, tof_range)| {
        //         let mut local_vec: Vec<(RawPeak, FH)> = vec![];
        //         self.query_peaks(
        //             *tof_range,
        //             precursor_mz_range,
        //             scan_range,
        //             frame_index_range,
        //             &mut |x| local_vec.push((RawPeak::from(x), *fh)),
        //         );

        //         local_vec
        //     })
        //     .collect()
    }

    fn add_query<O, AG: crate::Aggregator<RawPeak, Output = O>>(
        &self,
        fragment_query: &FragmentGroupIndexQuery<FH>,
        aggregator: &mut AG,
    ) {
        todo!();
        // let precursor_mz_range = (
        //     fragment_query.precursor_query.isolation_mz_range.0 as f64,
        //     fragment_query.precursor_query.isolation_mz_range.0 as f64,
        // );
        // let scan_range = Some(fragment_query.precursor_query.mobility_index_range);
        // let frame_index_range = Some(FrameRTTolerance::FrameIndex(
        //     fragment_query.precursor_query.frame_index_range,
        // ));

        // fragment_query
        //     .mz_index_ranges
        //     .iter()
        //     .for_each(|(fh, tof_range)| {
        //         self.query_peaks(
        //             *tof_range,
        //             precursor_mz_range,
        //             scan_range,
        //             frame_index_range,
        //             &mut |peak| aggregator.add(&(RawPeak::from(peak), *fh)),
        //         );
        //     })
    }

    fn add_query_multi_group<O, AG: crate::Aggregator<RawPeak, Output = O>>(
        &self,
        fragment_queries: &[FragmentGroupIndexQuery<FH>],
        aggregator: &mut [AG],
    ) {
        fragment_queries
            .par_iter()
            .zip(aggregator.par_iter_mut())
            .for_each(|(fragment_query, agg)| {
                let precursor_mz_range = (
                    fragment_query.precursor_query.isolation_mz_range.0 as f64,
                    fragment_query.precursor_query.isolation_mz_range.1 as f64,
                );
                assert!(precursor_mz_range.0 <= precursor_mz_range.1);
                assert!(precursor_mz_range.0 > 0.0);
                let scan_range = Some(fragment_query.precursor_query.mobility_index_range);
                let frame_index_range = fragment_query.precursor_query.frame_index_range;

                let local_quad_vec: Vec<SingleQuadrupoleSettingIndex> = get_matching_quad_settings(
                    &self.flat_quad_settings,
                    precursor_mz_range,
                    scan_range,
                )
                .collect();

                for (fh, tof_range) in fragment_query.mz_index_ranges.clone().into_iter() {
                    self.query_precursor_peaks(
                        &local_quad_vec,
                        tof_range,
                        scan_range,
                        frame_index_range,
                        &mut |peak| agg.add(&RawPeak::from(peak)),
                    );
                }
            });
    }
}

// ============================================================================

impl<FH: Copy + Clone + Serialize + Eq + Hash + Send + Sync>
    ToleranceAdapter<FragmentGroupIndexQuery<FH>, ElutionGroup<FH>> for ExpandedRawFrameIndex
{
    fn query_from_elution_group(
        &self,
        tol: &dyn crate::traits::tolerance::Tolerance,
        elution_group: &ElutionGroup<FH>,
    ) -> FragmentGroupIndexQuery<FH> {
        self.adapter.query_from_elution_group(tol, elution_group)
    }
}
