use crate::errors::Result;
use crate::models::adapters::FragmentIndexAdapter;
use crate::models::elution_group::ElutionGroup;
use crate::models::frames::expanded_frame::{
    par_read_and_expand_frames, ExpandedFrameSlice, FrameProcessingConfig, SortedState,
};
use crate::models::frames::peak_in_quad::PeakInQuad;
use crate::models::frames::raw_peak::RawPeak;
use crate::models::frames::single_quad_settings::{
    get_matching_quad_settings, matches_quad_settings, SingleQuadrupoleSetting,
    SingleQuadrupoleSettingIndex,
};
use crate::models::queries::FragmentGroupIndexQuery;
use crate::traits::aggregator::Aggregator;
use crate::traits::queriable_data::QueriableData;
use crate::ToleranceAdapter;
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashMap;
use std::hash::Hash;
use std::time::Instant;
use timsrust::converters::{Scan2ImConverter, Tof2MzConverter};
use timsrust::readers::{FrameReader, MetadataReader};
use tracing::info;
use tracing::instrument;

#[derive(Debug)]
pub struct ExpandedRawFrameIndex {
    bundled_ms1_frames: ExpandedSliceBundle,
    bundled_frames: HashMap<SingleQuadrupoleSettingIndex, ExpandedSliceBundle>,
    flat_quad_settings: Vec<SingleQuadrupoleSetting>,
    pub mz_converter: Tof2MzConverter,
    pub im_converter: Scan2ImConverter,
    adapter: FragmentIndexAdapter,
}

#[derive(Debug, Clone)]
pub struct ExpandedSliceBundle {
    slices: Vec<ExpandedFrameSlice<SortedState>>,
    frame_indices: Vec<usize>,
}

impl ExpandedSliceBundle {
    pub fn new(mut slices: Vec<ExpandedFrameSlice<SortedState>>) -> Self {
        slices.sort_unstable_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap());
        let frame_indices = slices.iter().map(|x| x.frame_index).collect();
        Self {
            slices,
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

    #[instrument(name = "ExpandedRawFrameIndex::from_path_centroided")]
    pub fn from_path_centroided(path: &str) -> Result<Self> {
        let config = FrameProcessingConfig::default_centroided();
        Self::from_path_base(path, config)
    }

    #[instrument(name = "ExpandedRawFrameIndex::from_path")]
    pub fn from_path(path: &str) -> Result<Self> {
        Self::from_path_base(path, FrameProcessingConfig::NotCentroided)
    }

    #[instrument(name = "ExpandedRawFrameIndex::from_path_base")]
    pub fn from_path_base(path: &str, centroid_config: FrameProcessingConfig) -> Result<Self> {
        info!(
            "Building ExpandedRawFrameIndex from path {} config {:?}",
            path, centroid_config,
        );
        let file_reader = FrameReader::new(path)?;

        let sql_path = std::path::Path::new(path).join("analysis.tdf");
        let meta_converters = MetadataReader::new(&sql_path)?;
        let centroid_config = centroid_config
            .with_converters(meta_converters.im_converter, meta_converters.mz_converter);

        let st = Instant::now();
        let centroided_split_frames = par_read_and_expand_frames(&file_reader, centroid_config)?;
        let centroided_elap = st.elapsed();
        info!("Reading + Centroiding took {:#?}", centroided_elap);

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
            mz_converter: meta_converters.mz_converter,
            im_converter: meta_converters.im_converter,
            adapter,
        };

        Ok(out)
    }
}

impl<FH: Eq + Hash + Copy + Serialize + Send + Sync>
    QueriableData<FragmentGroupIndexQuery<FH>, (RawPeak, FH)> for ExpandedRawFrameIndex
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

    fn add_query<A, O, AG>(&self, fragment_query: &FragmentGroupIndexQuery<FH>, aggregator: &mut AG)
    where
        A: From<(RawPeak, FH)> + Send + Sync + Clone + Copy,
        AG: Aggregator<Item = A, Output = O>,
    {
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

    fn add_query_multi_group<A, O, AG>(
        &self,
        fragment_queries: &[FragmentGroupIndexQuery<FH>],
        aggregator: &mut [AG],
    ) where
        A: From<(RawPeak, FH)> + Send + Sync + Clone + Copy,
        AG: Aggregator<Item = A, Output = O>,
    {
        // fragment_queries
        //     .par_iter()
        //     .zip(aggregator.par_iter_mut())
        //     .for_each(|(fragment_query, agg)| {
        //         let precursor_mz_range = (
        //             fragment_query.precursor_query.isolation_mz_range.0 as f64,
        //             fragment_query.precursor_query.isolation_mz_range.1 as f64,
        //         );
        //         assert!(precursor_mz_range.0 <= precursor_mz_range.1);
        //         assert!(precursor_mz_range.0 > 0.0);
        //         let scan_range = Some(fragment_query.precursor_query.mobility_index_range);
        //         let frame_index_range = fragment_query.precursor_query.frame_index_range;

        //         let local_quad_vec: Vec<SingleQuadrupoleSettingIndex> = get_matching_quad_settings(
        //             &self.flat_quad_settings,
        //             precursor_mz_range,
        //             scan_range,
        //         )
        //         .collect();

        //         for (fh, tof_range) in fragment_query.mz_index_ranges.clone().into_iter() {
        //             self.query_precursor_peaks(
        //                 &local_quad_vec,
        //                 tof_range,
        //                 scan_range,
        //                 frame_index_range,
        //                 &mut |peak| agg.add(&(RawPeak::from(peak), fh)),
        //             );
        //         }
        //     });
        let prec_mz_ranges = fragment_queries
            .iter()
            .map(|x| {
                (
                    x.precursor_query.isolation_mz_range.0 as f64,
                    x.precursor_query.isolation_mz_range.0 as f64,
                )
            })
            .collect::<Vec<_>>();

        let scan_ranges = fragment_queries
            .iter()
            .map(|x| Some(x.precursor_query.mobility_index_range))
            .collect::<Vec<_>>();

        let frame_index_ranges = fragment_queries
            .iter()
            .map(|x| x.precursor_query.frame_index_range)
            .collect::<Vec<_>>();

        for quad_setting in self.flat_quad_settings.iter() {
            let local_index = quad_setting.index;

            let tqi = self
                .bundled_frames
                .get(&local_index)
                .expect("Only existing quads should be queried.");

            // for i in 0..prec_mz_ranges.len() {
            //     if !matches_quad_settings(quad_setting, prec_mz_ranges[i], scan_ranges[i]) {
            //         continue;
            //     }

            //     for (fh, tof_range) in fragment_queries[i].mz_index_ranges.iter() {
            //         let mut local_lambda = |peak| aggregator[i].add(&(RawPeak::from(peak), *fh));
            //         tqi.query_peaks(
            //             *tof_range,
            //             scan_ranges[i],
            //             frame_index_ranges[i],
            //             &mut local_lambda,
            //         )
            //     }
            // }

            aggregator.par_iter_mut().enumerate().for_each(|(i, agg)| {
                if !matches_quad_settings(quad_setting, prec_mz_ranges[i], scan_ranges[i]) {
                    return;
                }

                for (fh, tof_range) in fragment_queries[i].mz_index_ranges.iter() {
                    let mut local_lambda = |peak| agg.add((RawPeak::from(peak), *fh));
                    tqi.query_peaks(
                        *tof_range,
                        scan_ranges[i],
                        frame_index_ranges[i],
                        &mut local_lambda,
                    );
                }
            });
        }
    }
}

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
