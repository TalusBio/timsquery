use super::peak_bucket::PeakBucket;
use super::quad_index::{
    FrameRTTolerance, PeakInQuad, TransposedQuadIndex, TransposedQuadIndexBuilder,
};

use crate::models::frames::raw_peak::RawPeak;
use crate::models::frames::single_quad_settings::expand_quad_settings;
use crate::models::frames::single_quad_settings::SingleQuadrupoleSetting;
use crate::models::queries::FragmentGroupIndexQuery;
use crate::models::queries::PrecursorIndexQuery;
use crate::traits::indexed_data::IndexedData;
use crate::utils::compress_explode::explode_vec;
use crate::utils::display::{glimpse_vec, GlimpseConfig};
use crate::ToleranceAdapter;

use log::{debug, info, trace};
use rayon::prelude::*;
use std::collections::hash_map::Entry;
use std::collections::HashMap;
use std::fmt::Display;
use std::sync::Arc;
use std::time::Instant;
use timsrust::converters::{
    ConvertableDomain, Frame2RtConverter, Scan2ImConverter, Tof2MzConverter,
};
use timsrust::readers::{FrameReader, MetadataReader};
use timsrust::Frame;
use timsrust::Metadata;
use timsrust::TimsRustError;

// TODO break this module apart ... its getting too big for my taste
// - JSP: 2024-11-19

pub struct QuadSplittedTransposedIndex {
    indices: HashMap<Arc<SingleQuadrupoleSetting>, TransposedQuadIndex>,
    flat_quad_settings: Vec<SingleQuadrupoleSetting>,
    rt_converter: Frame2RtConverter,
    mz_converter: Tof2MzConverter,
    im_converter: Scan2ImConverter,
    metadata: Metadata,
}

impl QuadSplittedTransposedIndex {
    pub fn query_peaks(
        &self,
        tof_range: (u32, u32),
        precursor_mz_range: (f64, f64),
        scan_range: Option<(usize, usize)>,
        rt_range: Option<FrameRTTolerance>,
    ) -> impl Iterator<Item = PeakInQuad> + '_ {
        trace!(
            "QuadSplittedTransposedIndex::query_peaks tof_range: {:?}, scan_range: {:?}, rt_range: {:?}, precursor_mz_range: {:?}",
            tof_range,
            scan_range,
            rt_range,
            precursor_mz_range,
        );
        let matching_quads: Vec<SingleQuadrupoleSetting> = self
            .get_matching_quad_settings(precursor_mz_range, scan_range)
            .collect();
        trace!("matching_quads: {:?}", matching_quads);

        let rt_range = rt_range.map(|x| x.to_frame_index_range(&self.rt_converter));

        matching_quads.into_iter().flat_map(move |qs| {
            let tqi = self.indices.get(&qs).unwrap();
            tqi.query_peaks(tof_range, scan_range, rt_range)
        })
    }

    fn get_matching_quad_settings(
        &self,
        precursor_mz_range: (f64, f64),
        scan_range: Option<(usize, usize)>,
    ) -> impl Iterator<Item = SingleQuadrupoleSetting> + '_ {
        self.flat_quad_settings
            .iter()
            .filter(move |qs| {
                (qs.ranges.isolation_low <= precursor_mz_range.1)
                    && (qs.ranges.isolation_high >= precursor_mz_range.0)
            })
            .filter(move |qs| {
                if let Some(scan_range) = scan_range {
                    // This is done for sanity tbh ... sometimes they get flipped
                    // bc the lowest scan is actually the highest 1/k0.
                    let min_scan = qs.ranges.scan_start.min(qs.ranges.scan_end);
                    let max_scan = qs.ranges.scan_start.max(qs.ranges.scan_end);
                    (min_scan <= scan_range.1) && (scan_range.0 <= max_scan)
                } else {
                    true
                }
            })
            .cloned()
    }

    fn queries_from_elution_elements_impl(
        &self,
        tol: &dyn crate::traits::tolerance::Tolerance,
        elution_group: &crate::models::elution_group::ElutionGroup,
    ) -> FragmentGroupIndexQuery {
        let rt_range = tol.rt_range(elution_group.rt_seconds);
        let mobility_range = tol.mobility_range(elution_group.mobility);
        let precursor_mz_range = tol.mz_range(elution_group.precursor_mz);
        let quad_range = tol.quad_range(elution_group.precursor_mz, elution_group.precursor_charge);

        let mz_index_range = (
            self.mz_converter.invert(precursor_mz_range.0) as u32,
            self.mz_converter.invert(precursor_mz_range.1) as u32,
        );
        let mobility_range = match mobility_range {
            Some(mobility_range) => mobility_range,
            None => (self.metadata.lower_im as f32, self.metadata.upper_im as f32),
        };
        let mobility_index_range = (
            self.im_converter.invert(mobility_range.0) as usize,
            self.im_converter.invert(mobility_range.1) as usize,
        );
        let rt_range = match rt_range {
            Some(rt_range) => rt_range,
            None => (self.metadata.lower_rt as f32, self.metadata.upper_rt as f32),
        };
        let frame_index_range = (
            self.rt_converter.invert(rt_range.0) as usize,
            self.rt_converter.invert(rt_range.1) as usize,
        );

        assert!(frame_index_range.0 <= frame_index_range.1);
        assert!(mz_index_range.0 <= mz_index_range.1);
        // Since mobilities get mixed up bc low scan ranges are high 1/k0, I
        // Just make sure they are sorted here.
        let mobility_index_range = (
            mobility_index_range.0.min(mobility_index_range.1),
            mobility_index_range.1.max(mobility_index_range.0),
        );

        let precursor_query = PrecursorIndexQuery {
            frame_index_range,
            mz_index_range,
            mobility_index_range,
            isolation_mz_range: quad_range,
        };

        // TODO: change this unwrap and use explicitly the lack of fragment mzs.
        // Does that mean its onlyt a precursor query?
        // Why is it an option?
        let fragment_mzs = elution_group.fragment_mzs.as_ref().unwrap();
        let mut fqs = Vec::with_capacity(fragment_mzs.len());
        for mz_range in fragment_mzs {
            let mz_range = tol.mz_range(*mz_range);
            fqs.push((
                self.mz_converter.invert(mz_range.0) as u32,
                self.mz_converter.invert(mz_range.1) as u32,
            ));
        }

        FragmentGroupIndexQuery {
            mz_index_ranges: fqs,
            precursor_query,
        }
    }
}

impl Display for QuadSplittedTransposedIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut disp_str = String::new();
        disp_str.push_str("QuadSplittedTransposedIndex\n");

        disp_str.push_str(&format!("rt_converter: ... not showing ...\n",));
        disp_str.push_str(&format!("mz_converter: {:?}\n", self.mz_converter));
        disp_str.push_str(&format!("im_converter: {:?}\n", self.im_converter));
        disp_str.push_str("flat_quad_settings: \n");
        disp_str.push_str(&glimpse_vec(
            &self.flat_quad_settings,
            Some(GlimpseConfig {
                max_items: 10,
                padding: 4,
                new_line: true,
            }),
        ));
        let mut num_shown = 0;
        for (qs, tqi) in self.indices.iter() {
            disp_str.push_str(&format!(" - {}: \n", qs));
            disp_str.push_str(&format!(" -- {}\n", tqi));
            num_shown += 1;
            if num_shown > 5 {
                disp_str.push_str(&format!(" ........ len = {}\n", self.indices.len()));
                break;
            }
        }

        write!(f, "{}", disp_str)
    }
}

impl QuadSplittedTransposedIndex {
    pub fn from_path(path: &str) -> Result<Self, TimsRustError> {
        let st = Instant::now();
        info!("Building transposed quad index from path {}", path);
        let tmp = QuadSplittedTransposedIndexBuilder::from_path(path)?;
        let out = tmp.build();
        let elapsed = st.elapsed();
        info!("Transposed quad index built in {:#?}", elapsed);
        debug!("{}", out);
        Ok(out)
    }
}

pub struct QuadSplittedTransposedIndexBuilder {
    indices: HashMap<SingleQuadrupoleSetting, TransposedQuadIndexBuilder>,
    rt_converter: Option<Frame2RtConverter>,
    mz_converter: Option<Tof2MzConverter>,
    im_converter: Option<Scan2ImConverter>,
    metadata: Option<Metadata>,
    // TODO use during build to make sure we
    // have a the right number of peaks in the end.
    // ... this means I need to implement len for TransposedQuadIndexBuilder
    added_peaks: u64,
}

impl QuadSplittedTransposedIndexBuilder {
    fn new() -> Self {
        Self {
            indices: HashMap::new(),
            rt_converter: None,
            mz_converter: None,
            im_converter: None,
            metadata: None,
            added_peaks: 0,
        }
    }

    fn add_frame(&mut self, frame: Frame) {
        let expanded_quad_settings = expand_quad_settings(&frame.quadrupole_settings);
        let exploded_scans = explode_vec(&frame.scan_offsets);

        for qs in expanded_quad_settings {
            // Add key if it doesnt exist ...

            if let Entry::Vacant(e) = self.indices.entry(qs) {
                let max_tof = frame.tof_indices.iter().max().unwrap();
                trace!(
                    "Adding new transposed quad index for qs {:?} with max tof {}",
                    qs,
                    max_tof
                );
                let new_index = TransposedQuadIndexBuilder::new(qs);
                e.insert(new_index);
            }

            let peak_start = frame.scan_offsets[qs.ranges.scan_start];
            let peak_end = frame.scan_offsets[qs.ranges.scan_end + 1];

            let int_slice = frame.intensities[peak_start..peak_end].to_vec();
            let tof_slice = frame.tof_indices[peak_start..peak_end].to_vec();
            let expanded_scan_slice = exploded_scans[peak_start..peak_end].to_vec();

            let frame_index = frame.index;
            let frame_rt = frame.rt;
            self.added_peaks += int_slice.len() as u64;

            self.indices.get_mut(&qs).unwrap().add_frame_slice(
                int_slice,
                tof_slice,
                expanded_scan_slice,
                frame_index,
                frame_rt,
            );
        }
    }

    fn from_path(path: &str) -> Result<Self, TimsRustError> {
        let file_reader = FrameReader::new(path)?;

        let sql_path = std::path::Path::new(path).join("analysis.tdf");
        let meta_converters = MetadataReader::new(&sql_path)?;

        let out_meta_converters = meta_converters.clone();
        let mut final_out = Self {
            indices: HashMap::new(),
            rt_converter: Some(meta_converters.rt_converter),
            mz_converter: Some(meta_converters.mz_converter),
            im_converter: Some(meta_converters.im_converter),
            metadata: Some(out_meta_converters),
            added_peaks: 0,
        };

        let out2: Result<Vec<Self>, TimsRustError> = file_reader
            .get_all()
            .into_par_iter()
            .chunks(100)
            .map(|frames| {
                let mut out = Self::new();
                for frame in frames {
                    let frame = frame?;
                    out.add_frame(frame);
                }
                Ok(out)
            })
            .collect();

        let out2 = out2?.into_iter().fold(Self::new(), |mut x, y| {
            x.fold(y);
            x
        });
        final_out.fold(out2);

        Ok(final_out)
    }

    fn fold(&mut self, other: Self) {
        for (qs, builder) in other.indices.into_iter() {
            self.indices
                .entry(qs)
                .and_modify(|bl| bl.fold(builder.clone()))
                .or_insert(builder);
        }
        self.added_peaks += other.added_peaks;
    }

    fn build(self) -> QuadSplittedTransposedIndex {
        let mut indices = HashMap::new();
        let mut flat_quad_settings = Vec::new();
        let built: Vec<(TransposedQuadIndex, SingleQuadrupoleSetting)> = self
            .indices
            .into_par_iter()
            .map(|(qs, builder)| (builder.build(), qs))
            .collect();

        for (qi, qs) in built.into_iter() {
            let qa: Arc<SingleQuadrupoleSetting> = Arc::new(qs);
            indices.insert(qa.clone(), qi);
            flat_quad_settings.push(qs);
        }

        flat_quad_settings.sort_by(|a, b| {
            a.ranges
                .isolation_mz
                .partial_cmp(&b.ranges.isolation_mz)
                .unwrap()
        });

        QuadSplittedTransposedIndex {
            indices,
            flat_quad_settings,
            rt_converter: self.rt_converter.unwrap(),
            mz_converter: self.mz_converter.unwrap(),
            im_converter: self.im_converter.unwrap(),
            metadata: self.metadata.unwrap(),
        }
    }
}

// TODO rewrite for the new data structure ... (btree instead of slice of optional)
fn display_opt_peak_bucket(opt_peak_bucket: &Option<PeakBucket>) -> String {
    match opt_peak_bucket {
        Some(peak_bucket) => format!("{}", peak_bucket),
        None => "None".to_string(),
    }
}

fn display_opt_peak_bucket_vec(opt_peak_buckets: &[Option<PeakBucket>]) -> String {
    let mut out = String::new();
    let num_none = opt_peak_buckets.iter().filter(|x| x.is_none()).count();
    let num_some = opt_peak_buckets.iter().filter(|x| x.is_some()).count();

    let mut max_peaks = 0;
    let mut tof_with_max = 0;
    for (i, opt_peak_bucket) in opt_peak_buckets.iter().enumerate() {
        if let Some(peak_bucket) = opt_peak_bucket {
            if max_peaks < peak_bucket.len() {
                max_peaks = peak_bucket.len();
                tof_with_max = i;
            }
        }
    }

    out.push_str(&format!(
        "PeakBuckets: num_none: {}, num_some: {}, max_peaks: {}, tof_with_max: {}\n",
        num_none, num_some, max_peaks, tof_with_max
    ));
    for (i, opt_peak_bucket) in opt_peak_buckets.iter().enumerate() {
        out.push_str(&format!(
            " - {}: {}\n",
            i,
            display_opt_peak_bucket(opt_peak_bucket)
        ));
        if i > 3 {
            out.push_str(&format!(" - ... len = {}\n", opt_peak_buckets.len()));
            break;
        }
    }

    if tof_with_max > 0 {
        out.push_str(&format!(
            " - Bucket with max tof: {} {}\n",
            tof_with_max,
            display_opt_peak_bucket(&opt_peak_buckets[tof_with_max])
        ));
    }

    out
}

impl IndexedData<FragmentGroupIndexQuery, RawPeak> for QuadSplittedTransposedIndex {
    fn query(&self, fragment_query: &FragmentGroupIndexQuery) -> Vec<RawPeak> {
        let precursor_mz_range = (
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
        );
        let scan_range = Some(fragment_query.precursor_query.mobility_index_range);
        let frame_index_range = Some(FrameRTTolerance::FrameIndex(
            fragment_query.precursor_query.frame_index_range,
        ));

        fragment_query
            .mz_index_ranges
            .iter()
            .flat_map(|tof_range| {
                self.query_peaks(
                    *tof_range,
                    precursor_mz_range,
                    scan_range,
                    frame_index_range,
                )
                .map(RawPeak::from)
            })
            .collect()
    }

    fn add_query<O, AG: crate::Aggregator<RawPeak, Output = O>>(
        &self,
        fragment_query: &FragmentGroupIndexQuery,
        aggregator: &mut AG,
    ) {
        let precursor_mz_range = (
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
        );
        let scan_range = Some(fragment_query.precursor_query.mobility_index_range);
        let frame_index_range = Some(FrameRTTolerance::FrameIndex(
            fragment_query.precursor_query.frame_index_range,
        ));

        fragment_query.mz_index_ranges.iter().for_each(|tof_range| {
            self.query_peaks(
                *tof_range,
                precursor_mz_range,
                scan_range,
                frame_index_range,
            )
            .for_each(|peak| aggregator.add(&RawPeak::from(peak)));
        })
    }

    fn add_query_multi_group<O, AG: crate::Aggregator<RawPeak, Output = O>>(
        &self,
        fragment_queries: &[FragmentGroupIndexQuery],
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
                let frame_index_range = Some(FrameRTTolerance::FrameIndex(
                    fragment_query.precursor_query.frame_index_range,
                ));

                for tof_range in fragment_query.mz_index_ranges.clone().into_iter() {
                    self.query_peaks(tof_range, precursor_mz_range, scan_range, frame_index_range)
                        .for_each(|peak| agg.add(&RawPeak::from(peak)));
                }
            });
    }
}

// Copy pasting for now ... TODO refactor or delete
// ============================================================================

impl IndexedData<FragmentGroupIndexQuery, (RawPeak, usize)> for QuadSplittedTransposedIndex {
    fn query(&self, fragment_query: &FragmentGroupIndexQuery) -> Vec<(RawPeak, usize)> {
        let precursor_mz_range = (
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
        );
        let scan_range = Some(fragment_query.precursor_query.mobility_index_range);
        let frame_index_range = Some(FrameRTTolerance::FrameIndex(
            fragment_query.precursor_query.frame_index_range,
        ));

        fragment_query
            .mz_index_ranges
            .iter()
            .enumerate()
            .flat_map(|(ind, tof_range)| {
                let out: Vec<(RawPeak, usize)> = self
                    .query_peaks(
                        *tof_range,
                        precursor_mz_range,
                        scan_range,
                        frame_index_range,
                    )
                    .map(|x| (RawPeak::from(x), ind))
                    .collect();

                out
            })
            .collect()
    }

    fn add_query<O, AG: crate::Aggregator<(RawPeak, usize), Output = O>>(
        &self,
        fragment_query: &FragmentGroupIndexQuery,
        aggregator: &mut AG,
    ) {
        let precursor_mz_range = (
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
        );
        let scan_range = Some(fragment_query.precursor_query.mobility_index_range);
        let frame_index_range = Some(FrameRTTolerance::FrameIndex(
            fragment_query.precursor_query.frame_index_range,
        ));

        fragment_query
            .mz_index_ranges
            .iter()
            .enumerate()
            .for_each(|(ind, tof_range)| {
                self.query_peaks(
                    *tof_range,
                    precursor_mz_range,
                    scan_range,
                    frame_index_range,
                )
                .for_each(|peak| aggregator.add(&(RawPeak::from(peak), ind)));
            })
    }

    fn add_query_multi_group<O, AG: crate::Aggregator<(RawPeak, usize), Output = O>>(
        &self,
        fragment_queries: &[FragmentGroupIndexQuery],
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
                let frame_index_range = Some(FrameRTTolerance::FrameIndex(
                    fragment_query.precursor_query.frame_index_range,
                ));

                for (ind, tof_range) in fragment_query
                    .mz_index_ranges
                    .clone()
                    .into_iter()
                    .enumerate()
                {
                    self.query_peaks(tof_range, precursor_mz_range, scan_range, frame_index_range)
                        .for_each(|peak| agg.add(&(RawPeak::from(peak), ind)));
                }
            });
    }
}

// ============================================================================

impl ToleranceAdapter<FragmentGroupIndexQuery> for QuadSplittedTransposedIndex {
    fn query_from_elution_group(
        &self,
        tol: &dyn crate::traits::tolerance::Tolerance,
        elution_group: &crate::models::elution_group::ElutionGroup,
    ) -> FragmentGroupIndexQuery {
        self.queries_from_elution_elements_impl(tol, elution_group)
    }
}
