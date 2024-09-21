use crate::models::frames::raw_peak::RawPeak;
use crate::models::frames::single_quad_settings::expand_quad_settings;
use crate::models::frames::single_quad_settings::SingleQuadrupoleSetting;
use crate::models::queries::FragmentGroupIndexQuery;
use crate::models::queries::PrecursorIndexQuery;
use crate::sort_by_indices_multi;
use crate::traits::indexed_data::IndexedData;
use crate::utils::compress_explode::{compress_vec, explode_vec};
use crate::utils::display::{glimpse_vec, GlimpseConfig};
use crate::utils::sorting::argsort_by;
use crate::ToleranceAdapter;
use core::panic;
use log::trace;
use log::{debug, info};
use std::collections::BTreeMap;
use std::time::Instant;
use timsrust::Metadata;

use rayon::prelude::*;

use std::collections::HashMap;
use std::fmt::Display;
use std::sync::Arc;
use timsrust::converters::{
    ConvertableDomain, Frame2RtConverter, Scan2ImConverter, Tof2MzConverter,
};
use timsrust::readers::{FrameReader, MetadataReader};
use timsrust::Frame;
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

#[derive(Debug, Clone, Copy)]
pub enum FrameRTTolerance {
    FrameIndex((usize, usize)),
    Seconds((f64, f64)),
}

impl FrameRTTolerance {
    pub fn to_frame_index_range(self, rt_converter: &Frame2RtConverter) -> FrameRTTolerance {
        match self {
            FrameRTTolerance::FrameIndex(_) => self,
            FrameRTTolerance::Seconds((low, high)) => {
                let low = rt_converter.invert(low).round() as usize;
                let high = rt_converter.invert(high).round() as usize;
                FrameRTTolerance::FrameIndex((low, high))
            }
        }
    }
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
        // let matching_quads = self.get_matching_quad_settings(precursor_mz_range, scan_range);
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
        // assert!(mobility_index_range.0 <= mobility_index_range.1);
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
}

impl QuadSplittedTransposedIndexBuilder {
    fn new() -> Self {
        Self {
            indices: HashMap::new(),
            rt_converter: None,
            mz_converter: None,
            im_converter: None,
            metadata: None,
        }
    }

    fn add_frame(&mut self, frame: Frame) {
        let expanded_quad_settings = expand_quad_settings(&frame.quadrupole_settings);
        let exploded_scans = explode_vec(&frame.scan_offsets);

        for qs in expanded_quad_settings {
            // Add key if it doesnt exist ...
            if !self.indices.contains_key(&qs) {
                let max_tof = frame.tof_indices.iter().max().unwrap();
                debug!(
                    "Adding new transposed quad index for qs {:?} with max tof {}",
                    qs, max_tof
                );
                let new_index = TransposedQuadIndexBuilder::new(qs.clone());
                self.indices.insert(qs, new_index);
            }

            let peak_start = frame.scan_offsets[qs.ranges.scan_start];
            let peak_end = frame.scan_offsets[qs.ranges.scan_end + 1];

            let int_slice = frame.intensities[peak_start..peak_end].to_vec();
            let tof_slice = frame.tof_indices[peak_start..peak_end].to_vec();
            let expanded_scan_slice = exploded_scans[peak_start..peak_end].to_vec();

            let frame_index = frame.index;
            let frame_rt = frame.rt;

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
        let file_reader = FrameReader::new(&path)?;

        let sql_path = std::path::Path::new(path).join("analysis.tdf");
        let meta_converters = MetadataReader::new(&sql_path)?;

        let mut final_out = Self {
            indices: HashMap::new(),
            rt_converter: Some(meta_converters.rt_converter.clone()),
            mz_converter: Some(meta_converters.mz_converter.clone()),
            im_converter: Some(meta_converters.im_converter.clone()),
            metadata: Some(meta_converters.clone()),
        };

        // If I implement fold ... I can do this in parallel.
        // By splitting the frames into batches and doing the folds in parallel.
        //
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

        // for frame in file_reader.get_all().into_iter().progress() {
        //     let frame = frame?;
        //     out.add_frame(frame);
        // }

        Ok(final_out)
    }

    fn fold(&mut self, other: Self) {
        for (qs, builder) in other.indices.into_iter() {
            self.indices
                .entry(qs)
                .and_modify(|bl| bl.fold(builder.clone()))
                .or_insert(builder);
        }
    }

    fn build(self) -> QuadSplittedTransposedIndex {
        let mut indices = HashMap::new();
        let mut flat_quad_settings = Vec::new();
        // for (qs, builder) in self.indices.into_iter().progress() {
        //     let tmp = Arc::new(qs);
        //     indices.insert(tmp.clone(), builder.build());
        //     flat_quad_settings.push(qs);
        // }

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

#[derive(Debug)]
struct TransposedQuadIndex {
    quad_settings: SingleQuadrupoleSetting,
    frame_rts: Vec<f64>,
    frame_indices: Vec<usize>,
    peak_buckets: BTreeMap<u32, PeakBucket>,
    // 72 bytes for peak Option<bucket> (24/vec)
    // Conservatvely ... 30_000_000 elements can be layed out in a vec.
    // 638425 is the max observed val for a tof index ... So we dont need an
    // additional bucketing.
}

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
            if max_peaks < peak_bucket.intensities.len() {
                max_peaks = peak_bucket.intensities.len();
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

impl Display for TransposedQuadIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "TransposedQuadIndex\n quad_settings: {}\n frame_indices: {}\n frame_rts: {}\n peak_buckets: NOT SHOWING\n",
            self.quad_settings,
            glimpse_vec(&self.frame_indices, Some(GlimpseConfig { max_items: 10, padding: 2, new_line: true })),
            glimpse_vec(&self.frame_rts, Some(GlimpseConfig { max_items: 10, padding: 2, new_line: true })),
            // display_opt_peak_bucket_vec(&self.peak_buckets),
        )
    }
}

#[derive(Debug, Clone, Copy)]
pub struct PeakInQuad {
    pub scan_index: usize,
    pub intensity: u32,
    pub retention_time: f32,
    pub tof_index: u32,
}

impl PeakInQuad {
    pub fn from_peak_in_bucket(peak_in_bucket: PeakInBucket, tof_index: u32) -> Self {
        Self {
            scan_index: peak_in_bucket.scan_index,
            intensity: peak_in_bucket.intensity,
            retention_time: peak_in_bucket.retention_time,
            tof_index,
        }
    }
}

impl TransposedQuadIndex {
    pub fn query_peaks(
        &self,
        tof_range: (u32, u32),
        scan_range: Option<(usize, usize)>,
        rt_range: Option<FrameRTTolerance>,
        // ) -> Vec<PeakInQuad> {
    ) -> impl Iterator<Item = PeakInQuad> + '_ {
        trace!(
            "TransposedQuadIndex::query_peaks tof_range: {:?}, scan_range: {:?}, rt_range: {:?}",
            tof_range,
            scan_range,
            rt_range
        );
        // TODO reimplement as an iterator ...
        // This version is not compatible with the borrow checker unless I collect the vec...
        // which will do for now for prototyping.
        let frame_index_range = self.convert_to_local_frame_range(rt_range);

        self.peak_buckets
            .range(tof_range.0..tof_range.1)
            .flat_map(move |(tof_index, pb)| {
                pb.query_peaks(scan_range, frame_index_range)
                    .map(move |p| PeakInQuad::from_peak_in_bucket(p, *tof_index))
            })
        // Coult I just return an Arc<[intensities]> + ...
        // If I made the peak buckets sparse, I could make it ... not be an option.
    }

    fn convert_to_local_frame_range(
        &self,
        rt_range: Option<FrameRTTolerance>,
    ) -> Option<(f32, f32)> {
        let frame_index_range = match rt_range {
            Some(FrameRTTolerance::Seconds((rt_low, rt_high))) => {
                Some((rt_low as f32, rt_high as f32))
            }
            Some(FrameRTTolerance::FrameIndex((frame_low, frame_high))) => {
                let frame_id_start = self
                    .frame_indices
                    .binary_search_by(|x| x.cmp(&frame_low))
                    .unwrap_or_else(|x| x);

                let frame_id_end = self
                    .frame_indices
                    .binary_search_by(|x| x.cmp(&frame_high))
                    .unwrap_or_else(|x| x);

                // TODO consider throwing a warning if we are
                // out of bounds here.
                Some((
                    self.frame_rts[frame_id_start.min(self.frame_rts.len() - 1)] as f32,
                    self.frame_rts[frame_id_end.min(self.frame_rts.len() - 1)] as f32,
                ))
            }
            None => None,
        };

        if cfg!(debug_assertions) {
            if let Some((low, high)) = frame_index_range {
                debug_assert!(low <= high);
            }
        }

        frame_index_range
    }
}

#[derive(Debug, Clone)]
struct TransposedQuadIndexBuilder {
    quad_settings: SingleQuadrupoleSetting,
    int_slices: Vec<Vec<u32>>,
    tof_slices: Vec<Vec<u32>>,
    scan_slices: Vec<Vec<usize>>,
    frame_indices: Vec<usize>,
    frame_rts: Vec<f64>,
}

impl TransposedQuadIndexBuilder {
    fn new(quad_settings: SingleQuadrupoleSetting) -> Self {
        Self {
            quad_settings,
            int_slices: Vec::new(),
            tof_slices: Vec::new(),
            scan_slices: Vec::new(),
            frame_indices: Vec::new(),
            frame_rts: Vec::new(),
        }
    }

    fn fold(&mut self, other: Self) {
        self.int_slices.extend(other.int_slices);
        self.tof_slices.extend(other.tof_slices);
        self.scan_slices.extend(other.scan_slices);
        self.frame_indices.extend(other.frame_indices);
        self.frame_rts.extend(other.frame_rts);
    }

    fn add_frame_slice(
        &mut self,
        int_slice: Vec<u32>,
        tof_slice: Vec<u32>,
        expanded_scan_slice: Vec<usize>,
        frame_index: usize,
        frame_rt: f64,
    ) {
        assert!(int_slice.len() == tof_slice.len());
        assert!(int_slice.len() == expanded_scan_slice.len());

        self.int_slices.push(int_slice);
        self.tof_slices.push(tof_slice);
        self.scan_slices.push(expanded_scan_slice);
        self.frame_indices.push(frame_index);
        self.frame_rts.push(frame_rt);
    }

    fn build(self) -> TransposedQuadIndex {
        let st = Instant::now();
        let max_tof = *self
            .tof_slices
            .iter()
            .map(|x| x.iter().max().unwrap())
            .max()
            .unwrap();
        let mut tof_counts = vec![0; max_tof as usize + 1];
        let mut tot_peaks = 0;
        for tof_slice in self.tof_slices.iter() {
            for tof in tof_slice.iter() {
                tof_counts[*tof as usize] += 1;
                tot_peaks += 1;
            }
        }

        let mut peak_buckets = HashMap::with_capacity((max_tof + 1).try_into().unwrap());
        for (tof, count) in tof_counts.clone().into_iter().enumerate() {
            if count > 0 {
                let peak_bucket = PeakBucketBuilder::new(count, tof as u32);
                peak_buckets.insert(tof as u32, peak_bucket);
            } else {
                continue;
            }
        }

        let out_rts = self.frame_rts.clone();
        let out_indices = self.frame_indices.clone();
        let mut added_peaks = 0;

        // TODO rewrite this as indices instead zipping the slices.
        // I could even try to concat + sort + pass ... I would need to implement
        // 'add peak slice'

        // let mut expanded_frame_rts =
        //     Vec::with_capacity(self.tof_slices.iter().map(|x| x.len()).sum());

        let aps = Instant::now();

        for slice_ind in 0..self.frame_indices.len() {
            let int_slice = &self.int_slices[slice_ind];
            let tof_slice = &self.tof_slices[slice_ind];
            let scan_slice = &self.scan_slices[slice_ind];
            // let frame_index = self.frame_indices[slice_ind];
            let frame_rt = self.frame_rts[slice_ind] as f32;

            for ((inten, tof), scan) in int_slice
                .into_iter()
                .zip(tof_slice.into_iter())
                .zip(scan_slice.into_iter())
            {
                peak_buckets
                    .get_mut(&tof)
                    .expect("Should have just added it...")
                    .add_peak(*scan, *inten, frame_rt);

                added_peaks += 1;
            }
        }
        let aps = aps.elapsed();

        if added_peaks != tot_peaks {
            println!("TransposedQuadIndex::add_frame_slice failed at peak count check, expected: {}, real: {}", tot_peaks, added_peaks);
            panic!("TransposedQuadIndex::add_frame_slice failed at peak count check");
        }

        for (tof, count) in tof_counts.into_iter().enumerate() {
            if count > 0 {
                let curr_bucket = peak_buckets.get(&(tof as u32)).unwrap();
                let real_count = curr_bucket.intensities.len();
                if real_count != count {
                    println!("TransposedQuadIndex::build failed at tof bucket count check, expected: {}, real: {}", count, real_count);
                    println!("Bucket -> {:?}", curr_bucket);

                    panic!("TransposedQuadIndex::build failed at tof bucket count check");
                };
            } else {
                continue;
            }
        }

        let bb_st = Instant::now();
        let peak_bucket = BTreeMap::from_par_iter(
            peak_buckets
                .into_par_iter()
                .map(|(tof, pb)| (tof, pb.build())),
        );

        let bbe = bb_st.elapsed();
        let elapsed = st.elapsed();
        info!(
            "TransposedQuadIndex::add_frame_slice adding peaks took {:#?} for {} peaks",
            aps, tot_peaks
        );
        info!(
            "TransposedQuadIndex::add_frame_slice building buckets took {:#?}",
            bbe
        );
        info!("TransposedQuadIndex::build took {:#?}", elapsed);
        TransposedQuadIndex {
            quad_settings: self.quad_settings,
            frame_rts: out_rts,
            frame_indices: out_indices,
            peak_buckets: peak_bucket,
        }
    }
}

pub struct PeakInBucket {
    pub scan_index: usize,
    pub intensity: u32,
    pub retention_time: f32,
}

#[derive(Debug, Clone, Copy)]
pub enum PeakBucketMode {
    Compressed,
    Sorted,
}

#[derive(Debug)]
struct PeakBucket {
    intensities: Vec<u32>,
    retention_times: Vec<f32>,
    scan_offsets: Vec<usize>,
    tof_index: u32,
    mode: PeakBucketMode,
}

impl Display for PeakBucket {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "PeakBucket tof={}: \n    len={},\n    retention_times={},\n    scan_offsets={},\n    intensities={}, \n    mode={:?}",
            self.tof_index,
            self.len(),
            glimpse_vec(
                &self.retention_times,
                Some(GlimpseConfig {
                    max_items: 10,
                    padding: 2,
                    new_line: false
                })
            ),
            glimpse_vec(
                &self.scan_offsets,
                Some(GlimpseConfig {
                    max_items: 10,
                    padding: 2,
                    new_line: false
                })
            ),
            glimpse_vec(
                &self.intensities,
                Some(GlimpseConfig {
                    max_items: 10,
                    padding: 2,
                    new_line: false
                })
            ),
            self.mode,
        )
    }
}

#[derive(Debug)]
struct PeakBucketBuilder {
    intensities: Vec<u32>,
    retention_times: Vec<f32>,
    scan_offsets: Vec<usize>,
    tof_index: u32,
}

impl PeakBucketBuilder {
    fn new(capacity: usize, tof_index: u32) -> Self {
        Self {
            intensities: Vec::with_capacity(capacity),
            retention_times: Vec::with_capacity(capacity),
            scan_offsets: Vec::with_capacity(capacity),
            tof_index,
        }
    }

    fn add_slices(&mut self, int_slice: &[u32], rt_slice: &[f32], scan_slice: &[usize]) {
        self.intensities.extend(int_slice);
        self.retention_times.extend(rt_slice);
        self.scan_offsets.extend(scan_slice);
    }

    fn add_peak(&mut self, scan_index: usize, intensity: u32, retention_time: f32) {
        self.intensities.push(intensity);
        self.retention_times.push(retention_time);
        self.scan_offsets.push(scan_index);
    }

    fn fold(&mut self, other: Self) {
        assert_eq!(self.tof_index, other.tof_index);
        self.intensities.extend(other.intensities);
        self.retention_times.extend(other.retention_times);
        self.scan_offsets.extend(other.scan_offsets);
    }

    fn build(mut self) -> PeakBucket {
        let mut indices = argsort_by(&self.scan_offsets, |x| *x);
        sort_by_indices_multi!(
            &mut indices,
            &mut self.scan_offsets,
            &mut self.retention_times,
            &mut self.intensities
        );
        // TODO consider if I really need to compress this.
        let out = if self.scan_offsets.len() > 1000 {
            let compressed = compress_vec(&self.scan_offsets);
            PeakBucket {
                intensities: self.intensities,
                retention_times: self.retention_times,
                scan_offsets: compressed,
                mode: PeakBucketMode::Compressed,
                tof_index: self.tof_index,
            }
        } else {
            PeakBucket {
                intensities: self.intensities,
                retention_times: self.retention_times,
                scan_offsets: self.scan_offsets,
                mode: PeakBucketMode::Sorted,
                tof_index: self.tof_index,
            }
        };
        out
    }
}

impl PeakBucket {
    pub fn len(&self) -> usize {
        self.intensities.len()
    }

    pub fn query_peaks(
        &self,
        scan_range: Option<(usize, usize)>,
        rt_range: Option<(f32, f32)>,
    ) -> impl Iterator<Item = PeakInBucket> + '_ {
        match self.mode {
            PeakBucketMode::Compressed => self
                .query_peaks_compressed(scan_range, rt_range)
                .collect::<Vec<PeakInBucket>>()
                .into_iter(),
            PeakBucketMode::Sorted => self
                .query_peaks_sorted(scan_range, rt_range)
                .collect::<Vec<PeakInBucket>>()
                .into_iter(),
        }
    }

    pub fn query_peaks_sorted(
        &self,
        scan_range: Option<(usize, usize)>,
        rt_range: Option<(f32, f32)>,
    ) -> impl Iterator<Item = PeakInBucket> + '_ {
        let (scan_min, scan_max) = match scan_range {
            Some((scan_low, scan_high)) => (scan_low, scan_high),
            None => (
                0,
                *(self
                    .scan_offsets
                    .last()
                    .expect("It should not be possible to build an empty PeakBucket")),
            ),
        };

        // TODO do a binary search if the data is large-ish.
        let mut start_min = 0;
        let mut end_max = self.len();
        while start_min < end_max && self.scan_offsets[start_min] < scan_min {
            start_min += 1;
        }
        while start_min < end_max && self.scan_offsets[end_max - 1] >= scan_max {
            end_max -= 1;
        }

        (start_min..end_max)
            .map(move |x| {
                let scan_index = self.scan_offsets[x];
                let retention_time = self.retention_times[x];
                match rt_range {
                    Some((low, high)) => {
                        if retention_time < low || retention_time > high {
                            return None;
                        }
                    }
                    None => {}
                }
                Some(PeakInBucket {
                    scan_index,
                    intensity: self.intensities[x],
                    retention_time,
                })
            })
            .filter_map(|x| x)
    }

    pub fn query_peaks_compressed(
        &self,
        scan_range: Option<(usize, usize)>,
        rt_range: Option<(f32, f32)>,
    ) -> impl Iterator<Item = PeakInBucket> + '_ {
        let scan_range = match scan_range {
            Some((scan_low, scan_high)) => {
                (scan_low..scan_high.min(self.scan_offsets.len() - 1)).into_iter()
            }
            None => 0..self.scan_offsets.len(),
        };
        scan_range
            .flat_map(move |scan_index| {
                let peak_index_start = self.scan_offsets[scan_index];
                let peak_index_end = self.scan_offsets[scan_index + 1];

                (peak_index_start..peak_index_end)
                    .into_iter()
                    .map(move |peak_index| {
                        let retention_time = self.retention_times[peak_index];

                        match rt_range {
                            Some((low, high)) => {
                                if retention_time < low || retention_time > high {
                                    return None;
                                }
                            }
                            None => {}
                        }

                        let intensity = self.intensities[peak_index];
                        Some(PeakInBucket {
                            scan_index,
                            intensity,
                            retention_time,
                        })
                    })
            })
            .filter_map(|x| x)
    }

    fn verify(&self) -> bool {
        if self.intensities.len() != self.retention_times.len() {
            println!("PeakBucket::verify failed at length check");
            return false;
        }
        match self.mode {
            PeakBucketMode::Compressed => {
                match self.scan_offsets.last() {
                    Some(last) => {
                        if *last != self.retention_times.len() {
                            println!(
                                "PeakBucket::verify failed at last scan check, in compressed mode."
                            );
                            return false;
                        }
                    }
                    None => {
                        println!("PeakBucket::verify failed at last scan check, no scans");
                        return false;
                    }
                }

                for i in 1..self.scan_offsets.len() {
                    if self.scan_offsets[i] < self.scan_offsets[i - 1] {
                        println!("PeakBucket::verify failed at scan order check");
                        return false;
                    }
                }
            }
            PeakBucketMode::Sorted => {
                if self.intensities.len() != self.scan_offsets.len() {
                    println!("PeakBucket::verify failed at length check on sorted mode");
                    return false;
                }
                for i in 1..self.scan_offsets.len() {
                    if self.scan_offsets[i] < self.scan_offsets[i - 1] {
                        println!("PeakBucket::verify failed at scan order check");
                        return false;
                    }
                }
            }
        };

        true
    }
}

impl From<PeakInQuad> for RawPeak {
    fn from(peak_in_quad: PeakInQuad) -> Self {
        RawPeak {
            scan_index: peak_in_quad.scan_index,
            tof_index: peak_in_quad.tof_index,
            intensity: peak_in_quad.intensity,
            retention_time: peak_in_quad.retention_time,
        }
    }
}

impl IndexedData<FragmentGroupIndexQuery, RawPeak> for QuadSplittedTransposedIndex {
    fn query(&self, fragment_query: &FragmentGroupIndexQuery) -> Option<Vec<RawPeak>> {
        let precursor_mz_range = (
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
        );
        let scan_range = Some(fragment_query.precursor_query.mobility_index_range);
        let frame_index_range = Some(FrameRTTolerance::FrameIndex(
            fragment_query.precursor_query.frame_index_range,
        ));

        let out = fragment_query
            .mz_index_ranges
            .iter()
            .flat_map(|tof_range| {
                self.query_peaks(
                    *tof_range,
                    precursor_mz_range,
                    scan_range,
                    frame_index_range,
                )
                .map(|peak| RawPeak::from(peak))
                .collect::<Vec<RawPeak>>()
            })
            .collect();

        // let mut out = Vec::new();
        // for tof_range in fragment_query.mz_index_ranges.clone().into_iter() {
        //     self.query_peaks(tof_range, precursor_mz_range, scan_range, frame_index_range)
        //         .for_each(|peak| out.push(RawPeak::from(peak)));
        // }
        Some(out)
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

        // for tof_range in fragment_query.mz_index_ranges.clone().into_iter() {
        //     self.query_peaks(tof_range, precursor_mz_range, scan_range, frame_index_range)
        //         .for_each(|peak| aggregator.add(&RawPeak::from(peak)));
        // }

        fragment_query
            .mz_index_ranges
            .iter()
            .map(|tof_range| {
                let mut out = Vec::new();
                self.query_peaks(
                    *tof_range,
                    precursor_mz_range,
                    scan_range,
                    frame_index_range,
                )
                .for_each(|peak| out.push(RawPeak::from(peak)));
                out
            })
            .flatten()
            .collect::<Vec<RawPeak>>()
            .into_iter()
            .for_each(|x| aggregator.add(&x));
    }

    fn add_query_multi_group<O, AG: crate::Aggregator<RawPeak, Output = O>>(
        &self,
        fragment_queries: &[FragmentGroupIndexQuery],
        aggregator: &mut [AG],
    ) {
        fragment_queries
            .iter()
            .zip(aggregator.iter_mut())
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

impl ToleranceAdapter<FragmentGroupIndexQuery> for QuadSplittedTransposedIndex {
    fn query_from_elution_group(
        &self,
        tol: &dyn crate::traits::tolerance::Tolerance,
        elution_group: &crate::models::elution_group::ElutionGroup,
    ) -> FragmentGroupIndexQuery {
        self.queries_from_elution_elements_impl(tol, elution_group)
    }
}
