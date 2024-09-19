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
use core::num;
use core::panic;
use log::trace;
use log::{debug, info};

use rayon::prelude::*;

use std::collections::{HashMap, HashSet};
use std::fmt::Display;
use std::sync::Arc;
use timsrust::converters::{
    ConvertableDomain, Frame2RtConverter, Scan2ImConverter, Tof2MzConverter,
};
use timsrust::readers::{FrameReader, MetadataReader};
use timsrust::TimsRustError;
use timsrust::{Frame, QuadrupoleSettings};

// TODO break this module apart ... its getting too big for my taste
// - JSP: 2024-11-19

pub struct QuadSplittedTransposedIndex {
    indices: HashMap<Arc<SingleQuadrupoleSetting>, TransposedQuadIndex>,
    flat_quad_settings: Vec<SingleQuadrupoleSetting>,
    rt_converter: Frame2RtConverter,
    mz_converter: Tof2MzConverter,
    im_converter: Scan2ImConverter,
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
            "QuadSplittedTransposedIndex::query_peaks tof_range: {:?}, scan_range: {:?}, rt_range: {:?}",
            tof_range,
            scan_range,
            rt_range,
        );
        let matching_quads = self.get_matching_quad_settings(precursor_mz_range, scan_range);
        let rt_range = match rt_range {
            Some(x) => Some(x.to_frame_index_range(&self.rt_converter)),
            None => None,
        };

        matching_quads
            .map(move |qs| {
                let tqi = self.indices.get(&qs).unwrap();
                tqi.query_peaks(tof_range, scan_range, rt_range.clone())
            })
            .flatten()
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
                    && (precursor_mz_range.0 <= qs.ranges.isolation_high)
            })
            .filter(move |qs| {
                if let Some(scan_range) = scan_range {
                    (qs.ranges.scan_start <= scan_range.1) && (scan_range.0 <= qs.ranges.scan_end)
                } else {
                    true
                }
            })
            .map(|x| x.clone())
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
        let mobility_index_range = (
            self.im_converter.invert(mobility_range.0) as usize,
            self.im_converter.invert(mobility_range.1) as usize,
        );
        let isolation_mz_range = (
            self.mz_converter
                .invert(precursor_mz_range.0 - quad_range.0 as f64) as f32,
            self.mz_converter
                .invert(precursor_mz_range.1 + quad_range.1 as f64) as f32,
        );
        let frame_index_range = (
            self.rt_converter.invert(rt_range.0) as usize,
            self.rt_converter.invert(rt_range.1) as usize,
        );

        assert!(frame_index_range.0 <= frame_index_range.1);
        assert!(mz_index_range.0 <= mz_index_range.1);
        assert!(isolation_mz_range.0 <= isolation_mz_range.1);
        // assert!(mobility_index_range.0 <= mobility_index_range.1);
        let mobility_index_range = (
            mobility_index_range.0.min(mobility_index_range.1) as usize,
            mobility_index_range.1.max(mobility_index_range.0) as usize,
        );

        let precursor_query = PrecursorIndexQuery {
            frame_index_range,
            mz_index_range,
            mobility_index_range,
            isolation_mz_range,
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
        disp_str.push_str(&"flat_quad_settings: \n");
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
        // 11.305 s on a normal file in my pc ... rlly not bad
        info!("Building transposed quad index from path {}", path);
        let tmp = QuadSplittedTransposedIndexBuilder::from_path(path)?;
        let out = tmp.build();
        info!("Transposed quad index built");
        debug!("{}", out);
        Ok(out)
    }
}

pub struct QuadSplittedTransposedIndexBuilder {
    indices: HashMap<SingleQuadrupoleSetting, TransposedQuadIndexBuilder>,
    rt_converter: Frame2RtConverter,
    mz_converter: Tof2MzConverter,
    im_converter: Scan2ImConverter,
}

impl QuadSplittedTransposedIndexBuilder {
    fn add_frame(&mut self, frame: Frame) {
        let expanded_quad_settings = expand_quad_settings(&frame.quadrupole_settings);

        for qs in expanded_quad_settings {
            // Add key if it doesnt exist ...
            if !self.indices.contains_key(&qs) {
                let max_tof = frame.tof_indices.iter().max().unwrap();
                debug!(
                    "Adding new transposed quad index for qs {:?} with max tof {}",
                    qs, max_tof
                );
                self.indices.insert(
                    qs,
                    TransposedQuadIndexBuilder::new(qs.clone(), *max_tof as usize),
                );
            }
            self.indices
                .get_mut(&qs)
                .unwrap()
                .add_frame_slice(&frame, (qs.ranges.scan_start, qs.ranges.scan_end));
        }
    }

    fn from_path(path: &str) -> Result<Self, TimsRustError> {
        let file_reader = FrameReader::new(&path)?;

        let sql_path = std::path::Path::new(path).join("analysis.tdf");
        let meta_converters = MetadataReader::new(&sql_path)?;

        let mut out = Self {
            indices: HashMap::new(),
            rt_converter: meta_converters.rt_converter,
            mz_converter: meta_converters.mz_converter,
            im_converter: meta_converters.im_converter,
        };
        for frame in file_reader.get_all() {
            let frame = frame?;
            out.add_frame(frame);
        }
        Ok(out)
    }

    fn build(self) -> QuadSplittedTransposedIndex {
        let mut indices = HashMap::new();
        let mut flat_quad_settings = Vec::new();
        for (qs, builder) in self.indices {
            let tmp = Arc::new(qs);
            indices.insert(tmp.clone(), builder.build());
            flat_quad_settings.push(qs.clone());
        }
        flat_quad_settings.sort_by(|a, b| {
            a.ranges
                .isolation_mz
                .partial_cmp(&b.ranges.isolation_mz)
                .unwrap()
        });

        let out = QuadSplittedTransposedIndex {
            indices: indices,
            flat_quad_settings,
            rt_converter: self.rt_converter,
            mz_converter: self.mz_converter,
            im_converter: self.im_converter,
        };

        out
    }
}

#[derive(Debug)]
struct TransposedQuadIndex {
    quad_settings: SingleQuadrupoleSetting,
    frame_index_rt_pairs: Vec<(usize, f64)>,
    peak_buckets: Vec<Option<PeakBucket>>,
    // 72 bytes for peak Option<bucket> (24/vec)
    // Conservatvely ... 30_000_000 elements can be layed out in a vec.
    // 638425 is the max observed val for a tof index ... So we dont need an
    // additional bucketing.

    // In some of my data a lot of the buckets are empty.
    // I could maybe use a sparse representation ... btree? hashmap?
    // num_none: 511271, num_some: 127161
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
            display_opt_peak_bucket(&opt_peak_buckets[tof_with_max as usize])
        ));
    }

    out
}

impl Display for TransposedQuadIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "TransposedQuadIndex\n quad_settings: {}\n frame_index_rt_pairs: {}\n peak_buckets: {}\n",
            self.quad_settings,
            glimpse_vec(&self.frame_index_rt_pairs, Some(GlimpseConfig { max_items: 10, padding: 2, new_line: true })),
            display_opt_peak_bucket_vec(&self.peak_buckets),
        )
    }
}

#[derive(Debug, Clone, Copy)]
pub struct PeakInQuad {
    pub scan_index: usize,
    pub intensity: u32,
    pub frame_index: usize,
    pub tof_index: u32,
}

impl PeakInQuad {
    pub fn from_peak_in_bucket(peak_in_bucket: PeakInBucket, tof_index: u32) -> Self {
        Self {
            scan_index: peak_in_bucket.scan_index,
            intensity: peak_in_bucket.intensity,
            frame_index: peak_in_bucket.local_frame_index as usize,
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

        // Coult I just return an Arc<[intensities]> + ...
        // If I made the peak buckets sparse, I could make it ... not be an option.
        (tof_range.0..tof_range.1)
            .filter(|tof_index| self.peak_buckets[*tof_index as usize].is_some())
            .map(move |tof_index| {
                self.peak_buckets[tof_index as usize]
                    .as_ref()
                    .unwrap()
                    .query_peaks(scan_range, frame_index_range)
                    .map(move |p| PeakInQuad::from_peak_in_bucket(p, tof_index.clone()))
            })
            .flatten()
    }

    fn convert_to_local_frame_range(
        &self,
        rt_range: Option<FrameRTTolerance>,
    ) -> Option<(LocalFrameIndex, LocalFrameIndex)> {
        let frame_index_range: Option<(LocalFrameIndex, LocalFrameIndex)> = match rt_range {
            Some(FrameRTTolerance::Seconds((rt_low, rt_high))) => {
                let frame_id_start = self
                    .frame_index_rt_pairs
                    .binary_search_by(|x| x.1.partial_cmp(&rt_low).unwrap())
                    .unwrap_or_else(|x| x);
                let frame_id_end = self
                    .frame_index_rt_pairs
                    .binary_search_by(|x| x.1.partial_cmp(&rt_high).unwrap())
                    .unwrap_or_else(|x| x);

                Some((
                    frame_id_start as LocalFrameIndex,
                    frame_id_end as LocalFrameIndex,
                ))
            }
            Some(FrameRTTolerance::FrameIndex((frame_low, frame_high))) => {
                let frame_id_start = self
                    .frame_index_rt_pairs
                    .binary_search_by(|x| x.0.partial_cmp(&frame_low).unwrap())
                    .unwrap_or_else(|x| x);

                let frame_id_end = self
                    .frame_index_rt_pairs
                    .binary_search_by(|x| x.0.partial_cmp(&frame_high).unwrap())
                    .unwrap_or_else(|x| x);

                Some((
                    frame_id_start as LocalFrameIndex,
                    frame_id_end as LocalFrameIndex,
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

    pub fn query_tof_index(&self, tof_index: u32) -> Option<&PeakBucket> {
        self.peak_buckets[tof_index as usize].as_ref()
    }
}

struct TransposedQuadIndexBuilder {
    quad_settings: SingleQuadrupoleSetting,
    peak_bucket_builders: Vec<Option<PeakBucketBuilder>>,
    frame_indices: HashMap<usize, u32>,
    frame_index_rt_pairs: Vec<(usize, f64)>,
    num_frames: u32,
    last_frame_index: usize,
}

impl TransposedQuadIndexBuilder {
    fn new(quad_settings: SingleQuadrupoleSetting, init_max_tof: usize) -> Self {
        let mut builder_vec = Vec::with_capacity(init_max_tof);
        for _ in 0..init_max_tof {
            builder_vec.push(None);
        }
        Self {
            quad_settings,
            peak_bucket_builders: builder_vec,
            frame_indices: HashMap::new(),
            frame_index_rt_pairs: Vec::new(),
            num_frames: 0,
            last_frame_index: 0,
        }
    }

    fn add_frame_slice(&mut self, frame: &Frame, scan_range: (usize, usize)) {
        let frame_index = frame.index;
        let frame_rt = frame.rt;
        self.add_frame_index_rt(frame_index, frame_rt);
        for scan_index in scan_range.0..(scan_range.1 - 1) {
            let peak_index_start = frame.scan_offsets[scan_index];
            let peak_index_end = frame.scan_offsets[scan_index + 1];

            for peak_index in peak_index_start..peak_index_end {
                let tof_index = frame.tof_indices[peak_index];
                let intensity = frame.intensities[peak_index];
            self.add_peak(tof_index, scan_index, intensity, frame_index);
            }
        }
    }


    pub fn add_frame_index_rt(&mut self, frame_index: usize, rt: f64) {
        self.num_frames += 1;
        // Make sure they are inserted in order.
        if self.last_frame_index > frame_index {
            panic!("Frame index out of order");
        }
        self.last_frame_index = frame_index;
        self.frame_index_rt_pairs.push((frame_index, rt));

        // Make sure frames are inserted only once.
        if self.frame_indices.contains_key(&frame_index) {
            panic!("Frame index already inserted");
            // TODO make this an error ...
        }
        self.frame_indices.insert(frame_index, self.num_frames);
    }

    fn add_peak(&mut self, tof_index: u32, scan_index: usize, intensity: u32, frame_index: usize) {
        if self.peak_bucket_builders.len() <= (tof_index + 1) as usize {
            let num_to_add = (tof_index as usize - self.peak_bucket_builders.len()) + 20;
            self.peak_bucket_builders
                .extend((0..num_to_add).map(|_| None));
        }
        if self.peak_bucket_builders[tof_index as usize].is_none() {
            self.peak_bucket_builders[tof_index as usize] = Some(PeakBucketBuilder::new(20));
        }
        self.peak_bucket_builders[tof_index as usize]
            .as_mut()
            .unwrap()
            .add_peak(scan_index, intensity, frame_index);
    }

    fn build(self) -> TransposedQuadIndex {
        let peak_buckets = self
            .peak_bucket_builders
            .into_iter()
            .map(|x| match x {
                Some(x) => Some(x.build(&self.frame_indices)),
                None => None,
            })
            .collect();

        TransposedQuadIndex {
            quad_settings: self.quad_settings,
            frame_index_rt_pairs: self.frame_index_rt_pairs,
            peak_buckets,
        }
    }
}

pub struct PeakInBucket {
    pub scan_index: usize,
    pub intensity: u32,
    pub local_frame_index: u32,
}

type LocalFrameIndex = u32;

#[derive(Debug)]
struct PeakBucket {
    intensities: Vec<u32>,
    // TODO: consider using rts instead of local_frame_indices
    local_frame_indices: Vec<LocalFrameIndex>,
    scan_offsets: Vec<usize>,
}

impl Display for PeakBucket {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "PeakBucket: \n    len={},\n    local_frame_indices={},\n    scan_offsets={},\n    intensities={}",
            self.len(),
            glimpse_vec(
                &self.local_frame_indices,
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
            )
        )
    }
}

impl PeakBucket {
    pub fn len(&self) -> usize {
        self.intensities.len()
    }

    pub fn query_peaks(
        &self,
        scan_range: Option<(usize, usize)>,
        local_index_range: Option<(LocalFrameIndex, LocalFrameIndex)>,
    ) -> impl Iterator<Item = PeakInBucket> + '_ {
        let scan_range = match scan_range {
            Some((scan_low, scan_high)) => {
                (scan_low.max(0)..scan_high.min(self.scan_offsets.len() - 1)).into_iter()
            }
            None => 0..self.scan_offsets.len(),
        };

        let tmp = scan_range
            .map(move |scan_index| {
                let peak_index_start = self.scan_offsets[scan_index];
                let peak_index_end = self.scan_offsets[scan_index + 1];

                let local_frame_slice = &self.local_frame_indices[peak_index_start..peak_index_end];
                let intensity_slice = &self.intensities[peak_index_start..peak_index_end];

                local_frame_slice
                    .iter()
                    .zip(intensity_slice)
                    .filter(move |(lfi, _inten)| match local_index_range {
                        Some((low, high)) => &low <= *lfi && *lfi <= &high,
                        None => true,
                    })
                    .map(move |(lfi, inten)| PeakInBucket {
                        scan_index,
                        intensity: *inten,
                        local_frame_index: *lfi,
                    })
            })
            .flatten();
        tmp
    }
}

#[derive(Debug)]
struct PeakBucketBuilder {
    intensities: Vec<u32>,
    frame_indices: Vec<usize>,
    scan_offsets: Vec<usize>,
}

impl PeakBucketBuilder {
    fn new(capacity: usize) -> Self {
        Self {
            intensities: Vec::with_capacity(capacity),
            frame_indices: Vec::with_capacity(capacity),
            scan_offsets: Vec::with_capacity(capacity),
        }
    }

    fn add_peak(&mut self, scan_index: usize, intensity: u32, frame_index: usize) {
        self.intensities.push(intensity);
        self.frame_indices.push(frame_index);
        self.scan_offsets.push(scan_index);
    }

    fn build(mut self, index_mapping: &HashMap<usize, u32>) -> PeakBucket {
        let mut indices = argsort_by(&self.scan_offsets, |x| *x);
        sort_by_indices_multi!(
            &mut indices,
            &mut self.scan_offsets,
            &mut self.frame_indices,
            &mut self.intensities
        );
        let compressed = compress_vec(&self.scan_offsets);
        let local_frame_indices = self
            .frame_indices
            .iter()
            .map(|x| index_mapping[x])
            .collect();

        let out = PeakBucket {
            intensities: self.intensities,
            local_frame_indices: local_frame_indices,
            scan_offsets: compressed,
        };
        debug_assert!(out.verify());
        out
    }
}

impl PeakBucket {
    fn verify(&self) -> bool {
        if self.intensities.len() != self.local_frame_indices.len() {
            return false;
        }
        match self.scan_offsets.last() {
            Some(last) => {
                if *last != self.local_frame_indices.len() {
                    return false;
                }
            }
            None => {
                return false;
            }
        }

        for i in 1..self.scan_offsets.len() {
            if self.scan_offsets[i] < self.scan_offsets[i - 1] {
                return false;
            }
        }

        true
    }
}

impl From<PeakInQuad> for RawPeak {
    fn from(peak_in_quad: PeakInQuad) -> Self {
        RawPeak {
            scan_index: peak_in_quad.scan_index,
            tof_index: peak_in_quad.tof_index,
            intensity: peak_in_quad.intensity,
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
            .par_iter()
            .map(|tof_range| {
                self.query_peaks(
                    *tof_range,
                    precursor_mz_range,
                    scan_range,
                    frame_index_range,
                )
                .map(|peak| RawPeak::from(peak))
                .collect::<Vec<RawPeak>>()
            })
            .flatten()
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
            .par_iter()
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
            .par_iter()
            .zip(aggregator.par_iter_mut())
            .for_each(|(fragment_query, agg)| {
                let precursor_mz_range = (
                    fragment_query.precursor_query.isolation_mz_range.0 as f64,
                    fragment_query.precursor_query.isolation_mz_range.0 as f64,
                );
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
