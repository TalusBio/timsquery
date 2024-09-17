use crate::models::frames::single_quad_settings::expand_quad_settings;
use crate::models::frames::single_quad_settings::SingleQuadrupoleSetting;
use crate::sort_by_indices_multi;
use crate::utils::compress_explode::{compress_vec, explode_vec};
use crate::utils::display::glimpse_vec;
use crate::utils::sorting::argsort_by;
use core::num;
use core::panic;
use log::{debug, info};
use std::collections::{HashMap, HashSet};
use std::fmt::Display;
use timsrust::readers::{FrameReader, MetadataReader};
use timsrust::TimsRustError;
use timsrust::{Frame, QuadrupoleSettings};

pub struct QuadSplittedTransposedIndex {
    indices: HashMap<SingleQuadrupoleSetting, TransposedQuadIndex>,
}

impl Display for QuadSplittedTransposedIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut disp_str = String::new();
        disp_str.push_str("QuadSplittedTransposedIndex\n");
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
}

impl QuadSplittedTransposedIndexBuilder {
    fn new() -> Self {
        Self {
            indices: HashMap::new(),
        }
    }

    fn add_frame(&mut self, frame: Frame) {
        let expanded_quad_settings = expand_quad_settings(&frame.quadrupole_settings);

        for qs in expanded_quad_settings {
            // Add key if it doesnt exist ...
            if !self.indices.contains_key(&qs) {
                let max_tof = frame.tof_indices.iter().max().unwrap();
                info!(
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

        let mut out = Self::new();
        for frame in file_reader.get_all() {
            let frame = frame?;
            out.add_frame(frame);
        }
        Ok(out)
    }

    fn build(self) -> QuadSplittedTransposedIndex {
        let mut out = QuadSplittedTransposedIndex {
            indices: HashMap::new(),
        };
        for (qs, builder) in self.indices {
            out.indices.insert(qs, builder.build());
        }
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

    out.push_str(&format!(
        "PeakBuckets: num_none: {}, num_some: {}\n",
        num_none, num_some
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

    out
}

impl Display for TransposedQuadIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "TransposedQuadIndex\n quad_settings: {}\n frame_index_rt_pairs: {}\n peak_buckets: {}\n",
            self.quad_settings,
            glimpse_vec(&self.frame_index_rt_pairs),
            display_opt_peak_bucket_vec(&self.peak_buckets),
        )
    }
}

struct PeakInQuad {
    scan_index: usize,
    intensity: u32,
    frame_index: usize,
    tof_index: u32,
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
        rt_range: Option<(f64, f64)>,
        // ) -> Vec<PeakInQuad> {
    ) -> impl Iterator<Item = PeakInQuad> + '_ {
        let frame_index_range: Option<(usize, usize)> = match rt_range {
            Some((rt_low, rt_high)) => {
                let frame_id_start = match self
                    .frame_index_rt_pairs
                    .binary_search_by(|x| x.1.partial_cmp(&rt_low).unwrap())
                {
                    Ok(x) => self.frame_index_rt_pairs[x].0,
                    Err(x) => self.frame_index_rt_pairs[x].0,
                };
                let frame_id_end = match self
                    .frame_index_rt_pairs
                    .binary_search_by(|x| x.1.partial_cmp(&rt_high).unwrap())
                {
                    Ok(x) => self.frame_index_rt_pairs[x].0,
                    Err(x) => self.frame_index_rt_pairs[x].0,
                };
                Some((frame_id_start, frame_id_end))
            }
            None => None,
        };

        // TODO reimplement as an iterator ...
        // This version is not compatible with the borrow checker unless I collect the vec...
        // which will do for now for prototyping.

        // Coult I just return an Arc<[intensities]> + ...
        (tof_range.0..tof_range.1)
            .filter(|tof_index| self.peak_buckets[*tof_index as usize].is_some())
            .map(move |tof_index| {
                self.peak_buckets[tof_index as usize]
                    .as_ref()
                    .unwrap()
                    .query_peaks(scan_range)
                    .filter(|p| match frame_index_range {
                        Some((frame_id_start, frame_id_end)) => {
                            p.local_frame_index >= frame_id_start as u32
                                && p.local_frame_index < frame_id_end as u32
                        }
                        None => true,
                    })
                    .map(|p| PeakInQuad::from_peak_in_bucket(p, tof_index.clone()))
                    .collect::<Vec<PeakInQuad>>()
            })
            .flatten()
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
        for scan_index in scan_range.0..scan_range.1 {
            let tof_index = frame.tof_indices[scan_index];
            let intensity = frame.intensities[scan_index];
            self.add_peak(tof_index, scan_index, intensity, frame_index);
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
        while self.peak_bucket_builders.len() <= tof_index as usize {
            self.peak_bucket_builders.push(None);
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

struct PeakInBucket {
    scan_index: usize,
    intensity: u32,
    local_frame_index: u32,
}

#[derive(Debug)]
struct PeakBucket {
    intensities: Vec<u32>,
    local_frame_indices: Vec<u32>,
    scan_offsets: Vec<usize>,
}

impl Display for PeakBucket {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "PeakBucket: len={}, local_frame_indices={}, scan_offsets={}, intensities={}",
            self.len(),
            glimpse_vec(&self.local_frame_indices),
            glimpse_vec(&self.scan_offsets),
            glimpse_vec(&self.intensities)
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
    ) -> impl Iterator<Item = PeakInBucket> + '_ {
        let scan_range = match scan_range {
            Some((scan_low, scan_high)) => (scan_low..scan_high).into_iter(),
            None => 0..self.len(),
        };

        let tmp = scan_range
            .map(move |scan_index| {
                let peak_index_start = self.scan_offsets[scan_index];
                let peak_index_end = self.scan_offsets[scan_index + 1];

                (peak_index_start..peak_index_end)
                    .into_iter()
                    .map(move |peak_index| {
                        let local_frame_index = self.local_frame_indices[peak_index];
                        let intensity = self.intensities[peak_index];
                        PeakInBucket {
                            scan_index,
                            intensity,
                            local_frame_index,
                        }
                    })
            })
            .flatten();
        tmp
    }
}

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
