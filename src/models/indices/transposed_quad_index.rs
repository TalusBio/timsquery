use crate::models::frames::single_quad_settings::expand_quad_settings;
use crate::models::frames::single_quad_settings::SingleQuadrupoleSetting;
use crate::sort_by_indices_multi;
use crate::utils::compress_explode::{compress_vec, explode_vec};
use crate::utils::sorting::argsort_by;
use core::panic;
use std::collections::{HashMap, HashSet};
use timsrust::{Frame, QuadrupoleSettings};

struct QuadSplittedTransposedIndex {
    indices: HashMap<SingleQuadrupoleSetting, TransposedQuadIndex>,
}

struct QuadSplittedTransposedIndexBuilder {
    indices: HashMap<SingleQuadrupoleSetting, TransposedQuadIndexBuilder>,
}

impl QuadSplittedTransposedIndexBuilder {
    fn new() -> Self {
        Self {
            indices: HashMap::new(),
        }
    }

    fn add_frame(&mut self, frame: Frame) {
        let frame_index = frame.index;
        let expanded_quad_settings = expand_quad_settings(&frame.quadrupole_settings);

        for qs in expanded_quad_settings {
            if self.indices.contains_key(&qs) {
                self.indices
                    .get_mut(&qs)
                    .unwrap()
                    .add_frame_slice(&frame, (qs.ranges.scan_start, qs.ranges.scan_end));
            }
        }
    }
}

struct TransposedQuadIndex {
    quad_settings: SingleQuadrupoleSetting,
    frame_index_rt_pairs: Vec<(usize, f64)>,
    peak_buckets: Vec<Option<PeakBucket>>,
    // 72 bytes for peak Option<bucket> (24/vec)
    // Conservatvely ... 30_000_000 elements can be layed out in a vec.
    // 638425 is the max observed val for a tof index ... So we dont need an
    // additional bucketing.
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
        if self.peak_bucket_builders[tof_index as usize].is_none() {
            self.peak_bucket_builders[tof_index as usize] =
                Some(PeakBucketBuilder::new(tof_index as usize));
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

struct PeakBucket {
    intensities: Vec<u32>,
    local_frame_indices: Vec<u32>,
    scan_offsets: Vec<usize>,
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
