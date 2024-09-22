use super::peak_bucket::PeakBucketBuilder;
use super::peak_bucket::{PeakBucket, PeakInBucket};
use crate::models::frames::raw_peak::RawPeak;
use crate::models::frames::single_quad_settings::SingleQuadrupoleSetting;
use crate::utils::display::{glimpse_vec, GlimpseConfig};

use log::info;
use log::trace;
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap};
use std::fmt::Display;
use std::time::Instant;
use timsrust::converters::{ConvertableDomain, Frame2RtConverter};

#[derive(Debug)]
pub struct TransposedQuadIndex {
    pub quad_settings: SingleQuadrupoleSetting,
    pub frame_rts: Vec<f64>,
    pub frame_indices: Vec<usize>,
    pub peak_buckets: BTreeMap<u32, PeakBucket>,
    // 72 bytes for peak Option<bucket> (24/vec)
    // Conservatvely ... 30_000_000 elements can be layed out in a vec.
    // 638425 is the max observed val for a tof index ... So we dont need an
    // additional bucketing.
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
        // TODO consider if I should allow only RT here, since it would in theory
        // force me to to the repreatable work beforehand.
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

#[derive(Debug, Clone, Copy)]
pub struct PeakInQuad {
    pub scan_index: usize,
    pub intensity: u32,
    pub retention_time: f32,
    pub tof_index: u32,
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

// Q: Do I use this? Should I just have it as the frame index as
// an option?
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

#[derive(Debug, Clone)]
pub struct TransposedQuadIndexBuilder {
    quad_settings: SingleQuadrupoleSetting,
    int_slices: Vec<Vec<u32>>,
    tof_slices: Vec<Vec<u32>>,
    scan_slices: Vec<Vec<usize>>,
    frame_indices: Vec<usize>,
    frame_rts: Vec<f64>,
}

impl TransposedQuadIndexBuilder {
    pub fn new(quad_settings: SingleQuadrupoleSetting) -> Self {
        Self {
            quad_settings,
            int_slices: Vec::new(),
            tof_slices: Vec::new(),
            scan_slices: Vec::new(),
            frame_indices: Vec::new(),
            frame_rts: Vec::new(),
        }
    }

    pub fn fold(&mut self, other: Self) {
        self.int_slices.extend(other.int_slices);
        self.tof_slices.extend(other.tof_slices);
        self.scan_slices.extend(other.scan_slices);
        self.frame_indices.extend(other.frame_indices);
        self.frame_rts.extend(other.frame_rts);
    }

    pub fn add_frame_slice(
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

    pub fn build(self) -> TransposedQuadIndex {
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

        let aps = Instant::now();

        for slice_ind in 0..self.frame_indices.len() {
            let int_slice = &self.int_slices[slice_ind];
            let tof_slice = &self.tof_slices[slice_ind];
            let scan_slice = &self.scan_slices[slice_ind];
            // let frame_index = self.frame_indices[slice_ind];
            let frame_rt = self.frame_rts[slice_ind] as f32;

            // Q: is it worth to sort and add peaks in slices?

            for ((inten, tof), scan) in int_slice
                .into_iter()
                .zip(tof_slice.into_iter())
                .zip(scan_slice.into_iter())
            {
                peak_buckets
                    .get_mut(tof)
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
                let real_count = curr_bucket.len();
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
