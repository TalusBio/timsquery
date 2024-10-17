use super::peak_bucket::PeakBucketBuilder;
use super::peak_bucket::{PeakBucket, PeakInBucket};
use crate::models::frames::expanded_frame::ExpandedFrameSlice;
use crate::models::frames::raw_peak::RawPeak;
use crate::models::frames::single_quad_settings::SingleQuadrupoleSetting;
use crate::sort_by_indices_multi;
use crate::utils::display::{glimpse_vec, GlimpseConfig};
use crate::utils::sorting::par_argsort_by;
use log::debug;
use log::info;
use log::trace;
use rayon::prelude::*;
use std::collections::{BTreeMap, HashMap};
use std::fmt::Display;
use std::time::Instant;
use timsrust::converters::{ConvertableDomain, Frame2RtConverter};

#[derive(Debug)]
pub struct TransposedQuadIndex {
    pub quad_settings: Option<SingleQuadrupoleSetting>,
    pub frame_rts: Vec<f64>,
    pub frame_indices: Vec<usize>,
    pub peak_buckets: BTreeMap<u32, PeakBucket>,
    // 72 bytes for peak Option<bucket> (24/vec)
    // Conservatvely ... 30_000_000 elements can be layed out in a vec.
    // 638425 is the max observed val for a tof index ... So we dont need an
    // additional bucketing.
}

fn display_peak_bucket_map(bucket_map: &BTreeMap<u32, PeakBucket>) -> String {
    let mut out = String::new();

    let mut max_peaks = 0;
    let mut tof_with_max = 0;
    let num_buckets = bucket_map.len();
    for (i, peak_bucket) in bucket_map.iter() {
        if max_peaks < peak_bucket.len() {
            max_peaks = peak_bucket.len();
            tof_with_max = *i;
        }
    }

    out.push_str(&format!(
        "PeakBuckets: num_buckets: {}, max_peaks: {}, tof_with_max: {}\n",
        num_buckets, max_peaks, tof_with_max
    ));
    for (i, bucket) in bucket_map.iter().take(3) {
        out.push_str(&format!(" - {}: {}\n", i, bucket));
    }
    out.push_str(&format!(" - ... len = {}\n", bucket_map.len()));

    if tof_with_max > 0 {
        out.push_str(&format!(
            " - Bucket with max tof: {} {}\n",
            tof_with_max, &bucket_map[&tof_with_max]
        ));
    }

    out
}

impl Display for TransposedQuadIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "TransposedQuadIndex\n quad_settings: {:?}\n frame_indices: {}\n frame_rts: {}\n peak_buckets: {}\n",
            self.quad_settings,
            glimpse_vec(&self.frame_indices, Some(GlimpseConfig { max_items: 10, padding: 2, new_line: true })),
            glimpse_vec(&self.frame_rts, Some(GlimpseConfig { max_items: 10, padding: 2, new_line: true })),
            display_peak_bucket_map(&self.peak_buckets),
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
    quad_settings: Option<SingleQuadrupoleSetting>,
    int_slices: Vec<Vec<u32>>,
    tof_slices: Vec<Vec<u32>>,
    scan_slices: Vec<Vec<usize>>,
    frame_indices: Vec<usize>,
    frame_rts: Vec<f64>,
}

impl TransposedQuadIndexBuilder {
    pub fn new(quad_settings: Option<SingleQuadrupoleSetting>) -> Self {
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

    pub fn add_frame_slice(&mut self, slice: ExpandedFrameSlice) {
        self.int_slices.push(slice.intensities);
        self.tof_slices.push(slice.tof_indices);
        self.scan_slices.push(slice.scan_numbers);
        self.frame_indices.push(slice.frame_index);
        self.frame_rts.push(slice.rt);
    }

    pub fn build(self) -> TransposedQuadIndex {
        // TODO: Refactor this function, its getting pretty large.
        let st = Instant::now();
        info!(
            "TransposedQuadIndex::build start ... quad_settings = {:?}",
            self.quad_settings
        );
        let max_tof = *self
            .tof_slices
            .iter()
            .map(|x| x.iter().max().unwrap())
            .max()
            .unwrap();
        let mut tof_counts = vec![0usize; max_tof as usize + 1];
        let mut tot_peaks = 0u64;
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
        let quad_settings = self.quad_settings;

        let aps = Instant::now();

        peak_buckets = if tot_peaks > 10_000_000 {
            self.batched_build_inner(peak_buckets, tot_peaks)
        } else {
            self.build_inner_ref(peak_buckets, tot_peaks)
        };

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

        // FYI par iter here makes no difference.
        let peak_bucket: BTreeMap<u32, PeakBucket> = peak_buckets
            .into_iter()
            .map(|(tof, pb)| (tof, pb.build()))
            .collect();

        let bbe = bb_st.elapsed();
        let elapsed = st.elapsed();
        info!(
            "TransposedQuadIndex::add_frame_slice adding peaks took {:#?} for {} peaks",
            aps.elapsed(),
            tot_peaks
        );
        info!(
            "TransposedQuadIndex::add_frame_slice building buckets took {:#?}",
            bbe
        );
        info!(
            "TransposedQuadIndex::build quad_settings={:?} took {:#?}",
            quad_settings, elapsed
        );
        TransposedQuadIndex {
            quad_settings,
            frame_rts: out_rts,
            frame_indices: out_indices,
            peak_buckets: peak_bucket,
        }
    }

    fn build_inner_ref(
        self,
        mut peak_buckets: HashMap<u32, PeakBucketBuilder>,
        tot_peaks: u64,
    ) -> HashMap<u32, PeakBucketBuilder> {
        let mut added_peaks = 0;

        for slice_ind in 0..self.frame_indices.len() {
            let int_slice = &self.int_slices[slice_ind];
            let tof_slice = &self.tof_slices[slice_ind];
            let scan_slice = &self.scan_slices[slice_ind];
            let frame_rt = self.frame_rts[slice_ind] as f32;

            // Q: is it worth to sort and add peaks in slices? If we do it has to be one level
            // higher, since within each slice, most tof indices would not repeat.
            // A: Nope ... this is faster.
            //
            // In theory there are 3 things to consider:
            // 1. Can we hold a vec that large? (u32::MAX in teory for a 32 bit system)
            //    - For According to the playground ... 268_435_455 is the max number of
            //    u32's that can be allocated ... in some of our data the MS1s have 500M peaks.
            //    so even if we do it, we would need to batch it ...
            //    as a note most MS2s are ~20M
            // 2. Is the sorting slower than querying the hash map that many times?
            //    - Sorting is slower ...
            // 3. I am assuming that sequential import would let us grow each vec only once
            //    per chunk addition instead of once per ... peak in the worst case (although vecs
            //    grow in chunks of +25% ish ...)
            //    - Might be true but we allocate the vecs with capacity, so they dont grow while
            //    we iterate
            //
            for ((inten, tof), scan) in int_slice
                .iter()
                .zip(tof_slice.iter())
                .zip(scan_slice.iter())
            {
                peak_buckets
                    .get_mut(tof)
                    .expect("Should have just added it...")
                    .add_peak(*scan, *inten, frame_rt);

                added_peaks += 1;
                if added_peaks % 500_000 == 0 {
                    debug!(
                        "TransposedQuadIndex::build quad_settings={:?} added_peaks={:?}/{:?}",
                        self.quad_settings, added_peaks, tot_peaks,
                    );
                }
            }
        }

        if added_peaks != tot_peaks {
            println!("TransposedQuadIndex::add_frame_slice failed at peak count check, expected: {}, real: {}", tot_peaks, added_peaks);
            panic!("TransposedQuadIndex::add_frame_slice failed at peak count check");
        }

        peak_buckets
    }

    fn batched_build_inner(
        self,
        mut peak_buckets: HashMap<u32, PeakBucketBuilder>,
        tot_peaks: u64,
    ) -> HashMap<u32, PeakBucketBuilder> {
        let info_prefix = format!("BatchedBuild: quad_settings={:?} ", self.quad_settings);
        info!("{} start", info_prefix);
        let num_slices = self.frame_indices.len();
        let mut added_peaks = 0;

        let mut peaks_in_chunk = 0;
        let mut start = 0;
        let mut end = 0;

        while start < self.frame_indices.len() {
            while peaks_in_chunk < 100_000_000 && end < self.frame_indices.len() {
                peaks_in_chunk += self.int_slices[end].len();
                end += 1;
            }

            let concat_st = Instant::now();
            let mut int_slice = self.int_slices[start..end].concat();
            let mut tof_slice = self.tof_slices[start..end].concat();
            let mut scan_slice = self.scan_slices[start..end].concat();
            let mut rt_slice: Vec<f32> = self.frame_rts[start..end]
                .iter()
                .zip(self.tof_slices[start..end].iter())
                .flat_map(|(rt, tofslice)| vec![*rt as f32; tofslice.len()])
                .collect();

            let concat_elapsed = concat_st.elapsed();

            let sorting_st = Instant::now();
            // let mut indices = argsort_by(&tof_slice, |x| *x);
            let mut indices = par_argsort_by(&tof_slice, |x| *x);

            let argsort_elapsed = sorting_st.elapsed();

            sort_by_indices_multi!(
                &mut indices,
                &mut tof_slice,
                &mut scan_slice,
                &mut int_slice,
                &mut rt_slice
            );

            let sorting_elapsed = sorting_st.elapsed();

            let insertion_st = Instant::now();
            let mut slice_start = 0;
            let mut slice_start_val = tof_slice[0];
            while slice_start < int_slice.len() {
                // // Binary search the current value + 1
                // let slice_end = tof_slice.binary_search(&(slice_start_val + 1));

                // // let local_slice_end = slice_end.unwrap_or_else(|x| x);
                // let local_slice_end = match slice_end {
                //     Ok(x) => {
                //         // If the value is found, we need to walk back to find the first instance
                //         // of the value ...
                //         let mut local_slice_end = x;
                //         while tof_slice[local_slice_end - 1] > slice_start_val {
                //             local_slice_end -= 1;
                //         }
                //         local_slice_end
                //     }
                //     Err(x) => x,
                // };
                //
                let mut local_slice_end = slice_start;
                while tof_slice[local_slice_end] == slice_start_val {
                    local_slice_end += 1;
                    // There has to be a way to add this to the conditional above ...
                    // But alas ... this works ...
                    if local_slice_end == tof_slice.len() {
                        break;
                    }
                }

                let range_use = slice_start..local_slice_end;

                peak_buckets
                    .get_mut(&slice_start_val)
                    .expect("Tof should have been added during the build")
                    .extend_peaks(
                        &scan_slice[range_use.clone()],
                        &int_slice[range_use.clone()],
                        &rt_slice[range_use.clone()],
                    );

                added_peaks += range_use.len() as u64;

                assert!(slice_start_val == tof_slice[slice_start]);

                if local_slice_end == tof_slice.len() {
                    break;
                }
                assert!(slice_start_val < tof_slice[local_slice_end]);

                slice_start = local_slice_end;
                slice_start_val = tof_slice[slice_start];
            }

            let insertion_elapsed = insertion_st.elapsed();
            info!(
                "BatchedBuild: quad_settings={:?} start={:?} end={:?}/{} peaks {}/{} concat took {:#?} sorting took arg: {:#?} / {:#?} insertion took {:#?}",
                self.quad_settings, start, end, num_slices, added_peaks, tot_peaks, concat_elapsed, argsort_elapsed, sorting_elapsed, insertion_elapsed,
            );
            start = end;
            peaks_in_chunk = 0;
        }

        if added_peaks != tot_peaks {
            println!("TransposedQuadIndex::add_frame_slice failed at peak count check, expected: {}, real: {}", tot_peaks, added_peaks);
            panic!("TransposedQuadIndex::add_frame_slice failed at peak count check");
        }

        peak_buckets
    }
}
