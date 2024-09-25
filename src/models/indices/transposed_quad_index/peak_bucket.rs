use crate::sort_by_indices_multi;
use crate::utils::compress_explode::compress_vec;
use crate::utils::display::{glimpse_vec, GlimpseConfig};
use crate::utils::sorting::argsort_by;
use std::fmt::Display;

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
pub struct PeakBucket {
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
pub struct PeakBucketBuilder {
    intensities: Vec<u32>,
    retention_times: Vec<f32>,
    scan_offsets: Vec<usize>,
    tof_index: u32,
}

impl PeakBucketBuilder {
    pub fn new(capacity: usize, tof_index: u32) -> Self {
        Self {
            intensities: Vec::with_capacity(capacity),
            retention_times: Vec::with_capacity(capacity),
            scan_offsets: Vec::with_capacity(capacity),
            tof_index,
        }
    }

    pub fn len(&self) -> usize {
        self.intensities.len()
    }

    pub fn is_empty(&self) -> bool {
        self.intensities.is_empty()
    }

    pub fn add_peak(&mut self, scan_index: usize, intensity: u32, retention_time: f32) {
        self.intensities.push(intensity);
        self.retention_times.push(retention_time);
        self.scan_offsets.push(scan_index);
    }

    pub fn build(mut self) -> PeakBucket {
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

        debug_assert!(out.verify(), "PeakBucket::build failed at verify");
        out
    }
}

impl PeakBucket {
    pub fn len(&self) -> usize {
        self.intensities.len()
    }

    pub fn is_empty(&self) -> bool {
        self.intensities.is_empty()
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

        (start_min..end_max).filter_map(move |x| {
            let scan_index = self.scan_offsets[x];
            let retention_time = self.retention_times[x];
            if let Some((low, high)) = rt_range {
                if retention_time < low || retention_time > high {
                    return None;
                }
            }
            Some(PeakInBucket {
                scan_index,
                intensity: self.intensities[x],
                retention_time,
            })
        })
    }

    pub fn query_peaks_compressed(
        &self,
        scan_range: Option<(usize, usize)>,
        rt_range: Option<(f32, f32)>,
    ) -> impl Iterator<Item = PeakInBucket> + '_ {
        let scan_range = match scan_range {
            Some((scan_low, scan_high)) => scan_low..scan_high.min(self.scan_offsets.len() - 1),
            None => 0..self.scan_offsets.len(),
        };
        scan_range
            .flat_map(move |scan_index| {
                let peak_index_start = self.scan_offsets[scan_index];
                let peak_index_end = self.scan_offsets[scan_index + 1];

                (peak_index_start..peak_index_end).map(move |peak_index| {
                    let retention_time = self.retention_times[peak_index];

                    if let Some((low, high)) = rt_range {
                        if retention_time < low || retention_time > high {
                            return None;
                        }
                    }

                    let intensity = self.intensities[peak_index];
                    Some(PeakInBucket {
                        scan_index,
                        intensity,
                        retention_time,
                    })
                })
            })
            .flatten()
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
