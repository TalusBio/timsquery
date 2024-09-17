use std::sync::Arc;

use super::single_quad_settings::{expand_quad_settings, SingleQuadrupoleSetting};
use timsrust::{AcquisitionType, Frame, MSLevel, QuadrupoleSettings};

use crate::sort_by_indices_multi;
use crate::utils::sorting::argsort_by;

/// A frame after expanding the mobility data and re-sorting it by tof.
pub struct ExpandedFrame {
    pub tof_indices: Vec<u32>,
    pub scan_numbers: Vec<usize>,
    pub intensities: Vec<u32>,
    pub frame_index: usize,
    pub rt: f64,
    pub acquisition_type: AcquisitionType,
    pub ms_level: MSLevel,
    pub quadrupole_settings: Arc<QuadrupoleSettings>,
    pub intensity_correction_factor: f64,
    pub window_group: u8,
}

pub struct ExpandedFrameSlice {
    pub tof_indices: Vec<u32>, // I could Arc<[u32]> if I didnt want to sort by it ...
    pub scan_numbers: Vec<usize>,
    pub intensities: Vec<u32>,
    pub frame_index: usize,
    pub rt: f64,
    pub acquisition_type: AcquisitionType,
    pub ms_level: MSLevel,
    pub quadrupole_settings: SingleQuadrupoleSetting,
    pub intensity_correction_factor: f64,
    pub window_group: u8,
    pub window_subindex: u8,
}

impl ExpandedFrameSlice {
    pub fn sort_by_tof(&mut self) {
        let mut indices = argsort_by(&self.tof_indices, |x| *x);
        sort_by_indices_multi!(
            &mut indices,
            &mut self.tof_indices,
            &mut self.scan_numbers,
            &mut self.intensities
        );
    }
}

fn expand_and_split_frame(frame: &Frame) -> Vec<ExpandedFrameSlice> {
    let quad_settings = &frame.quadrupole_settings;
    let expanded_quad_settings = expand_quad_settings(quad_settings);

    let mut out = Vec::with_capacity(quad_settings.scan_ends.len());
    for (i, qs) in expanded_quad_settings.into_iter().enumerate() {
        // This whole block is kind of ugly but I think its good enough for now ...
        let slice_scan_start = qs.ranges.scan_start;
        let slice_scan_end = qs.ranges.scan_end;
        let tof_index_index_slice_start = frame.scan_offsets[slice_scan_start];
        let tof_index_index_slice_end = frame.scan_offsets[slice_scan_end];

        let mut slice_tof_indices =
            frame.tof_indices[tof_index_index_slice_start..tof_index_index_slice_end].to_vec();
        let mut slice_scan_numbers = explode_scan_offsets(
            &frame.scan_offsets[slice_scan_start..slice_scan_end],
            slice_scan_start,
        );
        let mut slice_intensities =
            frame.intensities[tof_index_index_slice_start..tof_index_index_slice_end].to_vec();

        let mut indices = argsort_by(&slice_tof_indices, |x| *x);
        sort_by_indices_multi!(
            &mut indices,
            &mut slice_tof_indices,
            &mut slice_scan_numbers,
            &mut slice_intensities
        );

        let mut curr_slice = ExpandedFrameSlice {
            tof_indices: slice_tof_indices,
            scan_numbers: slice_scan_numbers,
            intensities: slice_intensities,
            frame_index: frame.index,
            rt: frame.rt,
            acquisition_type: frame.acquisition_type,
            ms_level: frame.ms_level,
            quadrupole_settings: qs,
            intensity_correction_factor: frame.intensity_correction_factor,
            window_group: frame.window_group,
            window_subindex: i as u8,
        };
        curr_slice.sort_by_tof();
        out.push(curr_slice);
    }
    out
}

fn explode_scan_offsets(scan_offsets: &[usize], offset_val: usize) -> Vec<usize> {
    let num_vals = scan_offsets.last().unwrap();
    let mut exploded_offsets = vec![0; *num_vals];
    for i in 0..scan_offsets.len() - 1 {
        let start = scan_offsets[i];
        let end = scan_offsets[i + 1];
        for j in start..end {
            exploded_offsets[j] = i + offset_val;
        }
    }
    exploded_offsets
}

impl ExpandedFrame {
    pub fn from_frame(frame: Frame) -> Self {
        let mut tof_indices = frame.tof_indices;
        let mut scan_numbers = explode_scan_offsets(&frame.scan_offsets, 0);
        let mut intensities = frame.intensities;

        let mut indices = argsort_by(&tof_indices, |x| *x);
        sort_by_indices_multi!(
            &mut indices,
            &mut tof_indices,
            &mut scan_numbers,
            &mut intensities
        );

        ExpandedFrame {
            tof_indices: tof_indices,
            scan_numbers: scan_numbers,
            intensities: intensities,
            frame_index: frame.index,
            rt: frame.rt,
            acquisition_type: frame.acquisition_type,
            ms_level: frame.ms_level,
            quadrupole_settings: frame.quadrupole_settings,
            intensity_correction_factor: frame.intensity_correction_factor,
            window_group: frame.window_group,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_explode_scan_offsets() {
        let scan_offsets = vec![0, 2, 4, 8];
        let exploded_offsets = explode_scan_offsets(&scan_offsets, 0);
        assert_eq!(exploded_offsets, vec![0, 0, 1, 1, 2, 2, 2, 2]);
    }

    #[test]
    fn test_expanded_frame_from_frame() {
        let frame = Frame {
            tof_indices: vec![1, 2, 3, 4, 1, 2],
            scan_offsets: vec![0, 2, 4, 6],
            intensities: vec![10, 20, 30, 40, 50, 60],
            index: 0,
            rt: 0.0,
            acquisition_type: AcquisitionType::DIAPASEF,
            ms_level: MSLevel::MS1,
            quadrupole_settings: Arc::new(QuadrupoleSettings::default()),
            intensity_correction_factor: 1.0,
            window_group: 0,
        };

        let expanded_frame = ExpandedFrame::from_frame(frame);
        assert_eq!(expanded_frame.tof_indices, vec![1, 1, 2, 2, 3, 4]);
        assert_eq!(expanded_frame.scan_numbers, vec![0, 2, 0, 2, 1, 1]);
        assert_eq!(expanded_frame.intensities, vec![10, 50, 20, 60, 30, 40]);
        assert_eq!(expanded_frame.frame_index, 0);
        assert_eq!(expanded_frame.rt, 0.0);
        assert_eq!(expanded_frame.acquisition_type, AcquisitionType::DIAPASEF);
        assert_eq!(expanded_frame.ms_level, MSLevel::MS1);
        assert_eq!(
            *expanded_frame.quadrupole_settings,
            QuadrupoleSettings::default()
        );
        assert_eq!(expanded_frame.intensity_correction_factor, 1.0);
        assert_eq!(expanded_frame.window_group, 0);
    }
}
