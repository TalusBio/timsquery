use std::sync::Arc;

use timsrust::{AcquisitionType, Frame, MSLevel, QuadrupoleSettings};

use crate::sort_by_indices_multi;
use crate::utils::sorting::argsort_by;

/// A frame after expanding the mobility data and re-sorting it by tof.
pub struct ExpandedFrame {
    pub tof_indices: Vec<u32>,
    pub scan_numbers: Vec<usize>,
    pub intensities: Vec<u32>,
    pub index: usize,
    pub rt: f64,
    pub acquisition_type: AcquisitionType,
    pub ms_level: MSLevel,
    pub quadrupole_settings: Arc<QuadrupoleSettings>,
    pub intensity_correction_factor: f64,
    pub window_group: u8,
}

fn explode_scan_offsets(scan_offsets: Vec<usize>) -> Vec<usize> {
    let num_vals = scan_offsets.last().unwrap();
    let mut exploded_offsets = vec![0; *num_vals];
    for i in 0..scan_offsets.len() - 1 {
        let start = scan_offsets[i];
        let end = scan_offsets[i + 1];
        for j in start..end {
            exploded_offsets[j] = i;
        }
    }
    exploded_offsets
}

impl ExpandedFrame {
    pub fn from_frame(frame: Frame) -> Self {
        let mut tof_indices = frame.tof_indices;
        let mut scan_numbers = explode_scan_offsets(frame.scan_offsets);
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
            index: frame.index,
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
        let exploded_offsets = explode_scan_offsets(scan_offsets);
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
        assert_eq!(expanded_frame.index, 0);
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
