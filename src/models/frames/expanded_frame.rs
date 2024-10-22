use std::collections::HashMap;
use std::sync::Arc;

use super::single_quad_settings::{
    expand_quad_settings, ExpandedFrameQuadSettings, SingleQuadrupoleSetting,
};
use itertools::Itertools;
use rayon::iter::IntoParallelIterator;
use rayon::prelude::*;
use timsrust::converters::{Scan2ImConverter, Tof2MzConverter};
use timsrust::{AcquisitionType, Frame, MSLevel, QuadrupoleSettings};

use crate::sort_by_indices_multi;
use crate::utils::compress_explode::explode_vec;
use crate::utils::frame_processing::lazy_centroid_weighted_frame;
use crate::utils::sorting::argsort_by;
use crate::utils::tolerance_ranges::{scan_tol_range, tof_tol_range};
use log::debug;

/// A frame after expanding the mobility data and re-sorting it by tof.
#[derive(Debug, Clone)]
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

#[derive(Debug, Clone)]
pub struct ExpandedFrameSlice {
    pub tof_indices: Vec<u32>, // I could Arc<[u32]> if I didnt want to sort by it ...
    pub scan_numbers: Vec<usize>,
    pub intensities: Vec<u32>,
    pub frame_index: usize,
    pub rt: f64,
    pub acquisition_type: AcquisitionType,
    pub ms_level: MSLevel,
    pub quadrupole_settings: Option<SingleQuadrupoleSetting>,
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

    pub fn len(&self) -> usize {
        self.tof_indices.len()
    }

    pub fn is_empty(&self) -> bool {
        self.tof_indices.is_empty()
    }
}

fn trim_scan_edges(scan_start: usize, scan_end: usize) -> (usize, usize) {
    // Poor man's fix to different quad windows bleeding onto each other ...
    // I will trim the smallest of 20 scans or 10% of the range size.
    //
    // old dumber implementation
    // let peak_start = frame.scan_offsets[qs.ranges.scan_start + 10];
    // let peak_end = frame.scan_offsets[qs.ranges.scan_end - 10];
    // (peak_start, peak_end)
    //
    let scan_range = scan_end - scan_start;
    assert!(scan_range > 0, "Expected scan range to be positive...");
    let offset_use = (scan_range / 10).min(20);

    let new_start = scan_start + offset_use;
    let new_end = scan_end - offset_use;

    assert!(new_start < new_end, "Expected new_start < new_end");

    (new_start, new_end)
}

fn expand_unfragmented_frame(frame: Frame) -> ExpandedFrameSlice {
    let scan_numbers = explode_vec(&frame.scan_offsets);
    let intensities = frame.intensities;
    let tof_indices = frame.tof_indices;
    let mut curr_slice = ExpandedFrameSlice {
        tof_indices,
        scan_numbers,
        intensities,
        frame_index: frame.index,
        rt: frame.rt,
        acquisition_type: frame.acquisition_type,
        ms_level: frame.ms_level,
        quadrupole_settings: None,
        intensity_correction_factor: frame.intensity_correction_factor,
        window_group: frame.window_group,
        window_subindex: 0u8,
    };
    curr_slice.sort_by_tof();
    curr_slice
}

fn expand_fragmented_frame(
    frame: Frame,
    quads: Vec<SingleQuadrupoleSetting>,
) -> Vec<ExpandedFrameSlice> {
    let mut out = Vec::with_capacity(quads.len());
    let exploded_scans = explode_vec(&frame.scan_offsets);
    for (i, qs) in quads.into_iter().enumerate() {
        // This whole block is kind of ugly but I think its good enough for now ...
        let (slice_scan_start, slice_scan_end) =
            trim_scan_edges(qs.ranges.scan_start, qs.ranges.scan_end);

        let tof_index_index_slice_start = frame.scan_offsets[slice_scan_start];
        let tof_index_index_slice_end = frame.scan_offsets[slice_scan_end];

        let slice_tof_indices =
            frame.tof_indices[tof_index_index_slice_start..tof_index_index_slice_end].to_vec();
        let slice_scan_numbers =
            exploded_scans[tof_index_index_slice_start..tof_index_index_slice_end].to_vec();
        let slice_intensities =
            frame.intensities[tof_index_index_slice_start..tof_index_index_slice_end].to_vec();

        let mut curr_slice = ExpandedFrameSlice {
            tof_indices: slice_tof_indices,
            scan_numbers: slice_scan_numbers,
            intensities: slice_intensities,
            frame_index: frame.index,
            rt: frame.rt,
            acquisition_type: frame.acquisition_type,
            ms_level: frame.ms_level,
            quadrupole_settings: Some(qs),
            intensity_correction_factor: frame.intensity_correction_factor,
            window_group: frame.window_group,
            window_subindex: i as u8,
        };
        // Q: is this sorting twice? since sort before creating the expanded frame slices.
        // TODO: Use the state type pattern to make sure only one sort happens.
        curr_slice.sort_by_tof();
        out.push(curr_slice);
    }
    out
}

pub fn expand_and_split_frame(frame: Frame) -> Vec<ExpandedFrameSlice> {
    let quad_settings = &frame.quadrupole_settings;
    // Q: Can I save a vec allocation if I make the expanded quad_settings
    // To return a vec of builders? or Enum(settings|frame_slice)?
    let expanded_quad_settings = expand_quad_settings(quad_settings);

    match expanded_quad_settings {
        ExpandedFrameQuadSettings::Unfragmented => vec![expand_unfragmented_frame(frame)],
        ExpandedFrameQuadSettings::Fragmented(quads) => expand_fragmented_frame(frame, quads),
    }
}

impl ExpandedFrame {
    pub fn from_frame(frame: Frame) -> Self {
        let mut tof_indices = frame.tof_indices;
        let mut scan_numbers = explode_vec(&frame.scan_offsets);
        let mut intensities = frame.intensities;

        let mut indices = argsort_by(&tof_indices, |x| *x);
        sort_by_indices_multi!(
            &mut indices,
            &mut tof_indices,
            &mut scan_numbers,
            &mut intensities
        );

        ExpandedFrame {
            tof_indices,
            scan_numbers,
            intensities,
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

pub fn expand_and_arrange_frames(
    frames: Vec<Frame>,
) -> HashMap<Option<SingleQuadrupoleSetting>, Vec<ExpandedFrameSlice>> {
    let mut out = HashMap::new();
    let split: Vec<ExpandedFrameSlice> = frames
        .into_par_iter()
        .flat_map(|frame| expand_and_split_frame(frame))
        .collect();
    for es in split {
        out.entry(es.quadrupole_settings)
            .or_insert(Vec::new())
            .push(es);
    }

    // Finally sort each of them internally by retention time.
    for (_, es) in out.iter_mut() {
        es.sort_unstable_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap());
    }
    out
}

pub fn par_expand_and_centroid_frames(
    frames: Vec<Frame>,
    ims_tol_pct: f64,
    mz_tol_ppm: f64,
    ims_converter: &Scan2ImConverter,
    mz_converter: &Tof2MzConverter,
) -> HashMap<Option<SingleQuadrupoleSetting>, Vec<ExpandedFrameSlice>> {
    let split_frames = expand_and_arrange_frames(frames);

    let out: HashMap<Option<SingleQuadrupoleSetting>, Vec<ExpandedFrameSlice>> = split_frames
        .into_iter()
        .map(|(qs, frameslices)| {
            // NOTE: Since the the centroiding runs in paralel over the windows, its ok if this
            // outer loop is done in series.
            let start_peaks: usize = frameslices.iter().map(|x| x.len()).sum();
            let centroided = par_lazy_centroid_frameslices(
                &frameslices,
                3,
                ims_tol_pct,
                mz_tol_ppm,
                ims_converter,
                mz_converter,
            );
            let end_peaks: usize = centroided.iter().map(|x| x.len()).sum();
            debug!(
                "Peak counts for quad {:?}: raw={}/centroid={}",
                qs, start_peaks, end_peaks
            );
            (qs, centroided)
        })
        .collect();

    out
}

fn centroid_frameslice_window(
    frameslices: &[ExpandedFrameSlice],
    ims_tol_pct: f64,
    mz_tol_ppm: f64,
    ims_converter: &Scan2ImConverter,
    mz_converter: &Tof2MzConverter,
) -> ExpandedFrameSlice {
    assert!(frameslices.len() > 1, "Expected at least 2 frameslices");
    let reference_index = frameslices.len() / 2;

    // this is A LOT of cloning ... look into whether it is needed.

    let mut tof_array: Vec<u32> = frameslices
        .iter()
        .flat_map(|x| x.tof_indices.clone())
        .collect();

    let mut ims_array: Vec<usize> = frameslices
        .iter()
        .flat_map(|x| x.scan_numbers.clone())
        .collect();
    let mut weight_array: Vec<u32> = frameslices
        .iter()
        .flat_map(|x| x.intensities.clone())
        .collect();
    let mut intensity_array: Vec<u32> = frameslices
        .iter()
        .enumerate()
        .flat_map(|(i, x)| {
            if i == reference_index {
                x.intensities.clone()
            } else {
                vec![0u32; x.tof_indices.len()]
            }
        })
        .collect();

    let mut tof_order = argsort_by(&tof_array, |x| *x);
    sort_by_indices_multi!(
        &mut tof_order,
        &mut tof_array,
        &mut ims_array,
        &mut weight_array,
        &mut intensity_array
    );

    let ((mzs, intensities), imss) = lazy_centroid_weighted_frame(
        &tof_array,
        &ims_array,
        &weight_array,
        &intensity_array,
        |&tof| tof_tol_range(tof, mz_tol_ppm, mz_converter),
        |&scan| scan_tol_range(scan, ims_tol_pct, ims_converter),
    );

    ExpandedFrameSlice {
        tof_indices: mzs,
        scan_numbers: imss,
        intensities,
        frame_index: frameslices[reference_index].frame_index,
        rt: frameslices[reference_index].rt,
        acquisition_type: frameslices[reference_index].acquisition_type,
        ms_level: frameslices[reference_index].ms_level,
        quadrupole_settings: frameslices[reference_index].quadrupole_settings,
        intensity_correction_factor: frameslices[reference_index].intensity_correction_factor,
        window_group: frameslices[reference_index].window_group,
        window_subindex: frameslices[reference_index].window_subindex,
    }
}

pub fn par_lazy_centroid_frameslices(
    frameslices: &[ExpandedFrameSlice],
    window_width: usize,
    ims_tol_pct: f64,
    mz_tol_ppm: f64,
    ims_converter: &Scan2ImConverter,
    mz_converter: &Tof2MzConverter,
) -> Vec<ExpandedFrameSlice> {
    assert!(
        frameslices
            .iter()
            .tuple_windows()
            .map(|(a, b)| { a.rt < b.rt && a.quadrupole_settings == b.quadrupole_settings })
            .all(|x| x),
        "All frames should be sorted by rt and have the same quad settings"
    );

    assert!(frameslices.len() > window_width);

    let local_lambda = |fss: &[ExpandedFrameSlice]| {
        centroid_frameslice_window(fss, ims_tol_pct, mz_tol_ppm, ims_converter, mz_converter)
    };

    frameslices
        .par_windows(window_width)
        .map(|window| local_lambda(window))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

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
