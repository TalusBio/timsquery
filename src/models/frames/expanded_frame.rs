use super::single_quad_settings::{
    expand_quad_settings, ExpandedFrameQuadSettings, SingleQuadrupoleSetting,
};
use itertools::Itertools;
use rayon::prelude::*;
use std::collections::HashMap;
use std::marker::PhantomData;
use std::sync::Arc;
use std::time::Instant;
use timsrust::converters::{Scan2ImConverter, Tof2MzConverter};
use timsrust::{AcquisitionType, Frame, MSLevel, QuadrupoleSettings};

use super::peak_in_quad::PeakInQuad;
use crate::sort_by_indices_multi;
use crate::sort_vecs_by_first;
use crate::utils::compress_explode::explode_vec;
use crate::utils::frame_processing::{lazy_centroid_weighted_frame, PeakArrayRefs};
use crate::utils::sorting::argsort_by;
use crate::utils::tolerance_ranges::{scan_tol_range, tof_tol_range};
use log::{info, trace, warn};
use timsrust::{
    readers::{FrameReader, FrameReaderError, MetadataReaderError},
    TimsRustError,
};

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

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum UnsortedState {}
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SortedState {}

pub trait SortingStateTrait {}
impl SortingStateTrait for UnsortedState {}
impl SortingStateTrait for SortedState {}

#[derive(Debug, Clone)]
pub struct ExpandedFrameSlice<S: SortingStateTrait> {
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
    _sorting_state: PhantomData<S>,
}

fn sort_by_tof2(
    tof_indices: Vec<u32>,
    scan_numbers: Vec<usize>,
    intensities: Vec<u32>,
) -> ((Vec<u32>, Vec<usize>), Vec<u32>) {
    let mut combined: Vec<((u32, usize), u32)> = tof_indices
        .iter()
        .zip(scan_numbers.iter())
        .zip(intensities.iter())
        .map(|((tof, scan), inten)| ((*tof, *scan), *inten))
        .collect();
    combined.sort_unstable_by(|a, b| a.0 .0.partial_cmp(&b.0 .0).unwrap());
    let ((tof_indices, scan_numbers), intensities) = combined.into_iter().unzip();
    ((tof_indices, scan_numbers), intensities)
}

// Example of how to use it with your specific types
fn sort_by_tof_macro(
    tof_indices: Vec<u32>,
    scan_numbers: Vec<usize>,
    intensities: Vec<u32>,
) -> (Vec<u32>, Vec<usize>, Vec<u32>) {
    sort_vecs_by_first!(tof_indices, scan_numbers, intensities)
}

impl<T: SortingStateTrait> ExpandedFrameSlice<T> {
    pub fn sort_by_tof(self) -> ExpandedFrameSlice<SortedState> {
        // let mut indices = argsort_by(&self.tof_indices, |x| *x);
        // sort_by_indices_multi!(
        //     &mut indices,
        //     &mut self.tof_indices,
        //     &mut self.scan_numbers,
        //     &mut self.intensities
        // );

        let (tof_indices, scan_numbers, intensities) =
            sort_by_tof_macro(self.tof_indices, self.scan_numbers, self.intensities);

        ExpandedFrameSlice {
            tof_indices,
            scan_numbers,
            intensities,
            frame_index: self.frame_index,
            rt: self.rt,
            acquisition_type: self.acquisition_type,
            ms_level: self.ms_level,
            quadrupole_settings: self.quadrupole_settings,
            intensity_correction_factor: self.intensity_correction_factor,
            window_group: self.window_group,
            window_subindex: self.window_subindex,
            _sorting_state: PhantomData,
        }
    }

    pub fn len(&self) -> usize {
        self.tof_indices.len()
    }

    pub fn is_empty(&self) -> bool {
        self.tof_indices.is_empty()
    }
}

impl ExpandedFrameSlice<SortedState> {
    pub fn query_peaks<F>(
        &self,
        tof_range: (u32, u32),
        scan_range: Option<(usize, usize)>,
        f: &mut F,
    ) where
        F: FnMut(PeakInQuad),
    {
        let peak_range = {
            let peak_ind_start = self.tof_indices.partition_point(|x| x < &tof_range.0);
            let peak_ind_end = self.tof_indices.partition_point(|x| x <= &tof_range.1);
            peak_ind_start..peak_ind_end
        };

        match scan_range {
            Some((min_scan, max_scan)) => {
                assert!(min_scan <= max_scan);
            }
            None => {}
        };

        for peak_ind in peak_range {
            let scan_index = self.scan_numbers[peak_ind];
            match scan_range {
                Some((min_scan, max_scan)) => {
                    if scan_index < min_scan || scan_index > max_scan {
                        continue;
                    }
                }
                None => {}
            };

            let intensity = self.intensities[peak_ind];
            let tof_index = self.tof_indices[peak_ind];
            let retention_time = self.rt;
            f(PeakInQuad {
                scan_index,
                intensity,
                tof_index,
                retention_time: retention_time as f32,
            });
        }
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

fn expand_unfragmented_frame(frame: Frame) -> ExpandedFrameSlice<SortedState> {
    let scan_numbers = explode_vec(&frame.scan_offsets);
    let intensities = frame.intensities;
    let tof_indices = frame.tof_indices;
    let curr_slice = ExpandedFrameSlice {
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
        _sorting_state: PhantomData::<UnsortedState>,
    };

    curr_slice.sort_by_tof()
}

fn expand_fragmented_frame(
    frame: Frame,
    quads: Vec<SingleQuadrupoleSetting>,
) -> Vec<ExpandedFrameSlice<SortedState>> {
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

        let curr_slice = ExpandedFrameSlice {
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
            _sorting_state: PhantomData::<UnsortedState>,
        };
        // Q: is this sorting twice? since sort before creating the expanded frame slices.
        // TODO: Use the state type pattern to make sure only one sort happens.
        out.push(curr_slice.sort_by_tof());
    }
    out
}

pub fn expand_and_split_frame(frame: Frame) -> Vec<ExpandedFrameSlice<SortedState>> {
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

pub fn par_expand_and_arrange_frames(
    frame_iter: impl ParallelIterator<Item = Frame>,
) -> HashMap<Option<SingleQuadrupoleSetting>, Vec<ExpandedFrameSlice<SortedState>>> {
    let start = Instant::now();
    let split = frame_iter.flat_map(|frame| expand_and_split_frame(frame));
    // let mut out = HashMap::new();
    // for es in split {
    //     out.entry(es.quadrupole_settings)
    //         .or_insert(Vec::new())
    //         .push(es);
    // }
    //
    // Attempted implementation so the folding occurs in parrallel.
    // Inspired by/Copied from: https://stackoverflow.com/a/70097253/4295016
    let mut out = split
        .fold(HashMap::new, |mut acc, x| {
            acc.entry(x.quadrupole_settings)
                .or_insert(Vec::new())
                .push(x);
            acc
        })
        .reduce_with(|mut acc1, acc2| {
            for (qs, frameslices) in acc2.into_iter() {
                acc1.entry(qs).or_insert(Vec::new()).extend(frameslices);
            }
            acc1
        })
        .expect("At least one frame is iterated over");

    // Finally sort each of them internally by retention time.
    for (_, es) in out.iter_mut() {
        es.par_sort_unstable_by(|a, b| a.rt.partial_cmp(&b.rt).unwrap());
    }
    let end = start.elapsed();
    info!("Expanding and arranging frames took {:#?}", end);
    out
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum FrameProcessingConfig {
    Centroided {
        ims_tol_pct: f64,
        mz_tol_ppm: f64,
        window_width: usize,
        max_ms1_peaks: usize,
        max_ms2_peaks: usize,
        ims_converter: Option<Scan2ImConverter>,
        mz_converter: Option<Tof2MzConverter>,
    },
    NotCentroided,
}

impl FrameProcessingConfig {
    pub fn with_converters(
        self,
        ims_converter: Scan2ImConverter,
        mz_converter: Tof2MzConverter,
    ) -> Self {
        match self {
            FrameProcessingConfig::Centroided {
                ims_tol_pct,
                mz_tol_ppm,
                window_width,
                max_ms1_peaks,
                max_ms2_peaks,
                ..
            } => FrameProcessingConfig::Centroided {
                ims_tol_pct,
                mz_tol_ppm,
                window_width,
                max_ms1_peaks,
                max_ms2_peaks,
                ims_converter: Some(ims_converter),
                mz_converter: Some(mz_converter),
            },
            FrameProcessingConfig::NotCentroided => FrameProcessingConfig::NotCentroided,
        }
    }

    pub fn default_centroided() -> Self {
        FrameProcessingConfig::Centroided {
            ims_tol_pct: 1.5,
            mz_tol_ppm: 15.0,
            window_width: 3,
            max_ms1_peaks: 100_000,
            max_ms2_peaks: 10_000,
            ims_converter: Default::default(),
            mz_converter: Default::default(),
        }
    }

    pub fn default_not_centroided() -> Self {
        FrameProcessingConfig::NotCentroided
    }
}

#[derive(Debug)]
pub enum DataReadingError {
    CentroidingError(FrameProcessingConfig),
    UnsupportedDataError(String),
    TimsRustError(TimsRustError), // Why doesnt timsrust error derive clone?
}

impl From<TimsRustError> for DataReadingError {
    fn from(e: TimsRustError) -> Self {
        DataReadingError::TimsRustError(e)
    }
}

impl From<MetadataReaderError> for DataReadingError {
    fn from(e: MetadataReaderError) -> Self {
        DataReadingError::TimsRustError(TimsRustError::MetadataReaderError(e))
    }
}

impl From<FrameReaderError> for DataReadingError {
    fn from(e: FrameReaderError) -> Self {
        DataReadingError::TimsRustError(TimsRustError::FrameReaderError(e))
    }
}

fn warn_and_skip_badframes(
    frame_iter: impl ParallelIterator<Item = Result<Frame, FrameReaderError>>,
) -> impl ParallelIterator<Item = Frame> {
    frame_iter.filter_map(|x| {
        // Log the info of the frame that broke ...
        match x {
            Ok(frame) => Some(frame),
            Err(e) => {
                warn!("Failed to read frame {:?}", e);
                None
            }
        }
    })
}

pub fn par_read_and_expand_frames(
    frame_reader: &FrameReader,
    centroiding_config: FrameProcessingConfig,
) -> Result<
    HashMap<Option<SingleQuadrupoleSetting>, Vec<ExpandedFrameSlice<SortedState>>>,
    DataReadingError,
> {
    let dia_windows = match frame_reader.get_dia_windows() {
        Some(dia_windows) => dia_windows,
        None => {
            return Err(DataReadingError::UnsupportedDataError(
                "No dia windows found".to_string(),
            ))
        }
    };

    let mut all_expanded_frames = HashMap::new();
    for dia_window in dia_windows.into_iter() {
        info!("Processing dia window: {:?}", dia_window);
        let curr_iter = frame_reader.parallel_filter(|x| x.quadrupole_settings == dia_window);
        let curr_iter = warn_and_skip_badframes(curr_iter);

        let expanded_frames = match centroiding_config {
            FrameProcessingConfig::Centroided {
                ims_tol_pct,
                mz_tol_ppm,
                window_width,
                max_ms1_peaks: _max_ms1_peaks,
                max_ms2_peaks,
                ims_converter,
                mz_converter,
            } => {
                let expanded_frames = par_expand_and_centroid_frames(
                    curr_iter,
                    ims_tol_pct,
                    mz_tol_ppm,
                    window_width,
                    max_ms2_peaks,
                    &ims_converter.unwrap(),
                    &mz_converter.unwrap(),
                );
                expanded_frames
            }
            FrameProcessingConfig::NotCentroided => {
                let expanded_frames = par_expand_and_arrange_frames(curr_iter);
                expanded_frames
            }
        };

        all_expanded_frames.extend(expanded_frames);
    }

    info!("Processing MS1 frames");
    let ms1_iter = frame_reader.parallel_filter(|x| x.ms_level == MSLevel::MS1);
    let ms1_iter = warn_and_skip_badframes(ms1_iter);
    let expanded_ms1_frames = match centroiding_config {
        FrameProcessingConfig::Centroided {
            ims_tol_pct,
            mz_tol_ppm,
            window_width,
            max_ms1_peaks,
            max_ms2_peaks: _max_ms2_peaks,
            ims_converter,
            mz_converter,
        } => {
            let expanded_frames = par_expand_and_centroid_frames(
                ms1_iter,
                ims_tol_pct,
                mz_tol_ppm,
                window_width,
                max_ms1_peaks,
                &ims_converter.unwrap(),
                &mz_converter.unwrap(),
            );
            expanded_frames
        }
        FrameProcessingConfig::NotCentroided => {
            let expanded_frames = par_expand_and_arrange_frames(ms1_iter);
            expanded_frames
        }
    };
    all_expanded_frames.extend(expanded_ms1_frames);
    info!("Done reading and expanding frames");

    Ok(all_expanded_frames)
}

pub fn par_expand_and_centroid_frames(
    frames: impl ParallelIterator<Item = Frame>,
    ims_tol_pct: f64,
    mz_tol_ppm: f64,
    window_width: usize,
    max_peaks: usize,
    ims_converter: &Scan2ImConverter,
    mz_converter: &Tof2MzConverter,
) -> HashMap<Option<SingleQuadrupoleSetting>, Vec<ExpandedFrameSlice<SortedState>>> {
    let split_frames = par_expand_and_arrange_frames(frames);

    let out: HashMap<Option<SingleQuadrupoleSetting>, Vec<ExpandedFrameSlice<SortedState>>> =
        split_frames
            .into_iter()
            .map(|(qs, frameslices)| {
                // NOTE: Since the the centroiding runs in paralel over the windows, its ok if this
                // outer loop is done in series.
                let start_peaks: usize = frameslices.iter().map(|x| x.len()).sum();
                let centroided = par_lazy_centroid_frameslices(
                    &frameslices,
                    window_width,
                    ims_tol_pct,
                    mz_tol_ppm,
                    max_peaks,
                    ims_converter,
                    mz_converter,
                );
                let end_peaks: usize = centroided.iter().map(|x| x.len()).sum();
                trace!(
                    "Peak counts for quad {:?}: raw={}/centroid={}",
                    qs,
                    start_peaks,
                    end_peaks
                );
                (qs, centroided)
            })
            .collect();

    out
}

fn centroid_frameslice_window2(
    frameslices: &[ExpandedFrameSlice<SortedState>],
    ims_tol_pct: f64,
    mz_tol_ppm: f64,
    max_peaks: usize,
    ims_converter: &Scan2ImConverter,
    mz_converter: &Tof2MzConverter,
) -> ExpandedFrameSlice<SortedState> {
    assert!(frameslices.len() > 1, "Expected at least 2 frameslices");
    let reference_index = frameslices.len() / 2;

    let peak_refs: Vec<PeakArrayRefs> = frameslices
        .iter()
        .map(|x| PeakArrayRefs::new(&x.tof_indices, &x.scan_numbers, &x.intensities))
        .collect();

    let ((tof_array, intensity_array), ims_array) = lazy_centroid_weighted_frame(
        &peak_refs,
        reference_index,
        max_peaks,
        |tof| tof_tol_range(tof, mz_tol_ppm, mz_converter),
        |scan| scan_tol_range(scan, ims_tol_pct, ims_converter),
    );

    ExpandedFrameSlice {
        tof_indices: tof_array,
        scan_numbers: ims_array,
        intensities: intensity_array,
        frame_index: frameslices[reference_index].frame_index,
        rt: frameslices[reference_index].rt,
        acquisition_type: frameslices[reference_index].acquisition_type,
        ms_level: frameslices[reference_index].ms_level,
        quadrupole_settings: frameslices[reference_index].quadrupole_settings,
        intensity_correction_factor: frameslices[reference_index].intensity_correction_factor,
        window_group: frameslices[reference_index].window_group,
        window_subindex: frameslices[reference_index].window_subindex,
        _sorting_state: PhantomData::<SortedState>,
    }
}

pub fn par_lazy_centroid_frameslices(
    frameslices: &[ExpandedFrameSlice<SortedState>],
    window_width: usize,
    ims_tol_pct: f64,
    mz_tol_ppm: f64,
    max_peaks: usize,
    ims_converter: &Scan2ImConverter,
    mz_converter: &Tof2MzConverter,
) -> Vec<ExpandedFrameSlice<SortedState>> {
    assert!(
        frameslices
            .iter()
            .tuple_windows()
            .map(|(a, b)| { a.rt < b.rt && a.quadrupole_settings == b.quadrupole_settings })
            .all(|x| x),
        "All frames should be sorted by rt and have the same quad settings"
    );

    assert!(frameslices.len() > window_width);

    let local_lambda = |fss: &[ExpandedFrameSlice<SortedState>]| {
        centroid_frameslice_window2(
            fss,
            ims_tol_pct,
            mz_tol_ppm,
            max_peaks,
            ims_converter,
            mz_converter,
        )
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
