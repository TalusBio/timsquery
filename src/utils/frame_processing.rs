use crate::sort_vecs_by_first;
use std::ops::RangeInclusive;
use tracing::{error, info, warn};

pub type TofIntensityVecs = (Vec<u32>, Vec<u32>);
pub type CentroidedVecs = (TofIntensityVecs, Vec<usize>);

fn sort_n_check(
    agg_intensity: Vec<u32>,
    agg_tof: Vec<u32>,
    agg_ims: Vec<usize>,
) -> ((Vec<u32>, Vec<u32>), Vec<usize>) {
    let (tof_array, ims_array, intensity_array) =
        sort_vecs_by_first!(agg_tof, agg_ims, agg_intensity);

    if let Some(x) = intensity_array.last() {
        assert!(
            *x > 0,
            "Expected all intensities to be positive non-zero, got {:?}",
            intensity_array
        );
        let max_tof = tof_array.iter().max().expect("At least one element");
        assert!(*max_tof < (u32::MAX - 1));
    }
    ((tof_array, intensity_array), ims_array)
}

pub struct PeakArrayRefs<'a> {
    pub tof_array: &'a [u32],
    pub ims_array: &'a [usize],
    pub intensity_array: &'a [u32],
}

impl<'a> PeakArrayRefs<'a> {
    pub fn new(
        tof_array: &'a [u32],
        ims_array: &'a [usize],
        intensity_array: &'a [u32],
    ) -> PeakArrayRefs<'a> {
        // TODO make the asserts optional at compile time.
        assert!(tof_array.len() == intensity_array.len());
        assert!(tof_array.len() == ims_array.len());
        assert!(
            tof_array.windows(2).all(|x| x[0] <= x[1]),
            "Expected tof array to be sorted"
        );
        Self {
            tof_array,
            ims_array,
            intensity_array,
        }
    }

    fn len(&self) -> usize {
        self.tof_array.len()
    }
}

// TODO: Refactor this function
pub fn lazy_centroid_weighted_frame<'a>(
    peak_refs: &'a [PeakArrayRefs<'a>],
    reference_index: usize,
    max_peaks: usize,
    tof_tol_range_fn: impl Fn(u32) -> RangeInclusive<u32>,
    ims_tol_range_fn: impl Fn(usize) -> RangeInclusive<usize>,
) -> CentroidedVecs {
    let slice_sizes: Vec<usize> = peak_refs.iter().map(|x| x.len()).collect();
    let tot_size: usize = slice_sizes.iter().sum();
    let ref_arrays = &peak_refs[reference_index];
    let num_arrays = peak_refs.len();
    assert!(num_arrays > 1);
    let arr_len = ref_arrays.len();
    if arr_len == 0 {
        error!("No peaks in reference array when centroiding");
        return ((Vec::new(), Vec::new()), Vec::new());
    }
    let initial_tot_intensity = ref_arrays
        .intensity_array
        .iter()
        .map(|x| *x as u64)
        .sum::<u64>();

    let mut touched = vec![false; tot_size];
    let mut global_num_touched = 0;
    let mut num_added = 0;

    // We will be iterating in decreasing order of intensity
    // Pre-calculate indices and intensities for sorting
    struct OrderItem {
        global_idx: usize,
        major_idx: usize,
        minor_idx: usize,
        intensity: u32,
    }

    let mut order: Vec<OrderItem> = Vec::with_capacity(tot_size);
    let mut offset = 0;
    for (major_idx, peak_ref) in peak_refs.iter().enumerate() {
        for (minor_idx, &intensity) in peak_ref.intensity_array.iter().enumerate() {
            order.push(OrderItem {
                global_idx: offset + minor_idx,
                major_idx,
                minor_idx,
                intensity,
            });
        }
        offset += peak_ref.len();
    }

    // Sort by intensity
    order.sort_unstable_by(|a, b| b.intensity.cmp(&a.intensity));
    assert!(order[0].intensity > order[tot_size - 1].intensity);

    let capacity = max_peaks.min(arr_len);
    let mut agg_tof = Vec::with_capacity(capacity);
    let mut agg_intensity = Vec::with_capacity(capacity);
    let mut agg_ims = Vec::with_capacity(capacity);

    for item in order {
        let major_idx = item.major_idx;
        let minor_idx = item.minor_idx;
        let this_intensity = item.intensity as u64;

        if touched[item.global_idx] {
            continue;
        }

        let tof = peak_refs[major_idx].tof_array[minor_idx];
        let ims = peak_refs[major_idx].ims_array[minor_idx];
        let tof_range = tof_tol_range_fn(tof);
        let ims_range = ims_tol_range_fn(ims);

        let mut curr_intensity = 0u64;
        let mut curr_weight = 0u64;
        // let mut curr_agg_tof = 0u64;
        // let mut curr_agg_ims = 0u64;

        // This is added to the index within the loop
        // To get the touched status of the peaks.
        let mut local_offset_touched = 0;
        for (ii, local_peak_refs) in peak_refs.iter().enumerate() {
            let ss_start = local_peak_refs
                .tof_array
                .partition_point(|x| x < tof_range.start());
            let ss_end = local_peak_refs
                .tof_array
                .partition_point(|x| x <= tof_range.end());
            for i in ss_start..ss_end {
                let ti = local_offset_touched + i;
                if !touched[ti] && ims_range.contains(&local_peak_refs.ims_array[i]) {
                    // Peaks are always weighted but not always intense!
                    let local_intensity = local_peak_refs.intensity_array[i] as u64;
                    if ii == reference_index {
                        curr_intensity += local_intensity;
                        global_num_touched += 1;
                    }
                    // curr_agg_tof += local_peak_refs.tof_array[i] as u64 * local_intensity;
                    // curr_agg_ims += local_peak_refs.ims_array[i] as u64 * local_intensity;
                    curr_weight += local_intensity;
                    touched[ti] = true;
                }
            }
            local_offset_touched += local_peak_refs.len();
        }

        // This means that at least 2 peaks need to be aggregated.
        if curr_weight > this_intensity && curr_intensity >= this_intensity {
            agg_intensity.push(u32::try_from(curr_intensity).expect("Expected to fit in u32"));
            // let calc_tof = (curr_agg_tof / curr_weight) as u32;
            // let calc_ims = (curr_agg_ims / curr_weight) as usize;
            let calc_tof = tof;
            let calc_ims = ims;
            debug_assert!(tof_range.contains(&calc_tof));
            debug_assert!(ims_range.contains(&calc_ims));
            agg_tof.push(calc_tof);
            agg_ims.push(calc_ims);
            num_added += 1;
            if num_added == max_peaks {
                break;
            }
        }

        if global_num_touched == arr_len {
            break;
        }
    }

    let out = sort_n_check(agg_intensity, agg_tof, agg_ims);

    // TODO:Make everything below this a separate function and accumulate it.
    let tot_final_intensity = out.0 .1.iter().map(|x| *x as u64).sum::<u64>();
    let inten_ratio = tot_final_intensity as f64 / initial_tot_intensity as f64;
    assert!(initial_tot_intensity >= tot_final_intensity);

    let output_len = out.0 .0.len();
    let compression_ratio = output_len as f64 / arr_len as f64;
    assert!(num_added == output_len);

    // 80% of the intensity being preserved sounds like a good cutoff.
    if output_len == max_peaks && inten_ratio < 0.80 {
        info!(
            "Frame trimmed to max peaks ({}), preserved intensity {}/{} intensity ratio: {} compression_ratio: {} if this is not acceptable consider increasing the parameter.",
            max_peaks, tot_final_intensity, initial_tot_intensity, inten_ratio, compression_ratio,
        );
    }

    if arr_len > 5000 && (output_len == arr_len) {
        warn!("Output length is the same as input length, this is probably a bug");
        warn!("Intensity ratio {:?}", inten_ratio);
        warn!("initial_tot_intensity: {:?}", initial_tot_intensity);
        warn!("tot_final_intensity: {:?}", tot_final_intensity);
        warn!(
            "First tof {} -> Range {:?}",
            out.0 .0[0],
            tof_tol_range_fn(out.0 .0[0])
        );
        panic!();
        // warn!("agg_intensity: {:?}", out.0 .0);
    }

    out
}

#[cfg(test)]
mod tests {
    use super::*;

    // Helper function to create a simple tolerance range for testing
    fn test_tof_tolerance(tof: u32) -> RangeInclusive<u32> {
        let tolerance = 2;
        (tof.saturating_sub(tolerance))..=tof.saturating_add(tolerance)
    }

    fn test_ims_tolerance(ims: usize) -> RangeInclusive<usize> {
        let tolerance = 1;
        (ims.saturating_sub(tolerance))..=ims.saturating_add(tolerance)
    }

    // Helper function to create PeakArrayRefs
    fn create_peak_refs<'a>(
        tof: &'a [u32],
        ims: &'a [usize],
        intensity: &'a [u32],
    ) -> PeakArrayRefs<'a> {
        PeakArrayRefs::new(tof, ims, intensity)
    }

    #[test]
    fn test_basic_centroid() {
        // Test case with two simple peaks that should be merged
        let peak_refs = vec![
            create_peak_refs(
                &[100, 101, 200, 201],   // tof
                &[1, 1, 2, 2],           // ims
                &[1000, 2000, 500, 600], // intensity
            ),
            create_peak_refs(
                &[101, 201], // tof
                &[1, 2],     // ims
                &[900, 600], // intensity
            ),
        ];

        let result = lazy_centroid_weighted_frame(
            &peak_refs,
            0,  // reference_index
            10, // max_peaks
            test_tof_tolerance,
            test_ims_tolerance,
        );

        let ((tof_array, intensity_array), ims_array) = result;

        assert_eq!(
            tof_array.len(),
            2,
            "Expected 2 centroids, got {:?}, {:?}, {:?}",
            tof_array,
            intensity_array,
            ims_array
        );
        assert_eq!(intensity_array.len(), 2);
        assert_eq!(ims_array.len(), 2);

        // Check that the peaks were properly merged and weighted
        assert!(tof_array[0] >= 100 && tof_array[0] <= 101);
        assert_eq!(ims_array[0], 1);
        assert!(intensity_array[0] == 3000); // Should be sum of intensities
        assert!(intensity_array[1] == 1100); // Should be sum of intensities
    }

    #[test]
    fn test_max_peaks_limit() {
        // Test that the function respects the max_peaks parameter
        let peak_refs = vec![
            create_peak_refs(&[100, 200, 300], &[1, 2, 3], &[1000, 900, 800]),
            create_peak_refs(&[101, 201, 301], &[1, 2, 3], &[950, 850, 750]),
        ];

        let result = lazy_centroid_weighted_frame(
            &peak_refs,
            0,
            2, // max_peaks set to 2
            test_tof_tolerance,
            test_ims_tolerance,
        );

        let ((tof_array, _), _) = result;
        assert_eq!(
            tof_array.len(),
            2,
            "Should only return max_peaks number of peaks"
        );
    }

    #[test]
    fn test_empty_input() {
        // Test handling of empty arrays
        let peak_refs = vec![
            create_peak_refs(&[], &[], &[]),
            create_peak_refs(&[], &[], &[]),
        ];

        let result =
            lazy_centroid_weighted_frame(&peak_refs, 0, 10, test_tof_tolerance, test_ims_tolerance);

        let ((tof_array, intensity_array), ims_array) = result;
        assert_eq!(tof_array.len(), 0);
        assert_eq!(intensity_array.len(), 0);
        assert_eq!(ims_array.len(), 0);
    }

    #[test]
    fn test_intensity_weighted_centroid() {
        // Test that centroids are properly weighted by intensity
        let peak_refs = vec![
            create_peak_refs(
                &[100, 200],
                &[1, 2],
                &[1000, 100], // High intensity for first peak
            ),
            create_peak_refs(
                &[102, 202],
                &[1, 2],
                &[100, 1000], // High intensity for second peak
            ),
        ];

        let result =
            lazy_centroid_weighted_frame(&peak_refs, 0, 10, test_tof_tolerance, test_ims_tolerance);

        let ((tof_array, _), ims_array) = result;

        // First centroid should be closer to 100 due to higher intensity
        assert!(tof_array[0] - 100 < 102 - tof_array[0]);
        assert_eq!(ims_array[0], 1);
    }

    #[test]
    #[should_panic]
    fn test_invalid_reference_index() {
        let peak_refs = vec![
            create_peak_refs(&[100], &[1], &[1000]),
            create_peak_refs(&[101], &[1], &[900]),
        ];

        // This should panic due to invalid reference index
        lazy_centroid_weighted_frame(
            &peak_refs,
            2, // Invalid index
            10,
            test_tof_tolerance,
            test_ims_tolerance,
        );
    }

    // Pretty decent test suggested by claude but im not sure if I really
    // need to suuport it ... since scan indices are usually less than 2_000
    // #[test]
    // fn test_large_input_values() {
    //     // Test handling of large values close to u32::MAX
    //     let peak_refs = vec![
    //         create_peak_refs(
    //             &[u32::MAX - 10, u32::MAX - 5],
    //             &[usize::MAX - 10, usize::MAX - 5],
    //             &[1000, 900],
    //         ),
    //         create_peak_refs(
    //             &[u32::MAX - 9, u32::MAX - 4],
    //             &[usize::MAX - 9, usize::MAX - 4],
    //             &[950, 850],
    //         ),
    //     ];

    //     let result =
    //         lazy_centroid_weighted_frame(&peak_refs, 0, 10, test_tof_tolerance, test_ims_tolerance);

    //     let ((tof_array, _), _) = result;
    //     assert!(!tof_array.is_empty(), "Should handle large values properly");
    // }
}
