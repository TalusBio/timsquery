use crate::sort_vecs_by_first;
use log::{info, warn};
use std::cmp::Ordering;
use std::ops::RangeInclusive;

pub fn squash_frame(
    mz_array: &[f32],
    intensity_array: &[f32],
    tol_ppm: f32,
) -> (Vec<f32>, Vec<f32>) {
    // Make sure the mz array is sorted
    assert!(mz_array.windows(2).all(|x| x[0] <= x[1]));

    let arr_len = mz_array.len();
    let mut touched = vec![false; arr_len];
    let mut global_num_touched = 0;

    let mut order: Vec<usize> = (0..arr_len).collect();
    order.sort_unstable_by(|&a, &b| {
        intensity_array[b]
            .partial_cmp(&intensity_array[a])
            .unwrap_or(Ordering::Equal)
    });

    let mut agg_mz = vec![0.0; arr_len];
    let mut agg_intensity = vec![0.0; arr_len];

    let utol = tol_ppm / 1e6;

    for &idx in &order {
        if touched[idx] {
            continue;
        }

        let mz = mz_array[idx];
        let da_tol = mz * utol;
        let left_e = mz - da_tol;
        let right_e = mz + da_tol;

        let ss_start = mz_array.partition_point(|&x| x < left_e);
        let ss_end = mz_array.partition_point(|&x| x <= right_e);

        let mut curr_intensity = 0.0;
        let mut curr_weighted_mz = 0.0;

        for i in ss_start..ss_end {
            if !touched[i] && intensity_array[i] > 0.0 {
                curr_intensity += intensity_array[i];
                curr_weighted_mz += mz_array[i] * intensity_array[i];
                touched[i] = true;
                global_num_touched += 1;
            }
        }

        if curr_intensity > 0.0 {
            curr_weighted_mz /= curr_intensity;

            agg_intensity[idx] = curr_intensity;
            agg_mz[idx] = curr_weighted_mz;

            touched[ss_start..ss_end].iter_mut().for_each(|x| *x = true);
        }

        if global_num_touched == arr_len {
            break;
        }
    }

    // Drop the zeros and sort
    let mut result: Vec<(f32, f32)> = agg_mz
        .into_iter()
        .zip(agg_intensity.into_iter())
        .filter(|&(mz, intensity)| mz > 0.0 && intensity > 0.0)
        .collect();

    result.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));

    result.into_iter().unzip()
}

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
        assert!(*x > 0);
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

/// Splits a global index into a major and minor that
/// matches elements in a collection.
///
/// As an example if there is a slice of slices, the major index
/// is the index of the slice and the minor index is the index
/// of the element in the slice.
///
/// Thus if we have      [[0,1,2], [3,4,5], [6,7,8]]
/// the major indice are   0,0,0    1,1,1    2,2,2
/// and the minor          0,1,2    0,1,2    0,1,2
///
///
fn index_split(index: usize, slice_sizes: &[usize]) -> (usize, usize) {
    let mut index = index;
    for (i, slice_size) in slice_sizes.iter().enumerate() {
        if index < *slice_size {
            return (i, index);
        }
        index -= *slice_size;
    }
    panic!("Index {} is out of bounds", index);
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
    let initial_tot_intensity = ref_arrays
        .intensity_array
        .iter()
        .map(|x| *x as u64)
        .sum::<u64>();

    const MIN_WEIGHT_PRESERVE: u64 = 50;

    let mut touched = vec![false; tot_size];
    let mut global_num_touched = 0;
    let mut num_added = 0;

    // We will be iterating in decreasing order of intensity
    let mut order: Vec<usize> = (0..tot_size).collect();
    order.sort_unstable_by(|&a, &b| {
        let (a_major, a_minor) = index_split(a, &slice_sizes);
        let (b_major, b_minor) = index_split(b, &slice_sizes);

        let a_intensity = peak_refs[a_major].intensity_array[a_minor];
        let b_intensity = peak_refs[b_major].intensity_array[b_minor];

        b_intensity
            .partial_cmp(&a_intensity)
            .unwrap_or(Ordering::Equal)
    });

    assert!({
        let (first_major, first_minor) = index_split(order[0], &slice_sizes);
        let (last_major, last_minor) = index_split(order[tot_size - 1], &slice_sizes);

        let first_intensity = peak_refs[first_major].intensity_array[first_minor];
        let last_intensity = peak_refs[last_major].intensity_array[last_minor];

        first_intensity >= last_intensity
    });

    // TODO explore if making vecs with capacity and appedning is better...
    let capacity = max_peaks.min(arr_len);
    let mut agg_tof = Vec::with_capacity(capacity);
    let mut agg_intensity = Vec::with_capacity(capacity);
    let mut agg_ims = Vec::with_capacity(capacity);

    for idx in order {
        let (major_idx, minor_idx) = index_split(idx, &slice_sizes);

        if touched[idx] {
            continue;
        }

        let tof = peak_refs[major_idx].tof_array[minor_idx];
        let ims = peak_refs[major_idx].ims_array[minor_idx];
        let this_intensity = peak_refs[major_idx].intensity_array[minor_idx] as u64;
        let tof_range = tof_tol_range_fn(tof);
        let ims_range = ims_tol_range_fn(ims);

        let mut curr_intensity = 0u64;
        let mut curr_weight = 0u64;
        let mut curr_agg_tof = 0u64;
        let mut curr_agg_ims = 0u64;

        let mut local_offset_touched = 0;
        for ii in 0..num_arrays {
            let local_peak_refs = &peak_refs[ii];
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
                    curr_agg_tof += local_peak_refs.tof_array[i] as u64 * local_intensity;
                    curr_agg_ims += local_peak_refs.ims_array[i] as u64 * local_intensity;
                    curr_weight += local_intensity;
                    touched[ti] = true;
                }
            }
            local_offset_touched += local_peak_refs.len();
        }

        if curr_intensity > this_intensity {
            agg_intensity.push(u32::try_from(curr_intensity).expect("Expected to fit in u32"));
            let calc_tof = (curr_agg_tof / curr_weight) as u32;
            let calc_ims = (curr_agg_ims / curr_weight) as usize;
            assert!(tof_range.contains(&calc_tof));
            assert!(ims_range.contains(&calc_ims));
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

    // Drop the zeros and sort by mz (tof)
    let out = sort_n_check(agg_intensity, agg_tof, agg_ims);
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
