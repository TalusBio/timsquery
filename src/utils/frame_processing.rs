use log::warn;
use std::cmp::Ordering;
use std::ops::RangeInclusive;
use timsrust::converters::{Scan2ImConverter, Tof2MzConverter};

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

        let slice_width = ss_end - ss_start;
        let local_num_touched = touched[ss_start..ss_end].iter().filter(|&&x| x).count();
        let local_num_untouched = slice_width - local_num_touched;

        if local_num_touched == slice_width {
            continue;
        }

        let mut curr_intensity = 0.0;
        let mut curr_weighted_mz = 0.0;

        for i in ss_start..ss_end {
            if !touched[i] && intensity_array[i] > 0.0 {
                curr_intensity += intensity_array[i];
                curr_weighted_mz += mz_array[i] * intensity_array[i];
            }
        }

        if curr_intensity > 0.0 {
            curr_weighted_mz /= curr_intensity;

            agg_intensity[idx] = curr_intensity;
            agg_mz[idx] = curr_weighted_mz;

            touched[ss_start..ss_end].iter_mut().for_each(|x| *x = true);
            global_num_touched += local_num_untouched;
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

pub fn lazy_centroid_weighted_frame(
    tof_array: &[u32],
    ims_array: &[usize],
    weight_array: &[u32],
    intensity_array: &[u32],
    tof_tol_range_fn: impl Fn(&u32) -> RangeInclusive<u32>,
    ims_tol_range_fn: impl Fn(&usize) -> RangeInclusive<usize>,
) -> CentroidedVecs {
    let arr_len = tof_array.len();
    const MIN_WEIGHT_PRESERVE: u64 = 50;

    // TODO make the asserts optional at compile time.
    assert!(
        tof_array.windows(2).all(|x| x[0] <= x[1]),
        "Expected tof array to be sorted"
    );
    assert_eq!(
        arr_len,
        intensity_array.len(),
        "Expected tof and intensity arrays to be the same length"
    );
    assert_eq!(
        arr_len,
        weight_array.len(),
        "Expected tof and weight arrays to be the same length"
    );
    assert_eq!(
        arr_len,
        ims_array.len(),
        "Expected tof and ims arrays to be the same length"
    );

    let mut touched = vec![false; arr_len];
    let mut global_num_touched = 0;

    // We will be iterating in decreasing order of intensity
    let mut order: Vec<usize> = (0..arr_len).collect();
    order.sort_unstable_by(|&a, &b| {
        weight_array[b]
            .partial_cmp(&weight_array[a])
            .unwrap_or(Ordering::Equal)
    });

    let mut agg_tof = vec![0; arr_len];
    let mut agg_intensity = vec![0; arr_len];
    // We will not be returning the weights.
    // let mut agg_weight = vec![0; arr_len];
    let mut agg_ims = vec![0; arr_len];

    for idx in order {
        if touched[idx] {
            continue;
        }

        let tof = tof_array[idx];
        let ims = ims_array[idx];
        let tof_range = tof_tol_range_fn(&tof);
        let ims_range = ims_tol_range_fn(&ims);

        let ss_start = tof_array.partition_point(|x| x < tof_range.start());
        let ss_end = tof_array.partition_point(|x| x <= tof_range.end());

        let mut curr_intensity = 0u64;
        let mut curr_weight = 0u64;
        let mut curr_agg_tof = 0u64;
        let mut curr_agg_ims = 0u64;

        for i in ss_start..ss_end {
            if !touched[i] && intensity_array[i] > 0 && ims_range.contains(&ims_array[i]) {
                let local_weight = weight_array[i] as u64;
                curr_intensity += intensity_array[i] as u64;
                curr_agg_tof += tof_array[i] as u64 * local_weight;
                curr_agg_ims += ims_array[i] as u64 * local_weight;
                curr_weight += local_weight;
            }
        }

        if curr_intensity > 0 && curr_weight > 0 {
            agg_intensity[idx] = u32::try_from(curr_intensity).expect("Expected to fit in u32");
            agg_tof[idx] = (curr_agg_tof / curr_weight) as u32;
            agg_ims[idx] = (curr_agg_ims / curr_weight) as usize;

            for i in ss_start..ss_end {
                touched[i] = true;
                global_num_touched += 1;
            }
        }

        if global_num_touched == arr_len {
            break;
        }
    }

    // Drop the zeros and sort by mz (tof)
    let mut res: Vec<((u32, u32), usize)> = agg_tof
        .into_iter()
        .zip(agg_intensity.into_iter())
        .zip(agg_ims.into_iter())
        .filter(|&((_, intensity), ims)| (intensity > 0) & (ims > 0))
        .collect();

    res.sort_unstable_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(Ordering::Equal));
    let output_len = res.len();
    if arr_len > 500 && (output_len == arr_len) {
        warn!("Output length is the same as input length, this is probably a bug");
    }

    res.into_iter().unzip()
}
