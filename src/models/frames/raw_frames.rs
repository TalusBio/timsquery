use timsrust::ms_data::{Frame, QuadrupoleSettings};

use super::raw_peak::RawPeak;

pub fn scans_matching_quad(
    quad_settings: &QuadrupoleSettings,
    quad_range: (f64, f64),
) -> Option<(usize, usize)> {
    let mut min_start = usize::MAX;
    let mut max_end = 0;

    for i in 0..quad_settings.scan_ends.len() {
        let half_width = quad_settings.isolation_width[i] / 2.0;
        let quad_start = quad_settings.isolation_mz[i] - half_width;
        let quad_end = quad_settings.isolation_mz[i] + half_width;

        if quad_start <= quad_range.1 && quad_end >= quad_range.0 {
            let start = quad_settings.scan_starts[i];
            let end = quad_settings.scan_ends[i];
            min_start = min_start.min(start);
            max_end = max_end.max(end);
        }
    }

    if min_start == usize::MAX {
        None
    } else {
        Some((min_start, max_end))
    }
}

pub fn frame_elems_matching<'a>(
    frame: &'a Frame,
    tof_range: (u32, u32),
    scan_range: (usize, usize),
    quad_range: Option<(f64, f64)>,
) -> impl Iterator<Item = RawPeak> + 'a {
    let quad_range = quad_range
        .and_then(|quad_range| scans_matching_quad(&frame.quadrupole_settings, quad_range));

    let (scan_ind_start, scan_ind_end): (usize, usize) = match quad_range {
        Some((start, end)) => (start.max(scan_range.0), end.min(scan_range.1)),
        None => (scan_range.0, scan_range.1),
    };
    (scan_ind_start..scan_ind_end)
        .map(move |scan_index| {
            let scan_is = frame.scan_offsets[scan_index];
            let scan_ie = frame.scan_offsets[scan_index + 1];

            let outs: Vec<RawPeak> = (scan_is..scan_ie)
                .map(move |peak_index| {
                    let tof_ind = frame.tof_indices[peak_index];
                    let intensity = frame.intensities[peak_index];

                    (tof_ind, intensity, scan_index)
                })
                .filter(|(tof_ind, _, _)| *tof_ind >= tof_range.0 && *tof_ind < tof_range.1)
                .map(|(tof_ind, intensity, scan_index)| RawPeak {
                    scan_index,
                    tof_index: tof_ind,
                    intensity,
                })
                .collect();

            outs
        })
        .flatten()
}
