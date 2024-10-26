use timsrust::converters::{ConvertableDomain, Scan2ImConverter, Tof2MzConverter};

use std::ops::RangeInclusive;

pub fn ppm_tol_range(elem: f64, tol_ppm: f64) -> RangeInclusive<f64> {
    let utol = elem * (tol_ppm / 1e6);
    let left_e = elem - utol;
    let right_e = elem + utol;
    left_e..=right_e
}

pub fn pct_tol_range(elem: f64, tol_pct: f64) -> RangeInclusive<f64> {
    let utol = elem * (tol_pct / 100.0);
    let left_e = elem - utol;
    let right_e = elem + utol;
    left_e..=right_e
}

pub fn tof_tol_range(tof: u32, tol_ppm: f64, converter: &Tof2MzConverter) -> RangeInclusive<u32> {
    let mz = converter.convert(tof);
    let mz_range = ppm_tol_range(mz, tol_ppm);
    RangeInclusive::new(
        converter.invert(*mz_range.start()).round() as u32,
        converter.invert(*mz_range.end()).round() as u32,
    )
}

pub fn scan_tol_range(
    scan: usize,
    tol_pct: f64,
    converter: &Scan2ImConverter,
) -> RangeInclusive<usize> {
    let im = converter.convert(scan as f64);
    let im_range = pct_tol_range(im, tol_pct);
    let scan_min = converter.invert(*im_range.start()).round() as usize;
    let scan_max = converter.invert(*im_range.end()).round() as usize;
    // Note I need to do this here bc the conversion between scan numbers and ion
    // mobilities is not monotonically increasing. IN OTHER WORDS, lower scan numbers
    // are higher 1/k0.... But im not sure if they are ALWAYS inversely proportional.
    if scan_min > scan_max {
        scan_max..=scan_min
    } else {
        scan_min..=scan_max
    }
}
