#[derive(Debug, Clone, Copy)]
pub struct RawPeak {
    pub scan_index: usize,
    pub tof_index: u32,
    pub intensity: u32,
}
