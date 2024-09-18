use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize)]
pub struct ElutionGroup {
    pub id: u64,
    pub mobility: f32,
    pub rt_seconds: f32,
    pub precursor_mz: f64,
    pub precursor_charge: u8,
    pub fragment_mzs: Option<Vec<f64>>,
    pub fragment_charges: Option<Vec<u8>>,
}
