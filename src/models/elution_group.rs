use crate::HasIntegerID;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::hash::Hash;

#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct ElutionGroup<T: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub id: u64,
    pub mobility: f32,
    pub rt_seconds: f32,
    pub precursor_mz: f64,
    pub precursor_charge: u8,
    pub fragment_mzs: HashMap<T, f64>,
}

impl<T: Clone + Eq + Serialize + Hash + Send + Sync> HasIntegerID for ElutionGroup<T> {
    fn get_id(&self) -> u64 {
        self.id
    }
}
