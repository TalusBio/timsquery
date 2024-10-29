use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::hash::Hash;

/// A struct that represents an elution group.
///
/// The elution group is a single precursor ion that is framented.
/// The fragments m/z values are stored in a hashmap where the key is
/// the generic type `T` and the value is the fragment m/z.
///
#[derive(Debug, Serialize, Deserialize, Clone)]
pub struct ElutionGroup<T: Clone + Eq + Serialize + Hash + Send + Sync> {
    pub id: u64,
    pub mobility: f32,
    pub rt_seconds: f32,
    pub precursor_mz: f64,
    pub precursor_charge: u8,
    pub fragment_mzs: HashMap<T, f64>,
}
