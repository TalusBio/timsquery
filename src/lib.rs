// Re-export main structures
pub use crate::models::elution_group::ElutionGroup;
pub use crate::queriable_tims_data::queriable_tims_data::QueriableTimsData;

// Re-export traits
pub use crate::traits::aggregator::Aggregator;
pub use crate::traits::indexed_data::IndexedData;
pub use crate::traits::tolerance::{Tolerance, ToleranceAdapter};

// Re-export utility functions
pub use crate::utils::sorting::sort_multiple_by;

// Declare modules
pub mod models;
pub mod queriable_tims_data;
pub mod traits;
pub mod utils;

// Any library-wide code, documentation, or additional exports can go here
