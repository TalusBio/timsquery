// Re-export main structures
pub use crate::models::{
    elution_group::ElutionGroup,
    indices::transposed_quad_index::QuadSplittedTransposedIndex,
};

// Re-export traits
pub use crate::traits::{
    aggregator::Aggregator,
    queriable_data::QueriableData,
    tolerance::{
        Tolerance,
        ToleranceAdapter,
    },
};

// Declare modules
pub mod errors;
pub mod models;
pub mod queriable_tims_data;
pub mod traits;
pub mod utils;

// Re-export errors
pub use crate::errors::{
    DataProcessingError,
    DataReadingError,
    TimsqueryError,
};
