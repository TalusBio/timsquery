pub mod base;
#[allow(clippy::module_inception)]
pub mod multi_chromatogram_agg;

pub use multi_chromatogram_agg::{
    MultiCMGStatsAgg,
    MultiCMGStatsFactory,
};
