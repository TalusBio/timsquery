use crate::models::elution_group::ElutionGroup;
use crate::models::queries::{FragmentGroupIndexQuery, PrecursorIndexQuery};
use crate::ToleranceAdapter;
use serde::Serialize;
use std::hash::Hash;
use timsrust::converters::ConvertableDomain;
use timsrust::Metadata;

#[derive(Debug, Default)]
pub struct FragmentIndexAdapter {
    pub metadata: Metadata,
}

impl From<Metadata> for FragmentIndexAdapter {
    fn from(metadata: Metadata) -> Self {
        Self { metadata }
    }
}

impl<FH: Copy + Clone + Serialize + Eq + Hash + Send + Sync>
    ToleranceAdapter<FragmentGroupIndexQuery<FH>, ElutionGroup<FH>> for FragmentIndexAdapter
{
    fn query_from_elution_group(
        &self,
        tol: &dyn crate::traits::tolerance::Tolerance,
        elution_group: &ElutionGroup<FH>,
    ) -> FragmentGroupIndexQuery<FH> {
        let rt_range = tol.rt_range(elution_group.rt_seconds);
        let mobility_range = tol.mobility_range(elution_group.mobility);
        let precursor_mz_range = tol.mz_range(elution_group.precursor_mz);
        let quad_range = tol.quad_range(elution_group.precursor_mz, elution_group.precursor_charge);

        let mz_index_range = (
            self.metadata.mz_converter.invert(precursor_mz_range.0) as u32,
            self.metadata.mz_converter.invert(precursor_mz_range.1) as u32,
        );
        let mobility_range = match mobility_range {
            Some(mobility_range) => mobility_range,
            None => (self.metadata.lower_im as f32, self.metadata.upper_im as f32),
        };
        let mobility_index_range = (
            self.metadata.im_converter.invert(mobility_range.0) as usize,
            self.metadata.im_converter.invert(mobility_range.1) as usize,
        );
        let rt_range = match rt_range {
            Some(rt_range) => rt_range,
            None => (self.metadata.lower_rt as f32, self.metadata.upper_rt as f32),
        };
        let frame_index_range = (
            self.metadata.rt_converter.invert(rt_range.0) as usize,
            self.metadata.rt_converter.invert(rt_range.1) as usize,
        );

        assert!(frame_index_range.0 <= frame_index_range.1);
        assert!(mz_index_range.0 <= mz_index_range.1);
        // Since mobilities get mixed up bc low scan ranges are high 1/k0, I
        // Just make sure they are sorted here.
        let mobility_index_range = (
            mobility_index_range.0.min(mobility_index_range.1),
            mobility_index_range.1.max(mobility_index_range.0),
        );

        let precursor_query = PrecursorIndexQuery {
            frame_index_range,
            rt_range_seconds: rt_range,
            mz_index_range,
            mobility_index_range,
            isolation_mz_range: quad_range,
        };

        let fqs = elution_group
            .fragment_mzs
            .iter()
            .map(|(k, v)| {
                let mz_range = tol.mz_range(*v);
                (
                    *k,
                    (
                        self.metadata.mz_converter.invert(mz_range.0) as u32,
                        self.metadata.mz_converter.invert(mz_range.1) as u32,
                    ),
                )
            })
            .collect();

        FragmentGroupIndexQuery {
            mz_index_ranges: fqs,
            precursor_query,
        }
    }
}
