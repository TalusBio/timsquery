use super::quad_index::{
    FrameRTTolerance, PeakInQuad, TransposedQuadIndex, TransposedQuadIndexBuilder,
};
use crate::models::elution_group::ElutionGroup;
use crate::models::frames::expanded_frame::{expand_and_split_frame, ExpandedFrameSlice};
use crate::models::frames::raw_peak::RawPeak;
use crate::models::frames::single_quad_settings::SingleQuadrupoleSetting;
use crate::models::queries::FragmentGroupIndexQuery;
use crate::models::queries::PrecursorIndexQuery;
use crate::traits::indexed_data::IndexedData;
use crate::utils::display::{glimpse_vec, GlimpseConfig};
use crate::ToleranceAdapter;
use log::{debug, info, trace};
use rayon::prelude::*;
use serde::Serialize;
use std::collections::HashMap;
use std::fmt::Debug;
use std::fmt::Display;
use std::hash::Hash;
use std::time::Instant;
use timsrust::converters::{
    ConvertableDomain, Frame2RtConverter, Scan2ImConverter, Tof2MzConverter,
};
use timsrust::readers::{FrameReader, MetadataReader};
use timsrust::Frame;
use timsrust::Metadata;
use timsrust::TimsRustError;

// TODO break this module apart ... its getting too big for my taste
// - JSP: 2024-11-19

pub struct QuadSplittedTransposedIndex {
    precursor_index: TransposedQuadIndex,
    fragment_indices: HashMap<SingleQuadrupoleSetting, TransposedQuadIndex>,
    flat_quad_settings: Vec<SingleQuadrupoleSetting>,
    rt_converter: Frame2RtConverter,
    pub mz_converter: Tof2MzConverter,
    pub im_converter: Scan2ImConverter,
    metadata: Metadata,
}

impl Debug for QuadSplittedTransposedIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "QuadSplittedTransposedIndex(num_quads: {})",
            self.flat_quad_settings.len()
        )
    }
}

impl QuadSplittedTransposedIndex {
    pub fn query_peaks(
        &self,
        tof_range: (u32, u32),
        precursor_mz_range: (f64, f64),
        scan_range: Option<(usize, usize)>,
        rt_range: Option<FrameRTTolerance>,
    ) -> impl Iterator<Item = PeakInQuad> + '_ {
        trace!(
            "QuadSplittedTransposedIndex::query_peaks tof_range: {:?}, scan_range: {:?}, rt_range: {:?}, precursor_mz_range: {:?}",
            tof_range,
            scan_range,
            rt_range,
            precursor_mz_range,
        );
        let matching_quads: Vec<SingleQuadrupoleSetting> = self
            .get_matching_quad_settings(precursor_mz_range, scan_range)
            .collect();
        trace!("matching_quads: {:?}", matching_quads);

        let rt_range = rt_range.map(|x| x.to_frame_index_range(&self.rt_converter));

        matching_quads.into_iter().flat_map(move |qs| {
            let tqi = self.fragment_indices.get(&qs).unwrap();
            tqi.query_peaks(tof_range, scan_range, rt_range)
        })
    }

    fn get_matching_quad_settings(
        &self,
        precursor_mz_range: (f64, f64),
        scan_range: Option<(usize, usize)>,
    ) -> impl Iterator<Item = SingleQuadrupoleSetting> + '_ {
        self.flat_quad_settings
            .iter()
            .filter(move |qs| {
                (qs.ranges.isolation_low <= precursor_mz_range.1)
                    && (qs.ranges.isolation_high >= precursor_mz_range.0)
            })
            .filter(move |qs| {
                // if let Some(scan_range) = scan_range {
                //     // This is done for sanity tbh ... sometimes they get flipped
                //     // bc the lowest scan is actually the highest 1/k0.
                // } else {
                //     true
                // }
                //
                match scan_range {
                    Some((min_scan, max_scan)) => {
                        assert!(qs.ranges.scan_start <= qs.ranges.scan_end);
                        assert!(min_scan <= max_scan);

                        // Above quad
                        // Quad                   [----------]
                        // Query                               [------]
                        let above_quad = qs.ranges.scan_end < min_scan;

                        // Below quad
                        // Quad                  [------]
                        // Query       [------]
                        let below_quad = qs.ranges.scan_start > max_scan;

                        if above_quad || below_quad {
                            // This quad is completely outside the scan range
                            false
                        } else {
                            true
                        }
                    }
                    None => true,
                }
            })
            .cloned()
    }

    fn queries_from_elution_elements_impl<FH: Clone + Eq + Serialize + Hash + Send + Sync>(
        &self,
        tol: &dyn crate::traits::tolerance::Tolerance,
        elution_group: &crate::models::elution_group::ElutionGroup<FH>,
    ) -> FragmentGroupIndexQuery<FH> {
        let rt_range = tol.rt_range(elution_group.rt_seconds);
        let mobility_range = tol.mobility_range(elution_group.mobility);
        let precursor_mz_range = tol.mz_range(elution_group.precursor_mz);
        let quad_range = tol.quad_range(elution_group.precursor_mz, elution_group.precursor_charge);

        let mz_index_range = (
            self.mz_converter.invert(precursor_mz_range.0) as u32,
            self.mz_converter.invert(precursor_mz_range.1) as u32,
        );
        let mobility_range = match mobility_range {
            Some(mobility_range) => mobility_range,
            None => (self.metadata.lower_im as f32, self.metadata.upper_im as f32),
        };
        let mobility_index_range = (
            self.im_converter.invert(mobility_range.0) as usize,
            self.im_converter.invert(mobility_range.1) as usize,
        );
        let rt_range = match rt_range {
            Some(rt_range) => rt_range,
            None => (self.metadata.lower_rt as f32, self.metadata.upper_rt as f32),
        };
        let frame_index_range = (
            self.rt_converter.invert(rt_range.0) as usize,
            self.rt_converter.invert(rt_range.1) as usize,
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
                    k.clone(),
                    (
                        self.mz_converter.invert(mz_range.0) as u32,
                        self.mz_converter.invert(mz_range.1) as u32,
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

impl Display for QuadSplittedTransposedIndex {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let mut disp_str = String::new();
        disp_str.push_str("QuadSplittedTransposedIndex\n");

        disp_str.push_str("rt_converter: ... not showing ...\n");
        disp_str.push_str(&format!("mz_converter: {:?}\n", self.mz_converter));
        disp_str.push_str(&format!("im_converter: {:?}\n", self.im_converter));
        disp_str.push_str("flat_quad_settings: \n");
        disp_str.push_str(&glimpse_vec(
            &self.flat_quad_settings,
            Some(GlimpseConfig {
                max_items: 10,
                padding: 4,
                new_line: true,
            }),
        ));

        disp_str.push_str("precursor_index: \n");
        disp_str.push_str(&format!(" -- {}\n", self.precursor_index));
        let mut num_shown = 0;
        for (qs, tqi) in self.fragment_indices.iter() {
            disp_str.push_str(&format!(" - {:?}: \n", qs));
            disp_str.push_str(&format!(" -- {}\n", tqi));
            num_shown += 1;
            if num_shown > 5 {
                disp_str.push_str(&format!(
                    " ........ len = {}\n",
                    self.fragment_indices.len()
                ));
                break;
            }
        }

        write!(f, "{}", disp_str)
    }
}

impl QuadSplittedTransposedIndex {
    pub fn from_path(path: &str) -> Result<Self, TimsRustError> {
        let st = Instant::now();
        info!("Building transposed quad index from path {}", path);
        let tmp = QuadSplittedTransposedIndexBuilder::from_path(path)?;
        let out = tmp.build();
        let elapsed = st.elapsed();
        info!("Transposed quad index built in {:#?}", elapsed);
        debug!("{}", out);
        Ok(out)
    }
}

#[derive(Debug, Clone, Default)]
pub struct QuadSplittedTransposedIndexBuilder {
    indices: HashMap<Option<SingleQuadrupoleSetting>, TransposedQuadIndexBuilder>,
    rt_converter: Option<Frame2RtConverter>,
    mz_converter: Option<Tof2MzConverter>,
    im_converter: Option<Scan2ImConverter>,
    metadata: Option<Metadata>,
    // TODO use during build to make sure we
    // have a the right number of peaks in the end.
    // ... this means I need to implement len for TransposedQuadIndexBuilder
    added_peaks: u64,
}

impl QuadSplittedTransposedIndexBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    fn add_frame(&mut self, frame: Frame) {
        let expanded_slices = expand_and_split_frame(frame);

        for es in expanded_slices.into_iter() {
            // Add key if it doesnt exist ...

            self.indices
                .entry(es.quadrupole_settings)
                .or_insert(TransposedQuadIndexBuilder::new(es.quadrupole_settings));

            self.added_peaks += es.len() as u64;
            self.add_frame_slice(es);
        }
    }

    fn add_frame_slice(&mut self, frame_slice: ExpandedFrameSlice) {
        self.indices
            .get_mut(&frame_slice.quadrupole_settings)
            .unwrap()
            .add_frame_slice(frame_slice);
    }

    fn from_path(path: &str) -> Result<Self, TimsRustError> {
        let file_reader = FrameReader::new(path)?;

        let sql_path = std::path::Path::new(path).join("analysis.tdf");
        let meta_converters = MetadataReader::new(&sql_path)?;

        let out_meta_converters = meta_converters.clone();
        let mut final_out = Self {
            indices: HashMap::new(),
            rt_converter: Some(meta_converters.rt_converter),
            mz_converter: Some(meta_converters.mz_converter),
            im_converter: Some(meta_converters.im_converter),
            metadata: Some(out_meta_converters),
            added_peaks: 0,
        };

        let out2: Result<Vec<Self>, TimsRustError> = (0..file_reader.len())
            .into_par_iter()
            .map(|id| file_reader.get(id))
            .chunks(100)
            .map(|frames| {
                let mut out = Self::new();
                for frame in frames {
                    let frame = frame?;
                    out.add_frame(frame);
                }
                Ok(out)
            })
            .collect();

        let out2 = out2?.into_iter().fold(Self::new(), |mut x, y| {
            x.fold(y);
            x
        });
        final_out.fold(out2);

        Ok(final_out)
    }

    pub fn fold(&mut self, other: Self) {
        for (qs, builder) in other.indices.into_iter() {
            self.indices
                .entry(qs)
                .and_modify(|bl| bl.fold(builder.clone()))
                .or_insert(builder);
        }
        self.added_peaks += other.added_peaks;
    }

    pub fn build(self) -> QuadSplittedTransposedIndex {
        let mut indices = HashMap::new();
        let mut flat_quad_settings = Vec::new();
        let built: Vec<(TransposedQuadIndex, Option<SingleQuadrupoleSetting>)> = self
            .indices
            .into_par_iter()
            .map(|(qs, builder)| (builder.build(), qs))
            .collect();

        let mut precursor_index: Option<TransposedQuadIndex> = None;
        for (qi, qs) in built.into_iter() {
            if let Some(qs) = qs {
                indices.insert(qs, qi);
                flat_quad_settings.push(qs);
            } else {
                precursor_index = Some(qi);
            }
        }

        flat_quad_settings.sort_by(|a, b| {
            a.ranges
                .isolation_mz
                .partial_cmp(&b.ranges.isolation_mz)
                .unwrap()
        });

        QuadSplittedTransposedIndex {
            precursor_index: precursor_index.expect("Precursor peaks should be present"),
            fragment_indices: indices,
            flat_quad_settings,
            rt_converter: self.rt_converter.unwrap(),
            mz_converter: self.mz_converter.unwrap(),
            im_converter: self.im_converter.unwrap(),
            metadata: self.metadata.unwrap(),
        }
    }
}

impl<FH: Hash + Copy + Clone + Serialize + Eq + Send + Sync>
    IndexedData<FragmentGroupIndexQuery<FH>, RawPeak> for QuadSplittedTransposedIndex
{
    fn query(&self, fragment_query: &FragmentGroupIndexQuery<FH>) -> Vec<RawPeak> {
        let precursor_mz_range = (
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
        );
        let scan_range = Some(fragment_query.precursor_query.mobility_index_range);
        let frame_index_range = Some(FrameRTTolerance::FrameIndex(
            fragment_query.precursor_query.frame_index_range,
        ));

        fragment_query
            .mz_index_ranges
            .iter()
            .flat_map(|(_, tof_range)| {
                self.query_peaks(
                    *tof_range,
                    precursor_mz_range,
                    scan_range,
                    frame_index_range,
                )
                .map(RawPeak::from)
            })
            .collect()
    }

    fn add_query<O, AG: crate::Aggregator<RawPeak, Output = O>>(
        &self,
        fragment_query: &FragmentGroupIndexQuery<FH>,
        aggregator: &mut AG,
    ) {
        let precursor_mz_range = (
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
        );
        let scan_range = Some(fragment_query.precursor_query.mobility_index_range);
        let frame_index_range = Some(FrameRTTolerance::FrameIndex(
            fragment_query.precursor_query.frame_index_range,
        ));

        fragment_query
            .mz_index_ranges
            .iter()
            .for_each(|(_, tof_range)| {
                self.query_peaks(
                    *tof_range,
                    precursor_mz_range,
                    scan_range,
                    frame_index_range,
                )
                .for_each(|peak| aggregator.add(&RawPeak::from(peak)));
            })
    }

    fn add_query_multi_group<O, AG: crate::Aggregator<RawPeak, Output = O>>(
        &self,
        fragment_queries: &[FragmentGroupIndexQuery<FH>],
        aggregator: &mut [AG],
    ) {
        fragment_queries
            .par_iter()
            .zip(aggregator.par_iter_mut())
            .for_each(|(fragment_query, agg)| {
                let precursor_mz_range = (
                    fragment_query.precursor_query.isolation_mz_range.0 as f64,
                    fragment_query.precursor_query.isolation_mz_range.1 as f64,
                );
                assert!(precursor_mz_range.0 <= precursor_mz_range.1);
                assert!(precursor_mz_range.0 > 0.0);
                let scan_range = Some(fragment_query.precursor_query.mobility_index_range);
                let frame_index_range = Some(FrameRTTolerance::FrameIndex(
                    fragment_query.precursor_query.frame_index_range,
                ));

                for (_, tof_range) in fragment_query.mz_index_ranges.clone().into_iter() {
                    self.query_peaks(tof_range, precursor_mz_range, scan_range, frame_index_range)
                        .for_each(|peak| agg.add(&RawPeak::from(peak)));
                }
            });
    }
}

// Copy pasting for now ... TODO refactor or delete
// ============================================================================

impl<FH: Eq + Hash + Copy + Serialize + Send + Sync>
    IndexedData<FragmentGroupIndexQuery<FH>, (RawPeak, FH)> for QuadSplittedTransposedIndex
{
    fn query(&self, fragment_query: &FragmentGroupIndexQuery<FH>) -> Vec<(RawPeak, FH)> {
        let precursor_mz_range = (
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
        );
        let scan_range = Some(fragment_query.precursor_query.mobility_index_range);
        let frame_index_range = Some(FrameRTTolerance::FrameIndex(
            fragment_query.precursor_query.frame_index_range,
        ));

        fragment_query
            .mz_index_ranges
            .iter()
            .flat_map(|(fh, tof_range)| {
                let out: Vec<(RawPeak, FH)> = self
                    .query_peaks(
                        *tof_range,
                        precursor_mz_range,
                        scan_range,
                        frame_index_range,
                    )
                    .map(|x| (RawPeak::from(x), *fh))
                    .collect();

                out
            })
            .collect()
    }

    fn add_query<O, AG: crate::Aggregator<(RawPeak, FH), Output = O>>(
        &self,
        fragment_query: &FragmentGroupIndexQuery<FH>,
        aggregator: &mut AG,
    ) {
        let precursor_mz_range = (
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
            fragment_query.precursor_query.isolation_mz_range.0 as f64,
        );
        let scan_range = Some(fragment_query.precursor_query.mobility_index_range);
        let frame_index_range = Some(FrameRTTolerance::FrameIndex(
            fragment_query.precursor_query.frame_index_range,
        ));

        fragment_query
            .mz_index_ranges
            .iter()
            .for_each(|(fh, tof_range)| {
                self.query_peaks(
                    *tof_range,
                    precursor_mz_range,
                    scan_range,
                    frame_index_range,
                )
                .for_each(|peak| aggregator.add(&(RawPeak::from(peak), *fh)));
            })
    }

    fn add_query_multi_group<O, AG: crate::Aggregator<(RawPeak, FH), Output = O>>(
        &self,
        fragment_queries: &[FragmentGroupIndexQuery<FH>],
        aggregator: &mut [AG],
    ) {
        fragment_queries
            .par_iter()
            .zip(aggregator.par_iter_mut())
            .for_each(|(fragment_query, agg)| {
                let precursor_mz_range = (
                    fragment_query.precursor_query.isolation_mz_range.0 as f64,
                    fragment_query.precursor_query.isolation_mz_range.1 as f64,
                );
                assert!(precursor_mz_range.0 <= precursor_mz_range.1);
                assert!(precursor_mz_range.0 > 0.0);
                let scan_range = Some(fragment_query.precursor_query.mobility_index_range);
                let frame_index_range = Some(FrameRTTolerance::FrameIndex(
                    fragment_query.precursor_query.frame_index_range,
                ));

                for (fh, tof_range) in fragment_query.mz_index_ranges.clone().into_iter() {
                    self.query_peaks(tof_range, precursor_mz_range, scan_range, frame_index_range)
                        .for_each(|peak| agg.add(&(RawPeak::from(peak), fh)));
                }
            });
    }
}

// ============================================================================

impl<FH: Copy + Clone + Serialize + Eq + Hash + Send + Sync>
    ToleranceAdapter<FragmentGroupIndexQuery<FH>, ElutionGroup<FH>>
    for QuadSplittedTransposedIndex
{
    fn query_from_elution_group(
        &self,
        tol: &dyn crate::traits::tolerance::Tolerance,
        elution_group: &ElutionGroup<FH>,
    ) -> FragmentGroupIndexQuery<FH> {
        self.queries_from_elution_elements_impl(tol, elution_group)
    }
}
