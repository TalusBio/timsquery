use crate::models::frames::expanded_frame::ExpandedFrame;
use crate::models::frames::expanded_frame::{
    expand_and_split_frame, par_expand_and_centroid_frames, ExpandedFrameSlice,
};
use crate::models::frames::single_quad_settings::SingleQuadrupoleSetting;
use log::debug;
use std::collections::HashMap;
use std::sync::Arc;
use std::time::Instant;
use timsrust::converters::{
    ConvertableDomain, Frame2RtConverter, Scan2ImConverter, Tof2MzConverter,
};
use timsrust::readers::{FrameReader, FrameReaderError, MetadataReader};
use timsrust::{QuadrupoleSettings, TimsRustError};

type QuadSettingsIndex = usize;

#[derive(Debug, Default)]
struct ExpandedRawFrameIndex {
    bundled_ms1_frames: Vec<ExpandedFrameSlice>,
    bundled_frames: HashMap<SingleQuadrupoleSetting, Vec<ExpandedFrameSlice>>,
    flat_quad_settings: Vec<SingleQuadrupoleSetting>,
    rt_converter: Frame2RtConverter,
    pub mz_converter: Tof2MzConverter,
    pub im_converter: Scan2ImConverter,
}

impl ExpandedRawFrameIndex {
    pub fn from_tdf_path(path: &str) -> Result<Self, TimsRustError> {
        let file_reader = FrameReader::new(path)?;

        let sql_path = std::path::Path::new(path).join("analysis.tdf");
        let meta_converters = MetadataReader::new(&sql_path)?;

        let st = Instant::now();
        let all_frames = file_reader.get_all().into_iter().flatten().collect();
        let read_elap = st.elapsed();
        debug!("Reading all frames took {:#?}", read_elap);
        let st = Instant::now();
        let centroided_split_frames = par_expand_and_centroid_frames(
            all_frames,
            1.5,
            15.0,
            &meta_converters.im_converter,
            &meta_converters.mz_converter,
        );
        let centroided_elap = st.elapsed();
        debug!("Centroiding took {:#?}", centroided_elap);

        let mut out_ms2_frames = HashMap::new();
        let mut out_ms1_frames: Option<Vec<ExpandedFrameSlice>> = None;

        centroided_split_frames
            .into_iter()
            .for_each(|(q, frameslices)| match q {
                None => {
                    out_ms1_frames = Some(frameslices);
                }
                Some(q) => {
                    out_ms2_frames.insert(q, frameslices);
                }
            });

        let mut flat_quad_settings = out_ms2_frames.keys().cloned().collect();

        let out = Self {
            bundled_ms1_frames: out_ms1_frames.expect("At least one ms1 frame should be present"),
            bundled_frames: out_ms2_frames,
            flat_quad_settings,
            rt_converter: meta_converters.rt_converter,
            mz_converter: meta_converters.mz_converter,
            im_converter: meta_converters.im_converter,
        };

        Ok(out)
    }
}
