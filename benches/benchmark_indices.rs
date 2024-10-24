use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use serde::Serialize;
use std::collections::HashMap;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};
use timsquery::{
    models::{
        aggregators::RawPeakIntensityAggregator,
        indices::{
            expanded_raw_index::ExpandedRawFrameIndex, raw_file_index::RawFileIndex,
            transposed_quad_index::QuadSplittedTransposedIndex,
        },
    },
    queriable_tims_data::queriable_tims_data::query_multi_group,
    traits::tolerance::{
        DefaultTolerance, MobilityTolerance, MzToleramce, QuadTolerance, RtTolerance,
    },
    ElutionGroup,
};

const NUM_ELUTION_GROUPS: usize = 1000;
const NUM_ITERATIONS: usize = 1;

#[derive(Debug, Serialize)]
struct BenchmarkIteration {
    iteration: usize,
    duration_seconds: f64,
}

#[derive(Debug, Serialize)]
struct BenchmarkResult {
    name: String,
    context: String,
    iterations: Vec<BenchmarkIteration>,
    mean_duration_seconds: f64,
    mean_duration_human_readable: String,
    note: Option<String>,
}

#[derive(Debug, Serialize)]
struct BenchmarkReport {
    settings: BenchmarkSettings,
    results: Vec<BenchmarkResult>,
    metadata: BenchmarkMetadata,
    full_benchmark_time_seconds: f64,
}

#[derive(Debug, Serialize)]
struct BenchmarkSettings {
    num_elution_groups: usize,
    num_iterations: usize,
}

#[derive(Debug, Serialize)]
struct BenchmarkMetadata {
    basename: String,
}

fn duration_to_seconds(duration: Duration) -> f64 {
    duration.as_secs() as f64 + duration.subsec_nanos() as f64 * 1e-9
}

fn get_file_from_env() -> (String, String) {
    let raw_file_path =
        std::env::var("TIMS_DATA_FILE").expect("TIMS_DATA_FILE environment variable not set");

    let basename = std::path::Path::new(&raw_file_path)
        .file_stem()
        .unwrap()
        .to_str()
        .unwrap()
        .to_string();

    (raw_file_path, basename)
}

fn build_elution_groups() -> Vec<ElutionGroup<u64>> {
    const NUM_FRAGMENTS: usize = 10;
    const MAX_RT: f32 = 22.0 * 60.0;
    const MAX_MOBILITY: f32 = 1.5;
    const MIN_MOBILITY: f32 = 0.5;
    const MAX_MZ: f64 = 1000.0;
    const MIN_MZ: f64 = 300.0;

    let mut out_egs = Vec::with_capacity(NUM_ELUTION_GROUPS);
    let mut rng = ChaCha8Rng::seed_from_u64(43u64);

    for i in 1..NUM_ELUTION_GROUPS {
        let rt = rng.gen::<f32>() * MAX_RT;
        let mobility = rng.gen::<f32>() * (MAX_MOBILITY - MIN_MOBILITY) + MIN_MOBILITY;
        let mz = rng.gen::<f64>() * (MAX_MZ - MIN_MZ) + MIN_MZ;

        let mut fragment_mzs = HashMap::with_capacity(NUM_FRAGMENTS);
        for ii in 0..NUM_FRAGMENTS {
            let fragment_mz = rng.gen::<f64>() * (MAX_MZ - MIN_MZ) + MIN_MZ;
            fragment_mzs.insert(ii as u64, fragment_mz);
        }

        out_egs.push(ElutionGroup {
            id: i as u64,
            rt_seconds: rt,
            mobility,
            precursor_mz: mz,
            precursor_charge: 2,
            fragment_mzs,
        });
    }
    out_egs
}

fn with_benchmark(name: &str, context: &str, f: impl Fn() -> Option<String>) -> BenchmarkResult {
    let mut iterations = Vec::with_capacity(NUM_ITERATIONS);
    let mut durations = Vec::with_capacity(NUM_ITERATIONS);
    let mut extras = None;
    let glob_start = Instant::now();
    let mut i = 0;
    while (glob_start.elapsed().as_millis() < 1000) || (i < NUM_ITERATIONS) {
        println!("{name} iteration {i}");
        let start = Instant::now();
        if let Some(out) = f() {
            extras = Some(out);
        }

        let duration = start.elapsed();

        iterations.push(BenchmarkIteration {
            iteration: i + 1,
            duration_seconds: duration_to_seconds(duration),
        });
        durations.push(duration);
        i += 1;
    }

    let mean = iterations.iter().map(|i| i.duration_seconds).sum::<f64>() / (i as f64);
    let avg_duration: Duration = durations.iter().sum::<Duration>() / i as u32;
    let mean_duration_human_readable = format!("{:?}", avg_duration);

    BenchmarkResult {
        name: name.to_string(),
        context: context.to_string(),
        iterations,
        mean_duration_seconds: mean,
        mean_duration_human_readable,
        note: extras,
    }
}

fn run_encoding_benchmark(raw_file_path: &str) -> Vec<BenchmarkResult> {
    let mut out = vec![];
    let rfi = with_benchmark("RawFileIndex", "Encoding", || {
        RawFileIndex::from_path(raw_file_path).unwrap();
        None
    });
    out.push(rfi);
    let erfic = with_benchmark("ExpandedRawFileIndexCentroided", "Encoding", || {
        ExpandedRawFrameIndex::from_path(raw_file_path).unwrap();
        None
    });
    out.push(erfic);
    let erfi = with_benchmark("ExpandedRawFileIndex", "Encoding", || {
        ExpandedRawFrameIndex::from_path(raw_file_path).unwrap();
        None
    });
    out.push(erfi);
    let tqi = with_benchmark("TransposedQuadIndex", "Encoding", || {
        QuadSplittedTransposedIndex::from_path(raw_file_path).unwrap();
        None
    });
    out.push(tqi);

    let tqic = with_benchmark("TransposedQuadIndexCentroided", "Encoding", || {
        QuadSplittedTransposedIndex::from_path_centroided(raw_file_path).unwrap();
        None
    });
    out.push(tqic);

    out
}

fn run_batch_access_benchmark(raw_file_path: &str) -> Vec<BenchmarkResult> {
    let mut out = vec![];
    let query_groups = build_elution_groups();
    let tolerance_with_rt = DefaultTolerance {
        ms: MzToleramce::Ppm((20.0, 20.0)),
        rt: RtTolerance::Absolute((5.0, 5.0)),
        mobility: MobilityTolerance::Pct((3.0, 3.0)),
        quad: QuadTolerance::Absolute((0.1, 0.1, 1)),
    };
    let tolerance_with_nort = DefaultTolerance {
        ms: MzToleramce::Ppm((20.0, 20.0)),
        rt: RtTolerance::None,
        mobility: MobilityTolerance::Pct((3.0, 3.0)),
        quad: QuadTolerance::Absolute((0.1, 0.1, 1)),
    };
    let tolerances = [
        (tolerance_with_rt, "narrow_rt"),
        (tolerance_with_nort, "full_rt"),
    ];

    // TODO: Refactor this ... there is a lot of code duplication but over different types ... int
    // ion thoty this is a good place for a macro.

    let rfi_index = RawFileIndex::from_path(raw_file_path).unwrap();
    for (tolerance, tol_name) in tolerances.clone() {
        if tol_name == "full_rt" {
            println!("Skipping full_rt");
            continue;
        }
        let rfi = with_benchmark("RawFileIndex", &format!("BatchAccess_{}", tol_name), || {
            let tmp = query_multi_group(
                &rfi_index,
                &rfi_index,
                &tolerance,
                &query_groups,
                &RawPeakIntensityAggregator::new,
            );
            let tot: u64 = tmp.into_iter().sum();
            let out = format!("RawFileIndex::query_multi_group aggregated {}", tot);
            println!("{}", out);
            Some(out)
        });
        out.push(rfi);
    }
    std::mem::drop(rfi_index);

    let erfi_index = ExpandedRawFrameIndex::from_path(raw_file_path).unwrap();
    for (tolerance, tol_name) in tolerances.clone() {
        let erfi = with_benchmark(
            "ExpandedRawFileIndex",
            &format!("BatchAccess_{}", tol_name),
            || {
                let tmp = query_multi_group(
                    &erfi_index,
                    &erfi_index,
                    &tolerance,
                    &query_groups,
                    &RawPeakIntensityAggregator::new,
                );
                let tot: u64 = tmp.into_iter().sum();
                let out = format!("ExpandedRawFileIndex::query_multi_group aggregated {}", tot);
                println!("{}", out);
                Some(out)
            },
        );
        out.push(erfi);
    }
    std::mem::drop(erfi_index);

    let erfic_index = ExpandedRawFrameIndex::from_path_centroided(raw_file_path).unwrap();
    for (tolerance, tol_name) in tolerances.clone() {
        let erfic = with_benchmark(
            "ExpandedRawFileIndexCentroided",
            &format!("BatchAccess_{}", tol_name),
            || {
                let tmp = query_multi_group(
                    &erfic_index,
                    &erfic_index,
                    &tolerance,
                    &query_groups,
                    &RawPeakIntensityAggregator::new,
                );
                let tot: u64 = tmp.into_iter().sum();
                let out = format!(
                    "ExpandedRawFileIndexCentroided::query_multi_group aggregated {}",
                    tot
                );
                println!("{}", out);
                Some(out)
            },
        );
        out.push(erfic);
    }
    std::mem::drop(erfic_index);

    let tqi_index = QuadSplittedTransposedIndex::from_path(raw_file_path).unwrap();
    for (tolerance, tol_name) in tolerances.clone() {
        let tqi = with_benchmark(
            "TransposedQuadIndex",
            &format!("BatchAccess_{}", tol_name),
            || {
                let tmp = query_multi_group(
                    &tqi_index,
                    &tqi_index,
                    &tolerance,
                    &query_groups,
                    &RawPeakIntensityAggregator::new,
                );
                let tot: u64 = tmp.into_iter().sum();
                let out = format!("TransposedQuadIndex::query_multi_group aggregated {}", tot);
                println!("{}", out);
                Some(out)
            },
        );
        out.push(tqi);
    }
    std::mem::drop(tqi_index);

    let tqic_index = QuadSplittedTransposedIndex::from_path_centroided(raw_file_path).unwrap();
    for (tolerance, tol_name) in tolerances.clone() {
        let tqi = with_benchmark(
            "TransposedQuadIndexCentroided",
            &format!("BatchAccess_{}", tol_name),
            || {
                let tmp = query_multi_group(
                    &tqic_index,
                    &tqic_index,
                    &tolerance,
                    &query_groups,
                    &RawPeakIntensityAggregator::new,
                );
                let tot: u64 = tmp.into_iter().sum();
                let out = format!("TransposedQuadIndex::query_multi_group aggregated {}", tot);
                println!("{}", out);
                Some(out)
            },
        );
        out.push(tqi);
    }
    std::mem::drop(tqic_index);
    out
}

fn write_results(
    report: BenchmarkReport,
    basename: &str,
    parent: &Path,
) -> std::io::Result<(PathBuf, String)> {
    let filepath = parent.join(format!("benchmark_results_{}.json", basename));
    let file = File::create(&filepath)?;
    serde_json::to_writer_pretty(file, &report)?;
    let out = serde_json::to_string_pretty(&report)?;
    Ok((filepath, out))
}

fn main() {
    env_logger::init();
    let st = Instant::now();
    let (raw_file_path, basename) = get_file_from_env();
    let file_parent = std::path::Path::new(&raw_file_path)
        .parent()
        .expect("Expected to find a parent directory");
    let mut all_results: Vec<BenchmarkResult> = vec![];

    println!(
        "Run Settings: Elution groups={}, Iterations={}, basename={}",
        NUM_ELUTION_GROUPS, NUM_ITERATIONS, basename
    );

    println!("Running encoding benchmarks...");
    let encoding_results = run_encoding_benchmark(&raw_file_path);
    all_results.extend(encoding_results);

    println!("Running batch access benchmarks...");
    let batch_results = run_batch_access_benchmark(&raw_file_path);
    all_results.extend(batch_results);
    let report = BenchmarkReport {
        settings: BenchmarkSettings {
            num_elution_groups: NUM_ELUTION_GROUPS,
            num_iterations: NUM_ITERATIONS,
        },
        results: all_results,
        metadata: BenchmarkMetadata {
            basename: basename.to_string(),
        },
        full_benchmark_time_seconds: st.elapsed().as_secs_f64(),
    };

    match write_results(report, &basename, file_parent) {
        Ok((fp, res)) => println!("Results written to {} \n{}", fp.to_string_lossy(), res),
        Err(e) => eprintln!("Error writing results: {}", e),
    }
}
