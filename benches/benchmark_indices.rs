use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;
use serde::Serialize;
use std::collections::HashMap;
use std::fs::File;
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
    traits::tolerance::DefaultTolerance,
    ElutionGroup,
};

const NUM_ELUTION_GROUPS: usize = 500;
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
}

#[derive(Debug, Serialize)]
struct BenchmarkReport {
    settings: BenchmarkSettings,
    results: Vec<BenchmarkResult>,
    metadata: BenchmarkMetadata,
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

fn with_benchmark(name: &str, context: &str, f: impl Fn()) -> BenchmarkResult {
    let mut iterations = Vec::with_capacity(NUM_ITERATIONS);
    for i in 0..NUM_ITERATIONS {
        println!("{name} iteration {i}");
        let start = Instant::now();
        f();
        let duration = start.elapsed();

        iterations.push(BenchmarkIteration {
            iteration: i + 1,
            duration_seconds: duration_to_seconds(duration),
        });
    }

    let mean = iterations.iter().map(|i| i.duration_seconds).sum::<f64>() / NUM_ITERATIONS as f64;
    BenchmarkResult {
        name: name.to_string(),
        context: context.to_string(),
        iterations,
        mean_duration_seconds: mean,
    }
}

fn run_encoding_benchmark(raw_file_path: &str) -> Vec<BenchmarkResult> {
    let rfi = with_benchmark("RawFileIndex", "Encoding", || {
        RawFileIndex::from_path(raw_file_path).unwrap();
    });
    let erfi = with_benchmark("ExpandedRawFileIndex", "Encoding", || {
        ExpandedRawFrameIndex::from_path(raw_file_path).unwrap();
    });
    let tqi = with_benchmark("TransposedQuadIndex", "Encoding", || {
        QuadSplittedTransposedIndex::from_path(raw_file_path).unwrap();
    });

    vec![tqi, rfi, erfi]
}

fn run_batch_access_benchmark(raw_file_path: &str) -> Vec<BenchmarkResult> {
    let query_groups = build_elution_groups();
    let tolerance = DefaultTolerance::default();

    let rfi_index = RawFileIndex::from_path(raw_file_path).unwrap();
    let rfi = with_benchmark("RawFileIndex", "BatchAccess", || {
        let _ = query_multi_group(
            &rfi_index,
            &rfi_index,
            &tolerance,
            &query_groups,
            &RawPeakIntensityAggregator::new,
        );
    });
    std::mem::drop(rfi_index);

    let erfi_index = ExpandedRawFrameIndex::from_path(raw_file_path).unwrap();
    let erfi = with_benchmark("ExpandedRawFileIndex", "BatchAccess", || {
        let _ = query_multi_group(
            &erfi_index,
            &erfi_index,
            &tolerance,
            &query_groups,
            &RawPeakIntensityAggregator::new,
        );
    });
    std::mem::drop(erfi_index);

    let tqi_index = QuadSplittedTransposedIndex::from_path(raw_file_path).unwrap();
    let tqi = with_benchmark("TransposedQuadIndex", "BatchAccess", || {
        let _ = query_multi_group(
            &tqi_index,
            &tqi_index,
            &tolerance,
            &query_groups,
            &RawPeakIntensityAggregator::new,
        );
    });
    std::mem::drop(tqi_index);

    vec![tqi, rfi, erfi]
}

fn write_results(results: Vec<BenchmarkResult>, basename: &str) -> std::io::Result<()> {
    let report = BenchmarkReport {
        settings: BenchmarkSettings {
            num_elution_groups: NUM_ELUTION_GROUPS,
            num_iterations: NUM_ITERATIONS,
        },
        results,
        metadata: BenchmarkMetadata {
            basename: basename.to_string(),
        },
    };

    let file = File::create(format!("benchmark_results_{}.json", basename))?;
    serde_json::to_writer_pretty(file, &report)?;
    Ok(())
}

fn main() {
    let (raw_file_path, basename) = get_file_from_env();
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

    match write_results(all_results, &basename) {
        Ok(_) => println!("Results written to benchmark_results_{}.json", basename),
        Err(e) => eprintln!("Error writing results: {}", e),
    }
}
