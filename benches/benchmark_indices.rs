use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use std::hint::black_box;
use timsquery::{
    models::{
        aggregators::RawPeakIntensityAggregator,
        indices::{
            raw_file_index::RawFileIndex, transposed_quad_index::QuadSplittedTransposedIndex,
        },
    },
    queriable_tims_data::queriable_tims_data::{query_indexed, query_multi_group},
    traits::tolerance::DefaultTolerance,
    ElutionGroup,
};

fn get_file_from_env() -> (String, String) {
    // Read from environment variable the raw file path to use
    // env var name: TIMS_DATA_FILE
    let raw_file_path = std::env::var("TIMS_DATA_FILE");
    let raw_file_path = match raw_file_path {
        Ok(path) => path,
        Err(_) => panic!("TIMS_DATA_FILE environment variable not set"),
    };

    let basename = std::path::Path::new(&raw_file_path)
        .file_stem()
        .unwrap()
        .to_str()
        .unwrap()
        .to_string();

    (raw_file_path, basename)
}

const NUM_ELUTION_GROUPS: usize = 500;
fn build_elution_groups(raw_file_path: String) -> Vec<ElutionGroup> {
    const NUM_FRAGMENTS: usize = 10;
    const MAX_RT: f32 = 22.0 * 60.0;
    const MAX_MOBILITY: f32 = 1.5;
    const MIN_MOBILITY: f32 = 0.5;
    const MAX_MZ: f64 = 1000.0;
    const MIN_MZ: f64 = 300.0;

    let mut out_egs: Vec<ElutionGroup> = Vec::with_capacity(NUM_ELUTION_GROUPS);
    let mut rng = ChaCha8Rng::seed_from_u64(43u64);

    for i in 1..NUM_ELUTION_GROUPS {
        // Rand f32/64 are number from 0-1
        let rt = rng.gen::<f32>() * MAX_RT;
        let mobility = rng.gen::<f32>() * (MAX_MOBILITY - MIN_MOBILITY) + MIN_MOBILITY;
        let mz = rng.gen::<f64>() * (MAX_MZ - MIN_MZ) + MIN_MZ;

        let mut fragment_mzs = Vec::with_capacity(NUM_FRAGMENTS);
        let mut fragment_charges = Vec::with_capacity(NUM_FRAGMENTS);

        for _ in 0..NUM_FRAGMENTS {
            let fragment_mz = rng.gen::<f64>() * (MAX_MZ - MIN_MZ) + MIN_MZ;
            // let fragment_charge = rng.gen::<u8>() * 3 + 1;
            let fragment_charge = 1;
            fragment_mzs.push(fragment_mz);
            fragment_charges.push(fragment_charge);
        }

        // rand u8 is number from 0-255 ... which is not amazing for us ...
        // let precursor_charge = rng.gen::<u8>() * 3 + 1;
        let precursor_charge = 2;

        out_egs.push(ElutionGroup {
            id: i as u64,
            rt_seconds: rt,
            mobility,
            precursor_mz: mz,
            precursor_charge,
            fragment_mzs: Some(fragment_mzs),
            fragment_charges: Some(fragment_charges),
        });
    }
    out_egs
}

fn criterion_benchmark(c: &mut Criterion) {
    let (raw_file_path, basename) = get_file_from_env();
    let mut group = c.benchmark_group("Encoding Time");

    group.sample_size(10);
    group.bench_function(
        BenchmarkId::new("TransposedQuadIndex", basename.clone()),
        |b| {
            b.iter_batched(
                || {},
                |()| {
                    QuadSplittedTransposedIndex::from_path(&black_box(raw_file_path.clone()))
                        .unwrap()
                },
                BatchSize::SmallInput,
            )
        },
    );
    group.bench_function(BenchmarkId::new("RawFileIndex", basename.clone()), |b| {
        b.iter_batched(
            || {},
            |()| RawFileIndex::from_path(&black_box(raw_file_path.clone())).unwrap(),
            BatchSize::SmallInput,
        )
    });

    group.finish();
}

macro_rules! add_bench_random {
    ($group:expr, $raw_file_path:expr, $basename:expr, $name:literal, $index_type:ty, $tolerance_type:ty, ) => {
        $group.bench_function(BenchmarkId::new($name, $basename.clone()), |b| {
            b.iter_batched(
                || {
                    (
                        <$index_type>::from_path(&($raw_file_path.clone())).unwrap(),
                        build_elution_groups($raw_file_path.clone()),
                        <$tolerance_type>::default(),
                    )
                },
                |(index, query_groups, tolerance)| {
                    let local_lambda = |elution_group: &ElutionGroup| {
                        query_indexed(
                            &index,
                            &RawPeakIntensityAggregator::new,
                            &index,
                            &tolerance,
                            &elution_group,
                        )
                    };
                    for elution_group in query_groups {
                        let query_results = local_lambda(&elution_group);
                        black_box((|_query_results| false)(query_results));
                    }
                },
                BatchSize::PerIteration,
            )
        });
    };
}
macro_rules! add_bench_optim {
    ($group:expr, $raw_file_path:expr, $basename:expr, $name:literal, $index_type:ty, $tolerance_type:ty, $query_func:expr,) => {
        $group.bench_function(BenchmarkId::new($name, $basename.clone()), |b| {
            b.iter_batched(
                || {
                    (
                        <$index_type>::from_path(&($raw_file_path.clone())).unwrap(),
                        build_elution_groups($raw_file_path.clone()),
                        <$tolerance_type>::default(),
                    )
                },
                |(index, query_groups, tolerance)| {
                    let qr = query_multi_group(
                        &index,
                        &index,
                        &tolerance,
                        &query_groups,
                        &RawPeakIntensityAggregator::new,
                    );
                    black_box((|_qr| false)(qr));
                },
                BatchSize::PerIteration,
            )
        });
    };
}

fn thoughput_benchmark_random(c: &mut Criterion) {
    let (raw_file_path, basename) = get_file_from_env();
    let mut group = c.benchmark_group("RandomAccessThroughput");
    group.significance_level(0.05).sample_size(10);
    group.throughput(criterion::Throughput::Elements(NUM_ELUTION_GROUPS as u64));

    add_bench_random!(
        group,
        raw_file_path,
        basename,
        "TransposedQuadIndex",
        QuadSplittedTransposedIndex,
        DefaultTolerance,
    );

    add_bench_random!(
        group,
        raw_file_path,
        basename,
        "RawFileIndex",
        RawFileIndex,
        DefaultTolerance,
    );

    group.finish();
}

fn thoughput_benchmark_optim(c: &mut Criterion) {
    env_logger::init();
    let (raw_file_path, basename) = get_file_from_env();
    let mut group = c.benchmark_group("BatchAccessThroughput");
    group.significance_level(0.05).sample_size(10);
    group.throughput(criterion::Throughput::Elements(NUM_ELUTION_GROUPS as u64));

    add_bench_optim!(
        group,
        raw_file_path,
        basename,
        "TransposedQuadIndex",
        QuadSplittedTransposedIndex,
        DefaultTolerance,
        query_multi_group,
    );

    add_bench_optim!(
        group,
        raw_file_path,
        basename,
        "RawFileIndex",
        RawFileIndex,
        DefaultTolerance,
        query_multi_group,
    );

    group.finish();
}

criterion_group!(
    benches,
    thoughput_benchmark_optim,
    criterion_benchmark,
    thoughput_benchmark_random,
);
criterion_main!(benches);
