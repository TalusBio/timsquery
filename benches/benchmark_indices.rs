use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion};
use rand::{Rng, SeedableRng};
use rand_chacha::ChaCha8Rng;

use std::hint::black_box;
use timsquery::{
    models::{
        indices::raw_file_index::{RawFileIndex, RawPeakIntensityAggregator},
        queries::{
            FragmentGroupIndexQuery, NaturalFragmentQuery, NaturalPrecursorQuery,
            PrecursorIndexQuery,
        },
    },
    queriable_tims_data::queriable_tims_data::{query_indexed, query_multi_group},
    traits::tolerance::DefaultTolerance,
    ElutionGroup, QueriableTimsData,
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

    for i in (1..NUM_ELUTION_GROUPS) {
        let rt = rng.gen::<f32>() * MAX_RT;
        let mobility = rng.gen::<f32>() * (MAX_MOBILITY - MIN_MOBILITY) + MIN_MOBILITY;
        let mz = rng.gen::<f64>() * (MAX_MZ - MIN_MZ) + MIN_MZ;

        let mut fragment_mzs = Vec::with_capacity(NUM_FRAGMENTS);
        let mut fragment_charges = Vec::with_capacity(NUM_FRAGMENTS);

        for _ in 0..NUM_FRAGMENTS {
            let fragment_mz = rng.gen::<f64>() * (MAX_MZ - MIN_MZ) + MIN_MZ;
            let fragment_charge = rng.gen::<i8>() * 3 + 1;
            fragment_mzs.push(fragment_mz);
            fragment_charges.push(fragment_charge);
        }

        let precursor_charge = rng.gen::<u8>() * 3 + 1;

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
    group.bench_function(BenchmarkId::new("RawFileIndex", basename.clone()), |b| {
        b.iter_batched(
            || {},
            |()| RawFileIndex::from_path(&black_box(raw_file_path.clone())).unwrap(),
            BatchSize::SmallInput,
        )
    });

    group.finish();
}

fn thoughput_benchmark_random(c: &mut Criterion) {
    let (raw_file_path, basename) = get_file_from_env();
    let mut group = c.benchmark_group("RandomAccessThroughput");
    group.significance_level(0.05).sample_size(10);
    group.throughput(criterion::Throughput::Elements(NUM_ELUTION_GROUPS as u64));
    group.bench_function(BenchmarkId::new("RawFileIndex", basename.clone()), |b| {
        b.iter_batched(
            || {
                (
                    RawFileIndex::from_path(&(raw_file_path.clone())).unwrap(),
                    build_elution_groups(raw_file_path.clone()),
                    DefaultTolerance::default(),
                )
            },
            |(raw_file_index, query_groups, tolerance)| {
                let aggregator_factory = |id| RawPeakIntensityAggregator { intensity: 0 };
                let local_lambda = |elution_group: &ElutionGroup| {
                    query_indexed(
                        &raw_file_index,
                        &aggregator_factory,
                        &raw_file_index,
                        &tolerance,
                        &elution_group,
                    )
                };
                for elution_group in query_groups {
                    let foo = local_lambda(&elution_group);
                    black_box((|foo| false)(foo));
                }
            },
            BatchSize::PerIteration,
        )
    });

    group.finish();
}

fn thoughput_benchmark_optim(c: &mut Criterion) {
    let (raw_file_path, basename) = get_file_from_env();
    let mut group = c.benchmark_group("BatchAccessThroughput");
    group.significance_level(0.05).sample_size(10);
    group.throughput(criterion::Throughput::Elements(NUM_ELUTION_GROUPS as u64));
    group.bench_function(BenchmarkId::new("RawFileIndex", basename.clone()), |b| {
        b.iter_batched(
            || {
                (
                    RawFileIndex::from_path(&(raw_file_path.clone())).unwrap(),
                    build_elution_groups(raw_file_path.clone()),
                    DefaultTolerance::default(),
                )
            },
            |(raw_file_index, query_groups, tolerance)| {
                let aggregator_factory = |id| RawPeakIntensityAggregator { intensity: 0 };
                let foo = query_multi_group(
                    &raw_file_index,
                    &raw_file_index,
                    &tolerance,
                    &query_groups,
                    &aggregator_factory,
                );
                black_box((|foo| false)(foo));
            },
            BatchSize::PerIteration,
        )
    });

    group.finish();
}

criterion_group!(
    benches,
    criterion_benchmark,
    thoughput_benchmark_optim,
    thoughput_benchmark_random
);
criterion_main!(benches);
