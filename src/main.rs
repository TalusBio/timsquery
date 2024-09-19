use timsquery::models::elution_group::ElutionGroup;
use timsquery::queriable_tims_data::queriable_tims_data::query_multi_group;
use timsquery::traits::tolerance::DefaultTolerance;
use timsquery::Aggregator;
use timsquery::{
    models::aggregators::RawPeakIntensityAggregator, models::indices::raw_file_index::RawFileIndex,
    models::indices::transposed_quad_index::QuadSplittedTransposedIndex,
};

use timsquery::traits::tolerance::{MobilityTolerance, MzToleramce, QuadTolerance, RtTolerance};

use clap::{Parser, Subcommand};
use serde::{Deserialize, Serialize};

// Read json with tolerance settings
// Read json with elution groups
// Load index
// Query index
// Serialize results

fn template_elution_groups(num: usize) -> Vec<ElutionGroup> {
    let mut egs = Vec::with_capacity(num);
    for i in 1..num {
        let rt = 300.0 + (i as f32 * 10.0);
        let mobility = 1.0 + (i as f32 * 0.01);
        let mz = 1000.0 + (i as f64 * 10.0);
        let precursor_charge = 2;
        let fragment_mzs = Some(vec![mz]);
        let fragment_charges = Some(vec![precursor_charge]);
        egs.push(ElutionGroup {
            id: i as u64,
            rt_seconds: rt,
            mobility,
            precursor_mz: mz,
            precursor_charge,
            fragment_mzs,
            fragment_charges,
        });
    }
    egs
}

fn template_tolerance_settings() -> DefaultTolerance {
    DefaultTolerance {
        ms: MzToleramce::Ppm((20.0, 20.0)),
        rt: RtTolerance::Absolute((5.0, 5.0)),
        mobility: MobilityTolerance::Pct((30.0, 30.0)),
        quad: QuadTolerance::Absolute((0.1, 0.1, 1)),
    }
}

#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, clap::ValueEnum)]
enum PossibleAggregator {
    #[default]
    RawPeakIntensityAggregator,
}

#[derive(Debug, Default, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, clap::ValueEnum)]
enum PossibleIndex {
    #[default]
    Unspecified,
    RawFileIndex,
    TransposedQuadIndex,
}

#[derive(Parser, Debug)]
struct QueryIndexArgs {
    /// The path to the raw file to query.
    #[arg(short, long)]
    raw_file_path: String,

    /// The path to the json file with the tolerance settings.
    #[arg(short, long)]
    tolerance_settings_path: String,

    /// The path to the json file with the elution groups.
    #[arg(short, long)]
    elution_groups_path: String,

    /// The path to the output files.
    #[arg(short, long)]
    output_path: String,

    // Whether the output json should be pretty printed.
    #[arg(short, long)]
    pretty: bool,

    // The aggregator to use.
    #[arg(short, long, default_value_t, value_enum)]
    aggregator: PossibleAggregator,

    // The index to use.
    #[arg(short, long, value_enum)]
    index: PossibleIndex,
}

#[derive(Parser, Debug)]
struct WriteTemplateArgs {
    /// The path to the output files.
    #[arg(short, long)]
    output_path: String,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Query the index.
    QueryIndex(QueryIndexArgs),
    WriteTemplate(WriteTemplateArgs),
}

#[derive(Debug, Serialize, Deserialize)]
struct ElutionGroupResults {
    elution_group: ElutionGroup,
    result: u64,
}

fn main() {
    env_logger::init();
    let args = Args::parse();

    match args.command {
        Some(Commands::QueryIndex(args)) => main_query_index(args),
        Some(Commands::WriteTemplate(args)) => main_write_template(args),
        None => {
            println!("No command provided");
        }
    }
}

fn main_write_template(args: WriteTemplateArgs) {
    let output_path = args.output_path;
    let egs = template_elution_groups(10);
    let tolerance = template_tolerance_settings();

    // Serialize both and write as files in the output path.
    // Do pretty serialization.
    let egs_json = serde_json::to_string_pretty(&egs).unwrap();
    let tolerance_json = serde_json::to_string_pretty(&tolerance).unwrap();

    let put_path = std::path::Path::new(&output_path);
    std::fs::create_dir_all(put_path).unwrap();
    println!("Writing to {}", put_path.display());
    let egs_json_path = put_path.join("elution_groups.json");
    let tolerance_json_path = put_path.join("tolerance_settings.json");
    std::fs::write(egs_json_path.clone(), egs_json).unwrap();
    std::fs::write(tolerance_json_path.clone(), tolerance_json).unwrap();
    println!(
        "use as `timsquery query-index --pretty --output-path '.' --raw-file-path 'your_file.d' --tolerance-settings-path {:#?} --elution-groups-path {:#?}`",
        tolerance_json_path,
        egs_json_path,
    );
}

fn main_query_index(args: QueryIndexArgs) {
    let raw_file_path = args.raw_file_path;
    let tolerance_settings_path = args.tolerance_settings_path;
    let elution_groups_path = args.elution_groups_path;
    let output_path = args.output_path;
    let index_use = args.index;

    let tolerance_settings: DefaultTolerance =
        serde_json::from_str(&std::fs::read_to_string(&tolerance_settings_path).unwrap()).unwrap();
    let elution_groups: Vec<ElutionGroup> =
        serde_json::from_str(&std::fs::read_to_string(&elution_groups_path).unwrap()).unwrap();

    let aggregator_factory = |_id| RawPeakIntensityAggregator { intensity: 0 };
    let foo = if (elution_groups.len() > 10) || index_use == PossibleIndex::TransposedQuadIndex {
        let index = QuadSplittedTransposedIndex::from_path(&(raw_file_path.clone())).unwrap();
        query_multi_group(
            &index,
            &index,
            &tolerance_settings,
            &elution_groups,
            &aggregator_factory,
        )
    } else {
        let index = RawFileIndex::from_path(&(raw_file_path.clone())).unwrap();
        query_multi_group(
            &index,
            &index,
            &tolerance_settings,
            &elution_groups,
            &aggregator_factory,
        )
    };

    let mut out = Vec::new();
    for (res, eg) in foo.into_iter().zip(elution_groups) {
        out.push(ElutionGroupResults {
            elution_group: eg,
            result: res.finalize(),
        });
    }

    let put_path = std::path::Path::new(&output_path).join("results.json");
    std::fs::create_dir_all(put_path.parent().unwrap()).unwrap();
    println!("Writing to {}", put_path.display());

    let serialized = if args.pretty {
        println!("Pretty printing enabled");
        serde_json::to_string_pretty(&out).unwrap()
    } else {
        serde_json::to_string(&out).unwrap()
    };
    std::fs::write(put_path, serialized).unwrap();
}

// fn main() {
//     println!("Hello, world!");
//     let qst_file_index = QuadSplittedTransposedIndex::from_path(&(raw_file_path.clone())).unwrap();
//     let tolerance = DefaultTolerance::default();
//     let aggregator_factory = |id| RawPeakIntensityAggregator { intensity: 0 };
//     let foo = query_multi_group(
//         &qst_file_index,
//         &qst_file_index,
//         &tolerance,
//         &query_groups,
//         &aggregator_factory,
//     );
// }
