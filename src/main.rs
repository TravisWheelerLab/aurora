mod alignment;
mod alphabet;
mod annotation;
mod assembly;
mod balanced_tree;
mod chunks;
mod confidence;
mod history_tracing;
mod join_estimation;
mod matrix;
mod pipeline;
mod score_params;
mod segment_groups;
mod segments;
mod statistics;
mod substitution_matrix;
mod support;

// Note: This will be used for components/archetectures additions soon...
#[allow(dead_code)]
mod union_find;

mod trace_statistics;
mod util;
mod viterbi;
mod viz;
mod windowed_scores;

use std::{
    collections::HashMap,
    fs::{self, create_dir_all, File},
    io::{BufRead, BufReader, BufWriter, Write},
    path::PathBuf,
};

use alignment::AlignmentData;
use chunks::ProximityGroup;

use anyhow::{Ok, Result};
use clap::{Args, Parser};
use itertools::Itertools;
use rayon::prelude::*;
use viz::VizConstraint;

use crate::{
    annotation::AmbiguousAnnotation,
    chunks::validate_groups,
    pipeline::{run_history_trace, run_naive_trace, NaiveTraceResults},
    trace_statistics::{trace_statistics, OccuranceCountingMode},
    viz::{
        stats::{write_family_statistics, write_inversion_statistics},
        write_index_file, ICON_SVG,
    },
};

#[cfg(not(target_env = "msvc"))]
use tikv_jemallocator::Jemalloc;

#[cfg(not(target_env = "msvc"))]
#[global_allocator]
static GLOBAL: Jemalloc = Jemalloc;

#[derive(Debug, Parser, Clone)]
#[command(name = "aurora")]
#[command(about = "stuff")]
pub struct AuroraArgs {
    /// The path to CAF formatted alignments
    #[arg()]
    alignments: String,

    /// The path to substitution matrices
    #[arg()]
    matrices: String,

    #[command(flatten)]
    #[clap(next_help_heading = "Annotation options")]
    pub annotation_args: AnnotationArgs,

    #[command(flatten)]
    #[clap(next_help_heading = "Performance options")]
    pub performance_args: PerformanceArgs,

    #[command(flatten)]
    #[clap(next_help_heading = "File I/O options")]
    pub io_args: IoArgs,

    #[command(flatten)]
    #[clap(next_help_heading = "Ultra options")]
    pub ultra_args: UltraArgs,

    #[command(flatten)]
    #[clap(next_help_heading = "Visualization options")]
    pub visualization_args: VisualizationArgs,
}

#[derive(Args, Debug, Clone, Default)]
pub struct PerformanceArgs {
    /// Number of threads to use durring processing.
    #[arg(
        short = 't',
        long = "threads",
        default_value_t = 1usize,
        value_name = "n"
    )]
    pub num_threads: usize,
}

#[derive(Args, Debug, Clone, Default)]
pub struct AnnotationArgs {
    /// The penalty of jumping between query models
    #[arg(
        short = 'J',
        long = "query-jump",
        default_value = "-127.0",
        value_name = "f"
    )]
    pub query_jump_penalty: f64,

    /// The number of skip loops that are
    /// equal to a jump between query models
    #[arg(
        short = 'L',
        long = "skip-loop",
        default_value = "30",
        value_name = "n"
    )]
    pub num_skip_loops_eq_to_jump: usize,

    /// The max distance across positions
    /// in the target (genome) at which a join is
    /// considered between compatible alignments
    #[arg(
        short = 'T',
        long = "target-join-distance",
        default_value = "10000",
        value_name = "n"
    )]
    pub target_join_distance: usize,

    /// Removes joins across positions
    /// in the target (genome) at which a join is
    /// less than this likely to not be generated
    /// at random.
    #[arg(
        long = "target-join-likelihood-threshold",
        default_value = "0.5",
        value_name = "f"
    )]
    pub target_distance_likelihood_threshold: f64,

    /// The maximum overlap in the consensus at which
    /// a join is considered between compatible alignments.
    #[arg(
        short = 'O',
        long = "consensus-join-overlap",
        default_value = "200",
        value_name = "n"
    )]
    pub consensus_join_overlap: isize,

    /// The maximum consensus seperation distance at which
    /// a join is considered between compatible alignments.
    #[arg(
        short = 'C',
        long = "consensus-join-distance",
        default_value = "3750",
        value_name = "n"
    )]
    pub consensus_join_distance: isize,

    /// The maximum seperation or overlap in nucleotides on both target and consensus
    /// for a join to be allowed between inverted alignments.
    #[arg(long = "inversion-distance", default_value = "50", value_name = "n")]
    pub inversion_distance: isize,

    /// The size of the window looked at to determine a single alignment score in nucleotides.
    #[arg(
        short = 'W',
        long = "window-size",
        default_value = "31",
        value_name = "n"
    )]
    pub score_window_size: usize,

    /// The size of the window looked at to determine alignment score of the background refrence score in nucleotides.
    #[arg(
        short = 'B',
        long = "background-window-size",
        default_value = "61",
        value_name = "n"
    )]
    pub background_window_size: usize,

    /// Apply an additional penalty to the skip state score.
    #[arg(
        short = 'S',
        long = "skip-state-penalty",
        default_value = "0.0",
        value_name = "f"
    )]
    pub skip_state_score_shift: f64,

    /// The minimum cost for keeping an alignment in a segment for history tracing.
    #[arg(
        long = "min-segment-confidence",
        default_value = "0.1",
        value_name = "f"
    )]
    pub min_block_confidence: f64,

    /// The max depth of the histories used for identifying joins.
    #[arg(long = "max-history-depth", default_value = "64", value_name = "n")]
    pub max_history_depth: usize,

    /// The total number of histories allowed in a single segment.
    /// Additional histories are removed, lowest scoring first.
    /// Set to 0 to disable.
    #[arg(long = "max-histories", default_value = "10000", value_name = "n")]
    pub max_histories_per_segment: usize,

    /// The lowest score a history can have before being pruned.
    /// This is relative to the best scoring history for a segment.
    /// Set to 0 or greater to disable.
    #[arg(long = "min-history-score", default_value = "-500.0", value_name = "f")]
    pub min_relative_history_score: f64,

    /// The amount of overlap between two joinable sequences in the consensus
    /// before a penalty starts being applied to the join.
    #[arg(long = "free-join-overlap", default_value = "4", value_name = "n")]
    pub free_join_consensus_overlap: usize,

    /// The amount of gap between two joinable sequences
    /// before a penalty starts being applied to the join.
    #[arg(long = "free-join-gap", default_value = "10", value_name = "n")]
    pub free_join_consensus_gap: usize,

    /// The amount of penalty to apply to a join at the maximum allowed consensus overlap
    /// A value of 1 means to apply a penalty equal to a query jump.
    /// The cost grows linearly to this value as the overlap increases.
    #[arg(
        long = "consensus-overlap-penalty",
        default_value = "1.0",
        value_name = "f"
    )]
    pub join_consensus_overlap_penalty: f64,

    /// The amount of penalty to apply to a join at the maximum allowed consensus gap
    /// A value of 1 means to apply a penalty equal to a query jump.
    /// The cost grows linearly to this value as the gap increases.
    #[arg(
        long = "consensus-gap-penalty",
        default_value = "0.5",
        value_name = "f"
    )]
    pub join_consensus_gap_penalty: f64,
}

#[derive(Args, Debug, Clone, Default)]
pub struct IoArgs {
    /// Specify path to save aurora annotations to.
    /// Defaults to sending results to standard output.
    #[arg(short = 'o', long = "output", value_name = "path")]
    pub output_path: Option<PathBuf>,
    /// Specify path to dump verbose annotations (with all ambiguous annotation options) to.
    /// Defaults to not saving verbose annotations.
    #[arg(short = 'a', long = "ambiguity-file", value_name = "path")]
    pub ambiguity_path: Option<PathBuf>,
    /// Produce a file that describes the regions
    #[arg(long = "regions", value_name = "path")]
    pub regions_path: Option<PathBuf>,
}

#[derive(Args, Debug, Clone)]
pub struct UltraArgs {
    /// The path to ULTRA output. Must be ULTRA's JSON output format.
    #[arg(short = 'U', long = "ultra-file", value_name = "path")]
    pub ultra_file_path: Option<PathBuf>,

    /// Don't adjudicate regions that are only made up of tandem repeats
    #[arg(short = 'X', long = "exclude-isolated-tr")]
    pub exclude_isolated_tandem_repeats: bool,
}

#[derive(Args, Debug, Clone, Default)]
pub struct VisualizationArgs {
    /// Produce visualization output for annotations
    #[arg(short = 'V', long = "viz")]
    pub viz: bool,

    /// The path to the directory to which
    /// visualization output will be written
    #[arg(long = "viz-out", default_value = "./viz", value_name = "path")]
    pub viz_output_path: PathBuf,

    /// A list of target names, starts, and ends
    /// that will constrain the visualization output
    #[arg(
        long = "viz-constraint",
        value_name = "\"target_name:start:end, ...\"",
        value_delimiter = ','
    )]
    pub viz_constraints: Vec<VizConstraint>,

    /// Enable output of per position scores to the visual.
    #[arg(long = "viz-enable-scores")]
    pub viz_enable_scores: bool,

    /// The path to the BED file that contains
    /// reference annotations for visualization
    #[arg(short = 'R', long = "viz-ref-bed", value_name = "path")]
    pub viz_reference_bed_path: Option<PathBuf>,

    /// Dump additional debug files to the visualization.
    #[arg(long = "debug")]
    pub debug: bool,

    /// Disable history tracing entirely, dumping only visuals.
    #[arg(long = "disable-tracing")]
    pub disable_tracing: bool,

    #[clap(skip)]
    pub viz_reference_bed_index: HashMap<String, usize>,
}

fn main() -> Result<()> {
    let mut args = AuroraArgs::parse();
    let viz_args = &mut args.visualization_args;

    if viz_args.viz {
        if let Result::Ok(metadata) = fs::metadata(&viz_args.viz_output_path) {
            if metadata.is_dir() {
                // TODO: real error
                panic!(
                    "directory: {} already exists",
                    viz_args.viz_output_path.to_str().unwrap()
                )
            }
        }

        create_dir_all(&viz_args.viz_output_path)?;
        viz_args.viz_output_path = viz_args.viz_output_path.canonicalize()?;

        if let Some(path) = &viz_args.viz_reference_bed_path {
            let file = File::open(path).expect("failed to open viz reference bed file");
            let reader = BufReader::new(file);

            let mut chrom_list = vec![String::from("sentinel")];
            let mut prev_start = 0usize;
            let mut index: HashMap<String, usize> = HashMap::new();
            reader
                .lines()
                .map(|l| l.unwrap())
                .enumerate()
                .for_each(|(line_num, line)| {
                    let tokens: Vec<&str> = line.split_whitespace().collect();
                    let chrom = tokens[0].to_string();
                    let start = tokens[1].parse::<usize>().expect("failed to parse int");

                    let last_chrom = chrom_list.last().expect("chrom list is empty");

                    if chrom == *last_chrom {
                        if prev_start > start {
                            panic!("bed file is unsorted");
                        }
                    } else if !chrom_list.contains(&chrom) {
                        chrom_list.push(chrom.clone());
                        index.insert(chrom, line_num);
                    } else {
                        panic!("bed file is unsorted");
                    }

                    prev_start = start;
                });

            viz_args.viz_reference_bed_index = index;
        }
    }

    let alignments_file = File::open(&args.alignments)?;
    let matrices_file = File::open(&args.matrices)?;

    let ultra_file = match args.ultra_args.ultra_file_path {
        Some(ref path) => Some(File::open(path)?),
        None => None,
    };

    let alignment_data =
        AlignmentData::from_caf_and_ultra_and_matrices(alignments_file, ultra_file, matrices_file)?;

    let proximity_groups = ProximityGroup::from_alignment_data(
        &alignment_data,
        args.annotation_args.target_join_distance,
    )
    .into_iter()
    .filter(|g| {
        if args.ultra_args.exclude_isolated_tandem_repeats {
            !g.alignments.is_empty()
        } else {
            true
        }
    })
    .collect_vec();

    if let Some(path) = &args.io_args.regions_path {
        let regions_file = File::create(path).unwrap();
        let mut regions_writer = BufWriter::new(regions_file);
        proximity_groups.iter().enumerate().for_each(|(idx, g)| {
            writeln!(
                &mut regions_writer,
                "{},{},{}:{},{}:{}",
                idx,
                alignment_data.target_name_map.get(g.target_id),
                g.target_start,
                g.target_end,
                g.line_start,
                g.line_end,
            )
            .expect("failed to write to regions file")
        });
    }

    if viz_args.viz {
        let mut index_file = File::create(viz_args.viz_output_path.join("index.html")).unwrap();

        write_index_file(
            &mut index_file,
            &alignment_data,
            &proximity_groups,
            &viz_args.viz_constraints,
        )
        .expect("failed to write to index.html");
    }

    debug_assert!(validate_groups(
        &proximity_groups,
        args.annotation_args.target_join_distance
    ));

    let mut output_file = if let Some(path) = &args.io_args.output_path {
        Some(File::create(path)?)
    } else {
        None
    };

    let mut ambiguity_file = if let Some(path) = &args.io_args.ambiguity_path {
        Some(File::create(path)?)
    } else {
        None
    };

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.performance_args.num_threads)
        .build_global()
        .unwrap();

    let mut naive_results = proximity_groups
        .par_iter()
        .panic_fuse()
        .enumerate()
        .map(|(region_idx, group)| run_naive_trace(group, &alignment_data, region_idx, &args))
        .collect::<Vec<NaiveTraceResults>>();
    naive_results.sort_by_key(|v| v.region_index);

    let trace_stats = trace_statistics(
        &naive_results,
        &alignment_data,
        OccuranceCountingMode::Segments,
    );

    let mut results: Vec<(usize, Vec<AmbiguousAnnotation>)> = proximity_groups
        .par_iter()
        .zip(naive_results)
        .map(|(group, mut naive_trace)| {
            (
                naive_trace.region_index,
                run_history_trace(
                    group,
                    &alignment_data,
                    &trace_stats,
                    &mut naive_trace,
                    &args,
                ),
            )
        })
        .collect();
    results.sort_by_key(|v| v.0);

    for (_region, annots) in results.iter() {
        if let Some(file_out) = output_file.as_mut() {
            AmbiguousAnnotation::write(annots, file_out, true)?;
        } else {
            AmbiguousAnnotation::write(annots, &mut std::io::stdout(), true)?;
        }

        if let Some(amb_file_out) = ambiguity_file.as_mut() {
            AmbiguousAnnotation::write(annots, amb_file_out, false)?;
        }
    }

    if args.visualization_args.viz {
        let mut family_stats_writer = File::create(
            args.visualization_args
                .viz_output_path
                .join("family_stats.html"),
        )?;
        write_family_statistics(&mut family_stats_writer, &results)?;
        let mut inv_stats_writer = File::create(
            args.visualization_args
                .viz_output_path
                .join("inversion_stats.html"),
        )?;
        write_inversion_statistics(&mut inv_stats_writer, &results)?;
        let mut icon_file = File::create(args.visualization_args.viz_output_path.join("icon.svg"))?;
        icon_file.write_all(ICON_SVG.as_bytes())?;
    }

    Ok(())
}
