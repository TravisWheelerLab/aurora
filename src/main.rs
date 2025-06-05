mod alignment;
mod alphabet;
mod annotation;
mod chunks;
mod collapse;
mod confidence;
mod matrix;
mod pipeline;
mod score_params;
mod split;
mod substitution_matrix;
mod support;
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

use anyhow::Result;
use clap::{Args, Parser};
use itertools::Itertools;
use rayon::prelude::*;
use viz::VizConstraint;

use crate::{chunks::validate_groups, pipeline::run_pipeline};

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

    /// The max distance across unaligned positions
    /// in the target (genome) at which a join is
    /// considered between compatible alignments
    #[arg(
        short = 'T',
        long = "target-join-distance",
        default_value = "10000",
        value_name = "n"
    )]
    pub target_join_distance: usize,

    /// The max consensus position difference at which
    /// a join is considered between compatible alignments.
    #[arg(
        short = 'C',
        long = "consensus-join-distance",
        default_value = "50",
        value_name = "n"
    )]
    pub consensus_join_distance: usize,

    /// The minimum length of an alignment fragment
    /// at which a join is considered between another
    /// alignment fragment.
    #[arg(
        short = 'M',
        long = "min-fragment-length",
        default_value = "10",
        value_name = "n"
    )]
    pub min_fragment_length: usize,

    /// The distance used to approximate various
    /// alignment overlap conditions.
    #[arg(
        short = 'F',
        long = "fudge-distance",
        default_value = "10",
        value_name = "n"
    )]
    pub fudge_distance: usize,

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
}

#[derive(Args, Debug, Clone, Default)]
pub struct IoArgs {
    /// Produce a file that describes the regions
    #[arg(long = "regions", value_name = "path")]
    pub regions_path: Option<PathBuf>,
}

#[derive(Args, Debug, Clone)]
pub struct UltraArgs {
    /// The path to ULTRA output
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

    /// Produce visualization output for potential join "assemblies"
    #[arg(long = "assembly-viz")]
    pub assembly_viz: bool,

    /// A list of target names, starts, and ends
    /// that will constrain the visualization output
    #[arg(
        long = "viz-constraint",
        value_name = "\"target_name:start:end, ...\"",
        value_delimiter = ','
    )]
    pub viz_constraints: Vec<VizConstraint>,

    /// The path to the BED file that contains
    /// reference annotations for visualization
    #[arg(short = 'R', long = "viz-ref-bed", value_name = "path")]
    pub viz_reference_bed_path: Option<PathBuf>,

    #[clap(skip)]
    pub viz_reference_bed_index: HashMap<String, usize>,
}

fn main() -> Result<()> {
    let mut args = AuroraArgs::parse();
    let vis_args = &mut args.visualization_args;

    if vis_args.viz {
        if let Ok(metadata) = fs::metadata(&vis_args.viz_output_path) {
            if metadata.is_dir() {
                // TODO: real error
                panic!(
                    "directory: {} already exists",
                    vis_args.viz_output_path.to_str().unwrap()
                )
            }
        }

        create_dir_all(&vis_args.viz_output_path)?;
        vis_args.viz_output_path = vis_args.viz_output_path.canonicalize()?;

        if let Some(path) = &vis_args.viz_reference_bed_path {
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

            vis_args.viz_reference_bed_index = index;
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

    if vis_args.viz || vis_args.assembly_viz {
        let error_msg = "failed to write to index.html";
        let index_file = File::create(vis_args.viz_output_path.join("index.html")).unwrap();
        let mut index_writer = BufWriter::new(index_file);

        if vis_args.viz {
            vis_args
                .viz_constraints
                .iter()
                .enumerate()
                .for_each(|(idx, c)| {
                    writeln!(
                        &mut index_writer,
                        "<a href={}-{}-{}.html>slice {} | {} {}:{}</a><br>",
                        c.target_name,
                        c.target_start,
                        c.target_end,
                        idx,
                        c.target_name,
                        c.target_start,
                        c.target_end,
                    )
                    .expect(error_msg);
                });
        }

        proximity_groups.iter().enumerate().for_each(|(idx, g)| {
            writeln!(
                &mut index_writer,
                "<h3>region {} | {} {}:{}</h3>\n<ul>",
                idx,
                alignment_data.target_name_map.get(g.target_id),
                g.target_start,
                g.target_end,
            )
            .expect(error_msg);

            if vis_args.viz {
                writeln!(
                    &mut index_writer,
                    "    <li><a href={}/index.html>Annotations</a></li>",
                    idx,
                )
                .expect(error_msg);
            }

            if vis_args.assembly_viz {
                writeln!(
                    &mut index_writer,
                    "    <li><a href={}/assembly_index.html>Assemblies</a></li>",
                    idx,
                )
                .expect(error_msg);
            }

            writeln!(&mut index_writer, "</ul>").expect(error_msg);
        });
    }

    debug_assert!(validate_groups(
        &proximity_groups,
        args.annotation_args.target_join_distance
    ));

    rayon::ThreadPoolBuilder::new()
        .num_threads(args.performance_args.num_threads)
        .build_global()
        .unwrap();

    proximity_groups
        .par_iter()
        .panic_fuse()
        // .inspect(|g| println!("{g:?}"))
        .enumerate()
        .for_each(|(region_idx, group)| {
            run_pipeline(group, &alignment_data, region_idx, args.clone());
        });

    Ok(())
}
