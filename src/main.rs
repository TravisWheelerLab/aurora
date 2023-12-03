mod alignment;
mod alphabet;
mod chunks;
mod collapse;
mod confidence;
mod matrix;
mod pipeline;
mod results;
mod score_params;
mod segments;
mod split;
mod substitution_matrix;
mod support;
mod viterbi;
mod viz;
mod windowed_scores;

use std::{
    collections::HashMap,
    fs::{create_dir_all, File},
    io::{BufRead, BufReader},
    path::PathBuf,
};

use alignment::AlignmentData;
use chunks::ProximityGroup;

use anyhow::Result;
use clap::Parser;

use crate::{chunks::validate_groups, pipeline::run_assembly_pipeline};

#[derive(Debug, Parser, Clone)]
#[command(name = "aurora")]
#[command(about = "stuff")]
pub struct Args {
    /// The path to CAF formatted alignments
    #[arg()]
    alignments: String,

    /// The path to substitution matrices
    #[arg()]
    matrices: String,

    /// The probability of jumping between query models
    #[arg(
        short = 'J',
        long = "query-jump",
        default_value = "1e-55",
        value_name = "f"
    )]
    pub query_jump_probability: f64,

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

    ///
    #[arg(
        short = 'F',
        long = "fudge-distance",
        default_value = "10",
        value_name = "n"
    )]
    pub fudge_distance: usize,

    /// The path to ULTRA output
    #[arg(short = 'U', long = "ultra-file", value_name = "path")]
    pub ultra_file_path: Option<PathBuf>,

    /// Produce visualization output
    #[arg(long = "viz")]
    pub viz: bool,

    /// Produce visualization output
    #[arg(long = "assembly-viz")]
    pub assembly_viz: bool,

    /// The path to the directory to which
    /// vizualization output will be written
    #[arg(long = "viz-out", default_value = "./viz", value_name = "path")]
    pub viz_output_path: PathBuf,

    /// The path to the BED file that contains
    /// reference annotations for visualization
    #[arg(long = "viz-ref-bed", value_name = "path")]
    pub viz_reference_bed_path: Option<PathBuf>,

    #[clap(skip)]
    pub viz_reference_bed_index: HashMap<String, usize>,
}

pub const SCORE_WINDOW_SIZE: usize = 31;
pub const BACKGROUND_WINDOW_SIZE: usize = 61;

fn main() -> Result<()> {
    let mut args = Args::parse();

    if args.viz {
        create_dir_all(&args.viz_output_path)?;
        args.viz_output_path = args.viz_output_path.canonicalize()?;

        if let Some(path) = &args.viz_reference_bed_path {
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

            args.viz_reference_bed_index = index;
        }
    }

    let alignments_file = File::open(&args.alignments)?;
    let matrices_file = File::open(&args.matrices)?;

    let ultra_file = match args.ultra_file_path {
        Some(ref path) => Some(File::open(path)?),
        None => None,
    };

    let alignment_data =
        AlignmentData::from_caf_and_ultra_and_matrices(alignments_file, ultra_file, matrices_file)?;

    let proximity_groups =
        ProximityGroup::from_alignment_data(&alignment_data, args.target_join_distance);

    debug_assert!(validate_groups(
        &proximity_groups,
        args.target_join_distance
    ));

    proximity_groups
        .iter()
        .inspect(|g| println!("{g:?}"))
        .enumerate()
        .for_each(|(region_idx, group)| {
            run_assembly_pipeline(group, &alignment_data, region_idx, args.clone());
        });
    Ok(())
}
