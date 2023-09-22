mod alignment;
mod alphabet;
mod chunks;
mod confidence;
mod matrix;
mod pipeline;
mod results;
mod score_params;
mod segments;
mod substitution_matrix;
mod support;
mod viterbi;
mod viz;
mod windowed_scores;

use std::fs::File;

use alignment::Alignment;
use chunks::chunks_by_query_distance;
use pipeline::run_pipeline;

use anyhow::Result;
use clap::Parser;
use substitution_matrix::SubstitutionMatrix;

#[derive(Debug, Parser)]
#[command(name = "aurora")]
#[command(about = "stuff")]
pub struct Cli {
    /// The path to CAF formatted alignments
    #[arg()]
    alignments: String,
    /// The path to substitution matrices
    #[arg()]
    matrices: String,
}

pub const DEFAULT_QUERY_JUMP_PROBABILITY: f64 = 1e-55;
pub const NUM_SKIP_LOOPS_EQ_TO_JUMP: usize = 30;

pub const MIN_FRAGMENT_LENGTH: usize = 10;
pub const TARGET_JOIN_DISTANCE: usize = 10_000;
pub const CONSENSUS_JOIN_DISTANCE: usize = 50;

pub const SCORE_WINDOW_SIZE: usize = 31;
pub const BACKGROUND_WINDOW_SIZE: usize = 61;

fn main() -> Result<()> {
    let cli = Cli::parse();
    let alignments_file = File::open(cli.alignments)?;
    let matrices_file = File::open(cli.matrices)?;

    let (mut alignments, query_names) = Alignment::from_caf(alignments_file)?;
    let substitution_matrices = SubstitutionMatrix::parse(matrices_file);

    let num_queries = query_names.len();

    let chunks = chunks_by_query_distance(&mut alignments, num_queries, TARGET_JOIN_DISTANCE);

    chunks
        .iter()
        // .inspect(|c| println!("{} - {}", c.start_idx, c.end_idx))
        .map(|c| &alignments[c.start_idx..=c.end_idx])
        .enumerate()
        .for_each(|(region_idx, ali_slice)| {
            run_pipeline(ali_slice, &query_names, &substitution_matrices, region_idx)
        });

    Ok(())
}
