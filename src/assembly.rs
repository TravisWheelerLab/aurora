use std::{collections::HashMap, hash::Hash};

use itertools::Itertools;

use crate::{
    alignment::{Alignment, Strand},
    score_params::ScoreParams,
    segments::SegmentedMatrix,
    AnnotationArgs,
};

#[derive(Hash, Eq, PartialEq, Clone, Copy, Debug)]
pub enum LinkType {
    Forward,
    Reverse,
    FRInversion,
    RFInversion,
}

#[derive(Clone, Copy, Debug)]
pub struct Edge {
    pub weight: f64,
    pub first_sparse_row: usize,
    #[allow(dead_code)]
    pub second_sparse_row: usize,
    #[allow(dead_code)]
    pub link_type: LinkType,
}

fn piecewise_linear_cost(
    neg_start: f64,
    pos_start: f64,
    neg_slope: f64,
    pos_slope: f64,
    value: f64,
) -> f64 {
    if value < neg_start {
        (value - neg_start).abs() * neg_slope
    } else if value > pos_slope {
        (value - pos_start).abs() * pos_slope
    } else {
        0.0
    }
}

fn get_link_cost(
    annotation_args: &AnnotationArgs,
    score_params: &ScoreParams,
    consensus_gap: f64,
    target_gap: f64,
) -> f64 {
    // Minimum cost (a query loop)
    let min_value = score_params.query_loop_score;
    let value_range = (score_params.query_loop_score - score_params.query_jump_score).abs();

    // Get overlap and gap ranges with free areas incorperated in, otherwise math is not quite right.
    let overlap_range = ((annotation_args.consensus_join_overlap as f64)
        - (annotation_args.free_join_consensus_overlap as f64))
        .abs()
        .max(1.0);
    let gap_range = ((annotation_args.consensus_join_distance as f64)
        - (annotation_args.free_join_consensus_gap as f64))
        .abs()
        .max(1.0);

    // Compute slopes....
    let lambda = -value_range
        * (annotation_args.join_target_gap_penalty
            / annotation_args.target_join_distance.max(1) as f64)
            .abs();
    let alpha =
        -value_range * (annotation_args.join_consensus_overlap_penalty / overlap_range).abs();
    let beta = -value_range * (annotation_args.join_consensus_gap_penalty / gap_range).abs();

    // Cost = linear consensus cost + linear target gap cost...
    min_value
        + piecewise_linear_cost(
            -(annotation_args.free_join_consensus_overlap as f64).abs(),
            (annotation_args.free_join_consensus_gap as f64).abs(),
            alpha,
            beta,
            consensus_gap,
        )
        + lambda * target_gap
}

fn link_assemblies(
    graph: &mut HashMap<(SegmentAndDenseRow, SegmentAndDenseRow), Edge>,
    compatable_blocks: impl Iterator<Item = (usize, usize)>,
    alignments: &[Alignment],
    segments: &SegmentedMatrix,
    score_params: &ScoreParams,
    args: &AnnotationArgs,
) {
    // this relies on the alignments being sorted by target start
    // note: this assertion iter will only run in debug mode
    let compatable_blocks = compatable_blocks.sorted().collect_vec();

    compatable_blocks.iter().enumerate().for_each(|(idx, a)| {
        compatable_blocks[idx + 1..].iter().for_each(|b| {
            // Same segment, don't allow merging...
            if a.0 == b.0 {
                return;
            }

            let a_block = &segments[a.0].blocks[a.1];
            let b_block = &segments[b.0].blocks[b.1];

            // If same alignment, and neighboring segments, don't join...
            if a_block.row_idx == b_block.row_idx && ((b.0 - 1) <= a.0) {
                return;
            }

            let target_distance = b_block.col_start as isize - a_block.col_end as isize;

            let a_length = a_block.query_end.abs_diff(a_block.query_start);
            let b_length = b_block.query_end.abs_diff(b_block.query_start);
            let min_length = a_length.min(b_length);

            // Query bounds are reversed for reverse sequences, so the start is actually greater than the end (Ex. start: 1510 -> end: 105)

            let (consensus_distance, link_type) = match (
                alignments[a_block.row_idx - 1].strand,
                alignments[b_block.row_idx - 1].strand,
            ) {
                (Strand::Forward, Strand::Forward) => (
                    b_block.query_start as isize - a_block.query_end as isize,
                    LinkType::Forward,
                ),
                (Strand::Reverse, Strand::Reverse) => (
                    a_block.query_end as isize - b_block.query_start as isize,
                    LinkType::Reverse,
                ),
                (Strand::Forward, Strand::Reverse) => (
                    b_block.query_end as isize - a_block.query_end as isize,
                    LinkType::FRInversion,
                ),
                (Strand::Reverse, Strand::Forward) => (
                    a_block.query_end as isize - b_block.query_end as isize,
                    LinkType::RFInversion,
                ),
                _ => panic!("Invalid strand types!"),
            };

            let within_target_distance_threshold = match link_type {
                LinkType::FRInversion | LinkType::RFInversion => {
                    target_distance.abs() < args.inversion_distance
                }
                _ => target_distance < args.target_join_distance as isize,
            };

            let consensus_is_colinear = match link_type {
                LinkType::FRInversion | LinkType::RFInversion => {
                    consensus_distance.abs() < args.inversion_distance
                }
                _ => {
                    consensus_distance > -args.consensus_join_overlap
                        && consensus_distance < args.consensus_join_distance
                }
            };

            // TODO: Hardcoded, change later...
            let is_significant =
                min_length >= 10 && -consensus_distance <= ((min_length / 2) as isize);

            let weight = get_link_cost(
                args,
                score_params,
                consensus_distance as f64,
                target_distance as f64,
            );

            // let not_reached_forward_limit = forward_count < args.max_forward_links;

            if within_target_distance_threshold && consensus_is_colinear && is_significant {
                graph.insert(
                    ((a.0, a_block.row_idx), (b.0, b_block.row_idx)),
                    Edge {
                        weight,
                        first_sparse_row: a.1,
                        second_sparse_row: b.1,
                        link_type,
                    },
                );
            }
        });
    });
}

type SegmentAndDenseRow = (usize, usize);

/// Represents graph of compatable alignments on the genome.
/// For each alignment, stores all alignments from the same query in front of it.
pub struct SegmentAssemblyGraph {
    #[allow(dead_code)]
    pub alignment_block_map: Vec<Vec<(usize, usize)>>, // Maps alignment to it's corresponding blocks...
    pub link_graph: HashMap<(SegmentAndDenseRow, SegmentAndDenseRow), Edge>,
}

impl SegmentAssemblyGraph {
    pub fn new(
        alignments: &[Alignment],
        segments: &SegmentedMatrix,
        score_params: &ScoreParams,
        annotation_args: &AnnotationArgs,
    ) -> Self {
        let mut alignment_block_map = vec![Vec::<SegmentAndDenseRow>::new(); alignments.len()];

        for (s_idx, segment) in segments.iter().enumerate() {
            for (b_idx, block) in segment.blocks.iter().enumerate() {
                if block.row_idx > 0 && block.row_idx <= alignments.len() {
                    alignment_block_map[block.row_idx - 1].push((s_idx, b_idx));
                }
            }
        }

        let mut query_ids: Vec<usize> = alignments.iter().map(|a| a.query_id).unique().collect();

        query_ids.sort();

        let mut link_graph = HashMap::new();

        query_ids
            .iter()
            // grab the alignments for this ID
            .map(|id| {
                alignments
                    .iter()
                    .enumerate()
                    .filter(|&(_, a)| a.query_id == *id)
                    .flat_map(|(a_idx, _)| alignment_block_map[a_idx].iter().copied())
            })
            .for_each(|compat_blocks| {
                link_assemblies(
                    &mut link_graph,
                    compat_blocks,
                    alignments,
                    segments,
                    score_params,
                    annotation_args,
                );
            });

        // Graph is constructed such that first block is always before the second block, this debug assert checks that...
        link_graph.keys().for_each(|(a, b)| debug_assert!(a < b));

        Self {
            alignment_block_map,
            link_graph,
        }
    }
}
