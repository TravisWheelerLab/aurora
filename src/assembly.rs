use std::{collections::HashMap, hash::Hash};

use itertools::Itertools;

use crate::{
    alignment::{Alignment, Strand},
    score_params::ScoreParams,
    segments::SegmentedMatrix,
    statistics::Distribution,
    trace_statistics::{QueryStatistics, RegionStatistics},
    AnnotationArgs,
};

#[derive(Hash, Eq, PartialEq, Clone, Copy, Debug)]
pub enum LinkType {
    Forward,
    Reverse,
    // Forward to reverse strand inversions...
    FRInversion1, // 1st sequence flipped.
    FRInversion2, // 2nd sequence flipped.
    // Reverse to forward strand inversions...
    RFInversion1, // 1st sequence flipped...
    RFInversion2, // 2nd sequence flipped...
}

/// The side of the sequence being referred to, in target (genome), space...
#[derive(Debug, Clone, Ord, PartialEq, PartialOrd, Eq, Copy)]
pub enum Side {
    Left,
    Right,
}

impl Side {
    pub fn flip(&self) -> Self {
        match self {
            Self::Left => Self::Right,
            Self::Right => Self::Left,
        }
    }

    pub fn to_index(&self) -> usize {
        match self {
            Self::Left => 0,
            Self::Right => 1,
        }
    }
}

impl LinkType {
    pub fn is_inversion(&self) -> bool {
        matches!(
            self,
            Self::FRInversion1 | Self::FRInversion2 | Self::RFInversion1 | Self::RFInversion2
        )
    }

    /// Get the linked sides of two segments. The first one is from the first sequence the genome, the second from the second one.
    pub fn get_linked_sides(&self) -> (Side, Side) {
        match self {
            Self::Forward | Self::Reverse => (Side::Right, Side::Left),
            Self::FRInversion1 | Self::RFInversion1 => (Side::Left, Side::Left),
            Self::FRInversion2 | Self::RFInversion2 => (Side::Right, Side::Right),
        }
    }

    #[allow(dead_code)]
    /// Get the unlinked, or still open sides of two segments. The first one is from the first sequence the genome, the second from the second one.
    pub fn get_open_sides(&self) -> (Side, Side) {
        let linked = self.get_linked_sides();
        (linked.0.flip(), linked.1.flip())
    }
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

fn link_assemblies<T: Distribution>(
    graph: &mut HashMap<(SegmentAndDenseRow, SegmentAndDenseRow), Edge>,
    compatable_blocks: impl Iterator<Item = (usize, usize)>,
    alignments: &[Alignment],
    segments: &SegmentedMatrix,
    query_statistics: &QueryStatistics<T>,
    region_statistics: &RegionStatistics,
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

            // We allow this now, otherwise inversions might not properly join...
            // If same alignment, and neighboring segments, don't join...
            //if a_block.row_idx == b_block.row_idx && ((b.0 - 1) <= a.0) {
            //    return;
            //}

            let target_distance = b_block.col_start as isize - a_block.col_end as isize - 1;

            let a_length = a_block.query_end.abs_diff(a_block.query_start) + 1;
            let b_length = b_block.query_end.abs_diff(b_block.query_start) + 1;
            let min_length = a_length.min(b_length);

            let select_closest = |prop1: (isize, LinkType), prop2: (isize, LinkType)| {
                if prop1.0.abs() < prop2.0.abs() {
                    prop1
                } else {
                    prop2
                }
            };

            // Query bounds are reversed for reverse sequences, so the start is actually greater than the end (Ex. start: 1510 -> end: 105)

            let (consensus_distance, link_type) = match (
                alignments[a_block.row_idx - 1].strand,
                alignments[b_block.row_idx - 1].strand,
            ) {
                (Strand::Forward, Strand::Forward) => (
                    b_block.query_start as isize - a_block.query_end as isize - 1,
                    LinkType::Forward,
                ),
                (Strand::Reverse, Strand::Reverse) => (
                    a_block.query_end as isize - b_block.query_start as isize - 1,
                    LinkType::Reverse,
                ),
                (Strand::Forward, Strand::Reverse) => select_closest(
                    (
                        a_block.query_start as isize - b_block.query_start as isize - 1,
                        LinkType::FRInversion1,
                    ),
                    (
                        b_block.query_end as isize - a_block.query_end as isize - 1,
                        LinkType::FRInversion2,
                    ),
                ),
                (Strand::Reverse, Strand::Forward) => select_closest(
                    (
                        b_block.query_start as isize - a_block.query_start as isize - 1,
                        LinkType::RFInversion1,
                    ),
                    (
                        a_block.query_end as isize - b_block.query_end as isize - 1,
                        LinkType::RFInversion2,
                    ),
                ),
                _ => panic!("Invalid strand types!"),
            };

            let within_target_distance_threshold =
                target_distance < args.target_join_distance as isize;

            let consensus_is_colinear = if link_type.is_inversion() {
                consensus_distance.abs() < args.inversion_distance
            } else {
                consensus_distance > -args.consensus_join_overlap
                    && consensus_distance < args.consensus_join_distance
            };

            // TODO: Hardcoded, change later...
            let is_significant =
                min_length >= 10 && -consensus_distance <= ((min_length / 2) as isize);

            let weight = if a_block.row_idx == b_block.row_idx && ((b.0 - 1) <= a.0) {
                score_params.query_loop_score
            } else {
                get_link_cost(
                    args,
                    score_params,
                    consensus_distance as f64,
                    target_distance as f64,
                )
            };

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
    pub fn new<T: Distribution>(
        alignments: &[Alignment],
        segments: &SegmentedMatrix,
        region_statistics: &RegionStatistics,
        query_statistics: &[QueryStatistics<T>],
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
                (
                    *id,
                    alignments
                        .iter()
                        .enumerate()
                        .filter(|&(_, a)| a.query_id == *id)
                        .flat_map(|(a_idx, _)| alignment_block_map[a_idx].iter().copied()),
                )
            })
            .for_each(|(id, compat_blocks)| {
                link_assemblies(
                    &mut link_graph,
                    compat_blocks,
                    alignments,
                    segments,
                    &query_statistics[id],
                    region_statistics,
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
