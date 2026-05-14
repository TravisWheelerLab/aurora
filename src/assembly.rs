use std::{collections::HashMap, hash::Hash};

use itertools::Itertools;

use crate::{
    alignment::{Alignment, Strand},
    join_estimation::{JoinEstimator, JoinStatisticsCollector},
    score_params::ScoreParams,
    segments::{Block, SegmentedMatrix, SegmentedMatrixView},
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
    consensus_gap: isize,
    join_prob: f64,
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
    let alpha =
        -value_range * (annotation_args.join_consensus_overlap_penalty / overlap_range).abs();
    let beta = -value_range * (annotation_args.join_consensus_gap_penalty / gap_range).abs();

    // Doing this as the expected value over the transition scores...
    let expected_score = join_prob * score_params.query_loop_score
        + (1.0 - join_prob) * score_params.query_jump_score;

    // Cost = linear consensus cost + linear target gap cost...
    min_value
        + piecewise_linear_cost(
            -(annotation_args.free_join_consensus_overlap as f64).abs(),
            (annotation_args.free_join_consensus_gap as f64).abs(),
            alpha,
            beta,
            consensus_gap as f64,
        )
        + expected_score
}

pub fn block_target_distance(first_block: &Block, second_block: &Block) -> isize {
    second_block.col_start as isize - first_block.col_end as isize - 1
}

pub fn block_consensus_distance(first_block: &Block, second_block: &Block) -> (isize, LinkType) {
    let select_closest = |prop1: (isize, LinkType), prop2: (isize, LinkType)| {
        if prop1.0.abs() < prop2.0.abs() {
            prop1
        } else {
            prop2
        }
    };

    match (first_block.strand, second_block.strand) {
        (Strand::Forward, Strand::Forward) => (
            second_block.query_start as isize - first_block.query_end as isize - 1,
            LinkType::Forward,
        ),
        (Strand::Reverse, Strand::Reverse) => (
            first_block.query_end as isize - second_block.query_start as isize - 1,
            LinkType::Reverse,
        ),
        (Strand::Forward, Strand::Reverse) => select_closest(
            (
                first_block.query_start as isize - second_block.query_start as isize - 1,
                LinkType::FRInversion1,
            ),
            (
                second_block.query_end as isize - first_block.query_end as isize - 1,
                LinkType::FRInversion2,
            ),
        ),
        (Strand::Reverse, Strand::Forward) => select_closest(
            (
                second_block.query_start as isize - first_block.query_start as isize - 1,
                LinkType::RFInversion1,
            ),
            (
                first_block.query_end as isize - second_block.query_end as isize - 1,
                LinkType::RFInversion2,
            ),
        ),
        _ => panic!("Invalid strand types!"),
    }
}

pub fn block_length_on_query(b: &Block) -> usize {
    b.query_end.abs_diff(b.query_start) + 1
}

fn is_joinable(
    target_distance: isize,
    consensus_distance: isize,
    link_type: LinkType,
    min_block_length: usize,
    args: &AnnotationArgs,
) -> bool {
    let within_target_distance_threshold =
        target_distance < args.target_join_distance as isize && target_distance >= 0;

    let consensus_is_colinear = if link_type.is_inversion() {
        consensus_distance.abs() < args.inversion_distance
    } else {
        consensus_distance > -args.consensus_join_overlap
            && consensus_distance < args.consensus_join_distance
    };

    // TODO: Hardcoded, change later...
    let is_significant =
        min_block_length >= 10 && -consensus_distance <= ((min_block_length / 2) as isize);

    within_target_distance_threshold && consensus_is_colinear && is_significant
}

fn new_alignment_to_blocks_map(
    segments: SegmentedMatrixView,
    alignments: &[Alignment],
) -> Vec<Vec<SegmentAndDenseRow>> {
    let mut alignment_block_map = vec![Vec::<SegmentAndDenseRow>::new(); alignments.len()];

    for (s_idx, segment) in segments.iter().enumerate() {
        for (b_idx, block) in segment.blocks.iter().enumerate() {
            if block.row_idx > 0 && block.row_idx <= alignments.len() {
                alignment_block_map[block.row_idx - 1].push((s_idx, b_idx));
            }
        }
    }

    alignment_block_map
}

pub fn gather_join_statistics<T: JoinStatisticsCollector>(
    alignments: &[Alignment],
    annotation_args: &AnnotationArgs,
) -> Vec<(usize, T)> {
    let mut query_ids: Vec<usize> = alignments.iter().map(|a| a.query_id).unique().collect();
    query_ids.sort();

    let mut query_stats: Vec<(usize, T)> = Vec::with_capacity(query_ids.len());

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
                    .map(|(i, a)| Block::from_alignment(a, i, 0.0, 0.0)),
            )
        })
        .for_each(|(id, compat_alignments)| {
            let mut new_stats = T::new();

            gather_join_statistics_single_family(
                compat_alignments,
                annotation_args,
                &mut new_stats,
            );

            query_stats.push((id, new_stats));
        });

    query_stats
}

fn gather_join_statistics_single_family<'a>(
    compatable_alignments: impl Iterator<Item = Block>,
    args: &AnnotationArgs,
    join_stats: &mut impl JoinStatisticsCollector,
) {
    let compatable_blocks = compatable_alignments
        .sorted_by_key(|a| a.col_start)
        .collect_vec();

    compatable_blocks
        .iter()
        .enumerate()
        .for_each(|(idx, a_block)| {
            compatable_blocks[idx + 1..]
                .iter()
                .enumerate()
                .for_each(|(idx2, b_block)| {
                    let (consensus_distance, link_type) =
                        block_consensus_distance(a_block, b_block);
                    let joinable = is_joinable(
                        block_target_distance(a_block, b_block),
                        consensus_distance,
                        link_type,
                        block_length_on_query(a_block).min(block_length_on_query(b_block)),
                        args,
                    );

                    join_stats.add(a_block, b_block, idx + 1 == idx2, joinable);
                })
        })
}

fn link_assemblies<T: JoinEstimator>(
    graph: &mut HashMap<(SegmentAndDenseRow, SegmentAndDenseRow), Edge>,
    compatable_blocks: impl Iterator<Item = (usize, usize)>,
    segments: &SegmentedMatrix,
    query_statistics: &QueryStatistics<T>,
    _region_statistics: &RegionStatistics,
    score_params: &ScoreParams,
    args: &AnnotationArgs,
) {
    // this relies on the alignments being sorted by target start
    let compatable_blocks = compatable_blocks.sorted().collect_vec();

    compatable_blocks.iter().enumerate().for_each(|(idx, a)| {
        compatable_blocks[idx + 1..].iter().for_each(|b| {
            // Same segment, don't allow merging...
            if a.0 == b.0 {
                return;
            }

            let a_block = &segments[a.0].blocks[a.1];
            let b_block = &segments[b.0].blocks[b.1];

            let target_distance = block_target_distance(a_block, b_block);
            let min_block_length =
                block_length_on_query(a_block).min(block_length_on_query(b_block));

            let (consensus_distance, link_type) = block_consensus_distance(a_block, b_block);

            if b_block.row_idx == 583 {
                println!("Block: {}", a_block.row_idx);
                println!(
                    "Score: {}",
                    query_statistics.estimator.predict(a_block, b_block, false)
                );

                println!(
                    "Is Joinable: {}",
                    is_joinable(
                        target_distance,
                        consensus_distance,
                        link_type,
                        min_block_length,
                        args,
                    )
                );

                println!(
                    "Weight: {}",
                    if a_block.row_idx == b_block.row_idx && ((b.0 - 1) <= a.0) {
                        score_params.query_loop_score
                    } else {
                        get_link_cost(
                            args,
                            score_params,
                            consensus_distance,
                            query_statistics.estimator.predict(a_block, b_block, false),
                        )
                    }
                );

                println!("Estimator: {:#?}", query_statistics.estimator);
                println!(
                    "Target Dist: {}, Div: {}, Cons Dist: {}",
                    target_distance,
                    (a_block.kimura80 - b_block.kimura80).abs(),
                    consensus_distance
                )
            }

            if is_joinable(
                target_distance,
                consensus_distance,
                link_type,
                min_block_length,
                args,
            ) {
                let join_prob = query_statistics.estimator.predict(a_block, b_block, false);

                if join_prob >= args.join_likelihood_threshold {
                    let mut weight = if a_block.row_idx == b_block.row_idx && ((b.0 - 1) <= a.0) {
                        score_params.query_loop_score
                    } else {
                        get_link_cost(args, score_params, consensus_distance, join_prob)
                    };

                    if b_block.query_id == Some(196) {
                        println!("Setting Weight to 1 for {}", b_block.row_idx);
                        weight = score_params.query_loop_score
                    }

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
    pub fn new<T: JoinEstimator>(
        alignments: &[Alignment],
        segments: &SegmentedMatrix,
        region_statistics: &RegionStatistics,
        query_statistics: &[QueryStatistics<T>],
        score_params: &ScoreParams,
        annotation_args: &AnnotationArgs,
    ) -> Self {
        let alignment_block_map = new_alignment_to_blocks_map(segments, alignments);
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
