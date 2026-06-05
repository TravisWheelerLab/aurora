use std::{collections::HashMap, hash::Hash};

use itertools::Itertools;

use crate::{
    alignment::{Alignment, Strand},
    chunks::ProximityGroup,
    join_estimation::{JoinEstimator, JoinStatisticsCollector, LinkInfo},
    segments::{Block, InitialSegments, SegmentedMatrix, SegmentedMatrixView},
    trace_statistics::{calculate_region_statistics, QueryStatistics, RegionStatistics},
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

#[allow(dead_code)]
pub enum ConsensusDistanceNormalization {
    Max,
    Min,
    Sum,
    WithLength(usize),
    WithUBAndLength(usize, usize),
}

pub fn relative_consensus_distance(
    first_block: &Block,
    second_block: &Block,
    mode: ConsensusDistanceNormalization,
) -> (f64, LinkType) {
    let (mut dist, link_type) = block_consensus_distance(first_block, second_block);

    if let ConsensusDistanceNormalization::WithUBAndLength(ub, _length) = mode {
        dist = if dist > 0 {
            dist.saturating_sub(ub as isize).max(0)
        } else {
            dist
        }
    }

    let div = match mode {
        ConsensusDistanceNormalization::Sum => {
            block_length_on_query(first_block) + block_length_on_query(second_block)
        }
        ConsensusDistanceNormalization::Max => {
            block_length_on_query(first_block).max(block_length_on_query(second_block))
        }
        ConsensusDistanceNormalization::Min => {
            block_length_on_query(first_block).min(block_length_on_query(second_block))
        }
        ConsensusDistanceNormalization::WithLength(length) => length,
        ConsensusDistanceNormalization::WithUBAndLength(_ub, length) => length,
    };
    (dist as f64 / div as f64, link_type)
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

fn calculate_unexplained_bases(
    segments: SegmentedMatrixView,
    region_statistics: &RegionStatistics,
    first_block_segment: usize,
    second_block_segment: usize,
    second_block_target_start: usize,
) -> usize {
    let ub = region_statistics.unexplained_bases[second_block_segment]
        .abs_diff(region_statistics.unexplained_bases[first_block_segment])
        + (second_block_target_start - segments[second_block_segment].start_col);
    ub
}

pub fn gather_join_statistics<T: JoinStatisticsCollector>(
    group: &ProximityGroup,
    initial_segments: &InitialSegments,
    query_lengths: &HashMap<usize, usize>,
    annotation_args: &AnnotationArgs,
) -> Vec<(usize, T)> {
    let alignments = group.alignments;

    let mut query_ids: Vec<usize> = alignments.iter().map(|a| a.query_id).unique().collect();
    query_ids.sort();

    let region_stats = calculate_region_statistics(initial_segments);
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
                    .map(|(i, a)| {
                        let b = Block::from_alignment(a, group.target_start, i, 0.0, 0.0);
                        let seg_i = initial_segments
                            .view_segments()
                            .partition_point(|v| v.start_col <= b.col_start)
                            .saturating_sub(1);
                        (seg_i, b)
                    }),
            )
        })
        .for_each(|(id, compat_alignments)| {
            let mut new_stats = T::new();

            gather_join_statistics_single_family(
                compat_alignments,
                *query_lengths
                    .get(&id)
                    .expect("Query length missing for alignment!"),
                initial_segments.view_segments(),
                &region_stats,
                annotation_args,
                &mut new_stats,
            );

            query_stats.push((id, new_stats));
        });

    query_stats
}

fn link_info(
    first_block: &Block,
    second_block: &Block,
    annotation_args: &AnnotationArgs,
    unexplained_bases: usize,
    consensus_length: usize,
    neighbors: bool,
) -> LinkInfo {
    let (consensus_distance, link_type) = block_consensus_distance(first_block, second_block);
    let joinable = is_joinable(
        block_target_distance(first_block, second_block),
        consensus_distance,
        link_type,
        block_length_on_query(first_block).min(block_length_on_query(second_block)),
        annotation_args,
    );

    LinkInfo {
        target_distance: block_target_distance(first_block, second_block),
        consensus_distance,
        link_type,
        consensus_length,
        unexplained_bases,
        neighbors,
        joinable,
    }
}

fn gather_join_statistics_single_family<'a>(
    compatable_alignments: impl Iterator<Item = (usize, Block)>,
    consensus_length: usize,
    segments: SegmentedMatrixView,
    region_statistics: &RegionStatistics,
    args: &AnnotationArgs,
    join_stats: &mut impl JoinStatisticsCollector,
) {
    let compatable_blocks = compatable_alignments
        .sorted_by_key(|(_u_b, a)| a.col_start)
        .collect_vec();

    compatable_blocks
        .iter()
        .enumerate()
        .for_each(|(idx, (a_segment_idx, a_block))| {
            compatable_blocks[idx + 1..].iter().enumerate().for_each(
                |(idx2, (b_segment_idx, b_block))| {
                    let link_info = &link_info(
                        a_block,
                        b_block,
                        args,
                        calculate_unexplained_bases(
                            segments,
                            region_statistics,
                            *a_segment_idx,
                            *b_segment_idx,
                            b_block.col_start,
                        ),
                        consensus_length,
                        idx + 1 == idx2,
                    );
                    join_stats.add(a_block, b_block, link_info);
                },
            )
        })
}

fn link_assemblies<T: JoinEstimator>(
    graph: &mut HashMap<(SegmentAndDenseRow, SegmentAndDenseRow), Edge>,
    compatable_blocks: impl Iterator<Item = (usize, usize)>,
    consensus_length: usize,
    segments: &SegmentedMatrix,
    query_statistics: &QueryStatistics<T>,
    region_statistics: &RegionStatistics,
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

            let link = link_info(
                a_block,
                b_block,
                args,
                calculate_unexplained_bases(
                    segments,
                    region_statistics,
                    a.0,
                    b.0,
                    b_block.col_start,
                ),
                consensus_length,
                a.0 + 1 == b.0,
            );

            if link.joinable {
                let join_prob = query_statistics
                    .estimator
                    .predict(a_block, b_block, &link, false);

                if join_prob >= args.join_likelihood_threshold {
                    let weight = if a_block.row_idx == b_block.row_idx && ((b.0 - 1) <= a.0) {
                        1.0
                    } else {
                        join_prob
                    };

                    graph.insert(
                        ((a.0, a_block.row_idx), (b.0, b_block.row_idx)),
                        Edge {
                            weight,
                            first_sparse_row: a.1,
                            second_sparse_row: b.1,
                            link_type: link.link_type,
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
        query_lengths: &HashMap<usize, usize>,
        segments: &SegmentedMatrix,
        region_statistics: &RegionStatistics,
        query_statistics: &[QueryStatistics<T>],
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
                    *query_lengths
                        .get(&id)
                        .expect("Unable to find query length for alignment!"),
                    segments,
                    &query_statistics[id],
                    region_statistics,
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
