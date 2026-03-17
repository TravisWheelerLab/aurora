use core::f64;
use std::{cmp::Ordering, fmt::Debug, iter::Fuse};

use crate::{
    assembly::SegmentAssemblyGraph, chunks::ProximityGroup, matrix::Matrix,
    score_params::ScoreParams, viterbi::TraceSegment, AnnotationArgs,
};
use itertools::Itertools;

// We use this to have a field ignored by comparisons...
#[repr(transparent)]
#[derive(Debug, Copy, Default, Clone)]
pub struct Unordered<T>(pub T);

impl<T> PartialEq for Unordered<T> {
    fn eq(&self, _: &Self) -> bool {
        true
    }
}

impl<T> PartialOrd for Unordered<T> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<T> Eq for Unordered<T> {}

impl<T> Ord for Unordered<T> {
    fn cmp(&self, _: &Self) -> std::cmp::Ordering {
        std::cmp::Ordering::Equal
    }
}

impl<T> std::hash::Hash for Unordered<T> {
    fn hash<H: std::hash::Hasher>(&self, _: &mut H) {}
}

impl<T> From<T> for Unordered<T> {
    fn from(value: T) -> Self {
        Self(value)
    }
}

#[derive(Debug, Eq, PartialEq, PartialOrd, Ord)]
pub enum BlockType {
    Skip,
    Alignment,
    TandemRepeat,
}

#[derive(Debug)]
pub struct Block {
    pub row_idx: usize,
    pub block_type: BlockType,
    pub query_id: Option<usize>,
    pub col_start: usize,
    pub col_end: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub avg_confidence: f64,
    pub alignment_score: f64,
    pub can_join_up_to: usize,
}

impl Ord for Block {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.query_id.cmp(&other.query_id) {
            Ordering::Equal => self.row_idx.cmp(&other.row_idx),
            val => val,
        }
    }
}

impl PartialOrd for Block {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Eq for Block {}

impl PartialEq for Block {
    fn eq(&self, other: &Self) -> bool {
        matches!(self.cmp(other), Ordering::Equal)
    }
}

impl Block {
    pub fn to_comparable(&self) -> (Option<usize>, usize) {
        (self.query_id, self.row_idx)
    }
}

#[derive(Debug)]
pub struct Segment {
    pub start_col: usize,
    pub end_col: usize,
    pub absolute_score_bound: f64,
    pub relative_score_bound: f64,
    pub blocks: Vec<Block>,
}

// type SegmentedMatrix = Vec<Segment>;
pub type SegmentedMatrix = Vec<Segment>;

#[derive(Copy, Clone, Debug)]
enum MergeEntry<T> {
    Start,
    Some(T),
    End,
}

impl<T> Ord for MergeEntry<T> {
    fn cmp(&self, other: &Self) -> Ordering {
        self.to_number().cmp(&other.to_number())
    }
}

impl<T> PartialOrd for MergeEntry<T> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl<T> Eq for MergeEntry<T> {}

impl<T> PartialEq for MergeEntry<T> {
    fn eq(&self, other: &Self) -> bool {
        self.to_number().eq(&other.to_number())
    }
}

impl<T> MergeEntry<T> {
    fn to_number(&self) -> u8 {
        match self {
            MergeEntry::Start => 0,
            MergeEntry::Some(_) => 1,
            MergeEntry::End => 2,
        }
    }

    fn is_start(&self) -> bool {
        matches!(self, Self::Start)
    }
}

impl<T> From<Option<T>> for MergeEntry<T> {
    fn from(value: Option<T>) -> Self {
        match value {
            Some(v) => MergeEntry::Some(v),
            None => MergeEntry::End,
        }
    }
}

impl<T> From<MergeEntry<T>> for Option<T> {
    fn from(value: MergeEntry<T>) -> Self {
        match value {
            MergeEntry::Some(v) => Some(v),
            _ => None,
        }
    }
}

pub struct MergeIterator<
    I: Iterator,
    J: Iterator<Item = I::Item>,
    F: Fn(&I::Item, &I::Item) -> Ordering,
> {
    iter1: Fuse<I>,
    iter2: Fuse<J>,
    val1: MergeEntry<I::Item>,
    val2: MergeEntry<I::Item>,
    prior_val: MergeEntry<I::Item>,
    comparator: F,
}

impl<I: Iterator, J: Iterator<Item = I::Item>, F: Fn(&I::Item, &I::Item) -> Ordering>
    MergeIterator<I, J, F>
where
    I::Item: Copy,
{
    pub fn new(iter1: I, iter2: J, comparator: F) -> Self {
        Self {
            iter1: iter1.fuse(),
            iter2: iter2.fuse(),
            val1: MergeEntry::Start,
            val2: MergeEntry::Start,
            prior_val: MergeEntry::Start,
            comparator,
        }
    }

    fn compare(&self, item1: &MergeEntry<I::Item>, item2: &MergeEntry<I::Item>) -> Ordering {
        match (item1, item2) {
            (MergeEntry::Some(val1), MergeEntry::Some(val2)) => (self.comparator)(val1, val2),
            (val1, val2) => val1.cmp(val2),
        }
    }
}

pub struct InitialSegments {
    segments: SegmentedMatrix,
    initial_trace_scores: Vec<f64>,
}

impl<I: Iterator, J: Iterator<Item = I::Item>, F: Fn(&I::Item, &I::Item) -> Ordering> Iterator
    for MergeIterator<I, J, F>
where
    I::Item: Copy,
{
    type Item = I::Item;

    fn next(&mut self) -> Option<Self::Item> {
        if let MergeEntry::End = self.prior_val {
            return None;
        }

        let mut next_val: MergeEntry<Self::Item> = self.prior_val;

        while next_val.is_start()
            || matches!(self.compare(&self.prior_val, &next_val), Ordering::Equal)
        {
            if matches!(
                self.compare(&self.val1, &self.val2),
                Ordering::Equal | Ordering::Less
            ) {
                next_val = self.val1;
                self.val1 = self.iter1.next().into();
            } else {
                next_val = self.val2;
                self.val2 = self.iter2.next().into();
            }
        }

        self.prior_val = next_val;
        next_val.into()
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let (min1, max1) = self.iter1.size_hint();
        let (min2, max2) = self.iter2.size_hint();
        (
            (min1.max(min2) > 0) as usize,
            match (max1, max2) {
                (Some(v1), Some(v2)) => Some(v1 + v2),
                _ => None,
            },
        )
    }
}

pub fn unique_merging_iterator<I: Iterator, J: Iterator<Item = I::Item>>(
    list1: I,
    list2: J,
) -> MergeIterator<I, J, impl Fn(&I::Item, &I::Item) -> Ordering>
where
    I::Item: Copy + Ord,
{
    MergeIterator::new(list1, list2, |a, b| a.cmp(b))
}

#[derive(Debug)]
struct SegmentInfo {
    can_join_back_to: usize,
    max_block_score: f64,
    first_pass_score: f64,
    max_resolution_segment: usize,
}

fn finalize_segments(
    segments: &mut SegmentedMatrix,
    initial_trace: &[TraceSegment],
    initial_trace_scores: &[f64],
    score_params: &ScoreParams,
    assembly_graph: &SegmentAssemblyGraph,
) {
    // There is some floating point error introduced for the absolute score bound...
    let epsilon = 1e-2;

    let mut visited_segment_info = vec![false; segments.len()];
    let mut segments_info: Vec<SegmentInfo> = Vec::with_capacity(segments.len());

    // Allow each alignment block to be bounded by the farthest segment it can be linked to...
    for ((first, second), edge) in assembly_graph.link_graph.iter() {
        let b = &mut segments[first.0].blocks[edge.first_sparse_row];
        b.can_join_up_to = b.can_join_up_to.max(second.0);
    }

    for (s_idx, seg) in segments.iter_mut().enumerate() {
        // Used for gathering segment info, this is used for computing a lower bound on valid history scores...
        let mut max_block_score = f64::NEG_INFINITY;
        let prior_score = segments_info
            .last()
            .map(|v| v.first_pass_score)
            .unwrap_or(0.0);
        let transition_score = if s_idx > 0 {
            let prior_row = initial_trace[s_idx - 1].row_idx;
            let current_row = initial_trace[s_idx].row_idx;
            score_params.transition(current_row == 0 || prior_row == 0, current_row != prior_row)
        } else {
            0.0
        };

        for b_idx in 0..seg.blocks.len() {
            let block = &seg.blocks[b_idx];
            max_block_score = max_block_score.max(block.alignment_score);
        }

        segments_info.push(SegmentInfo {
            can_join_back_to: s_idx, // Placeholder value...
            max_block_score,
            first_pass_score: prior_score + transition_score + initial_trace_scores[s_idx],
            max_resolution_segment: s_idx, // This is resolved in the next step...
        });
    }

    // Compute farthest back segment each segment can be joined to....
    for (first, second) in assembly_graph.link_graph.keys() {
        let seg_info = &mut segments_info[second.0];
        seg_info.can_join_back_to = seg_info.can_join_back_to.min(first.0);
    }

    // Run DFS-like algorithm to determine sections with joins that can be resolved seperately...
    for s_idx in (0..segments_info.len()).rev() {
        if visited_segment_info[s_idx] {
            continue;
        }
        visited_segment_info[s_idx] = true;

        let mut farthest_back = segments_info[s_idx].can_join_back_to;
        let mut i = s_idx;

        while i > 0 && i > farthest_back {
            visited_segment_info[i - 1] = true;
            segments_info[i - 1].max_resolution_segment = s_idx;
            farthest_back = farthest_back.min(segments_info[i - 1].can_join_back_to);
            i -= 1;
        }
    }

    // Compute lower bounds per segment...
    let mut prior_abs_score = f64::NEG_INFINITY;
    let mut prior_rel_score = f64::NEG_INFINITY;
    let transition_min_max = (0..4)
        .map(|v| score_params.transition(v / 2 == 1, v % 2 == 1))
        .minmax_by(f64::total_cmp)
        .into_option()
        .unwrap_or((0.0, 0.0));
    let max_transition_gap = transition_min_max.1 - transition_min_max.0;

    for s_idx in (0..segments_info.len()).rev() {
        let seg = &segments_info[s_idx];
        let resolved_abs_score = if seg.max_resolution_segment == s_idx {
            seg.first_pass_score
        } else {
            prior_abs_score
        };

        let resolved_rel_score = if seg.max_resolution_segment == s_idx {
            0.0
        } else {
            prior_rel_score
        };

        // Don't remove best score until the next segment...
        segments[s_idx].absolute_score_bound = resolved_abs_score - epsilon;
        prior_abs_score = resolved_abs_score - seg.max_block_score - transition_min_max.1;

        segments[s_idx].relative_score_bound = resolved_rel_score;
        prior_rel_score = resolved_rel_score - max_transition_gap;
    }
}

pub fn segments_from_matrix_trace(
    group: &ProximityGroup,
    trace_segments: &[TraceSegment],
    confidence_matrix: &Matrix<f64>,
    score_params: &ScoreParams,
    annotation_args: &AnnotationArgs,
) -> InitialSegments {
    // Matrix should always have at least 1 row (for the skip state)...
    debug_assert!(confidence_matrix.def.num_rows > 0);
    debug_assert!(
        confidence_matrix.def.num_rows == group.alignments.len() + group.tandem_repeats.len() + 1
    );

    let matrix_definition = confidence_matrix.def;
    let mut segments: SegmentedMatrix = Vec::with_capacity(trace_segments.len());

    // Monitor alignment scores, note we'll preallocate for performance...
    let mut row_scores: Vec<f64> = vec![0.0; matrix_definition.num_rows];
    let mut row_conf_sum: Vec<f64> = vec![0.0; matrix_definition.num_rows];
    let mut row_valid_cell_count: Vec<usize> = vec![0; matrix_definition.num_rows];
    let mut trace_row_scores = Vec::with_capacity(trace_segments.len());
    // This tracks the first segment each alignment is found in.
    let mut alignment_segment_bounds: Vec<Option<(usize, usize)>> =
        vec![None; matrix_definition.num_rows];
    // Tracks, for each alignment, if it existed in the prior row...
    let mut prior_val: Vec<usize> = vec![0; matrix_definition.num_rows];

    for (s_idx, seg) in trace_segments.iter().enumerate() {
        // Initialize offsets...
        for i in 0..matrix_definition.num_rows {
            row_scores[i] = 0.0; // ln(1)
            row_conf_sum[i] = 0.0;
            row_valid_cell_count[i] = 0;
            prior_val[i] = i; // This causes skip state cost to be calculated correctly for the start of a segment...
        }

        // Identify alignments actually in this segment, computations are restricted to these values.
        let valid_rows = matrix_definition
            .col_range_by_logical_row
            .iter()
            .enumerate()
            .filter_map(|(i, (s, e))| {
                if *s <= seg.col_end && *e >= seg.col_start {
                    Some(i)
                } else {
                    None
                }
            })
            .collect_vec();

        // Compute scores and start/end points for all rows in this segment....
        // TODO: This isn't fully correct, if we want it to be we need to track if this segment starts in, and if the segment ends in a skip state to compute correct transitions for history tracing...
        for column in seg.col_start..=seg.col_end {
            let rows = &matrix_definition.active_rows_by_col[column];
            let row_iter = rows
                .iter()
                .enumerate()
                .map(|(score_idx, &row_idx)| (row_idx, Unordered(score_idx)));
            let all_row_iter = valid_rows.iter().map(|&v| (v, Unordered(0)));

            for (row_idx, Unordered(score_idx)) in unique_merging_iterator(row_iter, all_row_iter) {
                let trans_cost = score_params.transition(
                    score_idx == 0 || prior_val[row_idx] == 0,
                    (prior_val[row_idx] > 0) != (score_idx > 0),
                );

                let confidence_value = confidence_matrix.get_sparse(score_idx, column);

                row_scores[row_idx] += trans_cost + confidence_value.ln();
                if score_idx != 0 {
                    row_conf_sum[row_idx] += confidence_value;
                    row_valid_cell_count[row_idx] += 1;
                }

                // Set for the next column...
                prior_val[row_idx] = score_idx;
            }
        }

        // Compute total confidence of all entries added together for this block (in log space)...
        let best_conf = valid_rows
            .iter()
            .map(|&row_idx| row_conf_sum[row_idx] / (row_valid_cell_count[row_idx].max(1) as f64))
            .max_by(f64::total_cmp)
            .unwrap_or(0.0);
        let min_confidence = best_conf * annotation_args.min_block_confidence;

        let row_filter = |&&row_idx: &&usize| {
            (row_idx == 0)
                || (seg.row_idx != 0  // If this is a skip state in the first trace, we force it to be one in the second history trace...
                    && (row_conf_sum[row_idx] / (row_valid_cell_count[row_idx].max(1) as f64))
                        > min_confidence)
        };

        valid_rows.iter().filter(row_filter).for_each(|&row_idx| {
            let new_bound = match alignment_segment_bounds[row_idx] {
                Some((first, last)) => Some((first.min(segments.len()), last.max(segments.len()))),
                None => Some((segments.len(), segments.len())),
            };

            alignment_segment_bounds[row_idx] = new_bound;
        });

        trace_row_scores.push(row_scores[trace_segments[s_idx].row_idx]);

        let mut new_segment = Segment {
            start_col: seg.col_start,
            end_col: seg.col_end,
            absolute_score_bound: f64::NEG_INFINITY,
            relative_score_bound: f64::NEG_INFINITY,
            blocks: valid_rows
                .iter()
                .filter(row_filter)
                .map(|&row_idx| {
                    let start = seg
                        .col_start
                        .max(matrix_definition.col_range_by_logical_row[row_idx].0);
                    let end = seg
                        .col_end
                        .min(matrix_definition.col_range_by_logical_row[row_idx].1);

                    let block_type = if row_idx == 0 {
                        BlockType::Skip
                    } else if row_idx <= group.alignments.len() {
                        BlockType::Alignment
                    } else {
                        BlockType::TandemRepeat
                    };

                    let query_id = match block_type {
                        BlockType::Alignment => Some(group.alignments[row_idx - 1].query_id),
                        _ => None,
                    };

                    Block {
                        row_idx,
                        block_type,
                        query_id,
                        col_start: start,
                        col_end: end,
                        query_start: confidence_matrix.consensus_position(row_idx, start),
                        query_end: confidence_matrix.consensus_position(row_idx, end),
                        avg_confidence: row_conf_sum[row_idx]
                            / (row_valid_cell_count[row_idx].max(1) as f64),
                        alignment_score: row_scores[row_idx],
                        can_join_up_to: s_idx,
                    }
                })
                .collect_vec(),
        };
        // Order blocks by query id, then row... This order allows for really fast intersection checks in history code...
        new_segment.blocks.sort_unstable();

        segments.push(new_segment);
    }

    InitialSegments {
        segments,
        initial_trace_scores: trace_row_scores,
    }
}

pub fn assemble_and_link_segments<'a>(
    proximity_group: &ProximityGroup,
    initial_segments: &'a mut InitialSegments,
    trace_segments: &[TraceSegment],
    score_params: &ScoreParams,
    annotation_args: &AnnotationArgs,
) -> (&'a SegmentedMatrix, SegmentAssemblyGraph) {
    let assembly_graph = SegmentAssemblyGraph::new(
        proximity_group.alignments,
        &initial_segments.segments,
        score_params,
        annotation_args,
    );
    finalize_segments(
        &mut initial_segments.segments,
        trace_segments,
        &initial_segments.initial_trace_scores,
        score_params,
        &assembly_graph,
    );

    (&initial_segments.segments, assembly_graph)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Type used for identifying ordering when values are the same...
    #[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
    #[allow(dead_code)]
    struct LComp(u32, Unordered<u32>);

    fn lcomp(v1: u32, v2: u32) -> LComp {
        LComp(v1, Unordered(v2))
    }

    #[test]
    fn test_merge_iterator() {
        let empty = [0; 0];
        /* Two empty iterators return an empty iterator... */
        assert!(unique_merging_iterator(empty.iter(), empty.iter()).eq(empty.iter()));
        /* One empty iterator just returns the other iterator... */
        assert!(unique_merging_iterator(empty.iter(), [1, 2, 3].iter()).eq([1, 2, 3].iter()));
        assert!(unique_merging_iterator([1].iter(), empty.iter()).eq([1].iter()));

        /* Two sorted lists should merge, with only unique elements returned... */
        assert!(
            unique_merging_iterator([1, 2, 2, 3, 4, 4, 4, 10].iter(), [3, 4, 5, 8].iter())
                .eq([1, 2, 3, 4, 5, 8, 10].iter())
        );

        /* When two lists have the same value, values are taken from the first iterator. */
        assert!(unique_merging_iterator(
            [lcomp(1, 500), lcomp(1, 20), lcomp(2, 4)].iter(),
            [
                lcomp(1, 15),
                lcomp(1, 15),
                lcomp(1, 15),
                lcomp(1, 15),
                lcomp(1, 15),
                lcomp(2, 10),
                lcomp(2, 20),
                lcomp(2, 30)
            ]
            .iter()
        )
        .eq([lcomp(1, 500), lcomp(2, 4)].iter()));

        assert!(
            unique_merging_iterator([lcomp(1, 500)].iter(), [lcomp(1, 15)].iter())
                .eq([lcomp(1, 500)].iter())
        );
    }
}
