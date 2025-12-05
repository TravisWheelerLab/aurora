use std::{fmt::Debug, iter::Fuse};

use crate::{
    assembly::AssemblyGraph, chunks::ProximityGroup, matrix::Matrix, score_params::ScoreParams,
    viterbi::TraceSegment, AnnotationArgs,
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
    pub target_start: usize,
    pub target_end: usize,
    pub confidence: f64,
    pub can_join_up_to: usize,
}

#[derive(Debug)]
pub struct Segment {
    pub start_col: usize,
    pub end_col: usize,
    pub blocks: Vec<Block>,
}

// type SegmentedMatrix = Vec<Segment>;
pub type SegmentedMatrix = Vec<Segment>;

#[derive(Eq, Ord, PartialEq, PartialOrd, Copy, Clone, Debug)]
enum MergeEntry<T> {
    Start,
    Some(T),
    End,
}

impl<T> MergeEntry<T> {
    fn is_start(&self) -> bool {
        matches!(self, Self::Start)
    }

    fn is_end(&self) -> bool {
        matches!(self, Self::End)
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

struct MergeIterator<I: Iterator, J: Iterator<Item = I::Item>> {
    iter1: Fuse<I>,
    iter2: Fuse<J>,
    val1: MergeEntry<I::Item>,
    val2: MergeEntry<I::Item>,
    prior_val: MergeEntry<I::Item>,
}

impl<I: Iterator, J: Iterator<Item = I::Item>> MergeIterator<I, J>
where
    I::Item: Copy,
{
    pub fn new(iter1: I, iter2: J) -> Self {
        Self {
            iter1: iter1.fuse(),
            iter2: iter2.fuse(),
            val1: MergeEntry::Start,
            val2: MergeEntry::Start,
            prior_val: MergeEntry::Start,
        }
    }
}

impl<I: Iterator, J: Iterator<Item = I::Item>> Iterator for MergeIterator<I, J>
where
    I::Item: Copy + Ord,
{
    type Item = I::Item;

    fn next(&mut self) -> Option<Self::Item> {
        if let MergeEntry::End = self.prior_val {
            return None;
        }

        let mut next_val: MergeEntry<Self::Item> = self.prior_val;

        while next_val.is_start() || next_val == self.prior_val {
            if self.val1 <= self.val2 {
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

fn unique_merging_iterator<I: Iterator, J: Iterator<Item = I::Item>>(
    list1: I,
    list2: J,
) -> MergeIterator<I, J>
where
    I::Item: Copy,
{
    MergeIterator::new(list1, list2)
}

fn logsumexp(a: f64, b: f64) -> f64 {
    let max = a.max(b);
    let min = a.min(b);
    max + (min - max).exp().ln_1p()
}

pub fn segments_from_matrix_trace(
    group: &ProximityGroup,
    trace_segments: &[TraceSegment],
    confidence_matrix: &Matrix<f64>,
    score_params: &ScoreParams,
    assembly_graph: &AssemblyGraph,
    annotation_args: &AnnotationArgs,
) -> SegmentedMatrix {
    // Matrix should always have at least 1 row (for the skip state)...
    debug_assert!(confidence_matrix.def.num_rows > 0);
    debug_assert!(
        confidence_matrix.def.num_rows == group.alignments.len() + group.tandem_repeats.len() + 1
    );

    let matrix_definition = confidence_matrix.def;
    let mut segments: SegmentedMatrix = Vec::with_capacity(trace_segments.len());

    // Monitor alignment scores, note we'll preallocate for performance...
    let mut row_scores: Vec<f64> = vec![0.0; matrix_definition.num_rows];
    // This tracks the first segment each alignment is found in.
    let mut first_segment_seen: Vec<usize> = vec![0; matrix_definition.num_rows];
    let mut last_segment_seen: Vec<usize> = vec![0; matrix_definition.num_rows];
    // Tracks, for each alignment, if it existed in the prior row...
    let mut prior_val: Vec<usize> = vec![0; matrix_definition.num_rows];

    for (s_idx, seg) in trace_segments.iter().enumerate() {
        // Initialize offsets...
        for i in 0..matrix_definition.num_rows {
            row_scores[i] = 0.0; // ln(1)
                                 // This causes skip state cost to be calculated correctly for the start of a segment...
            prior_val[i] = i;
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
        for column in seg.col_start..=seg.col_end {
            let rows = matrix_definition.active_rows_by_col[column];
            let row_iter = rows
                .iter()
                .enumerate()
                .map(|(score_idx, &row_idx)| (row_idx, Unordered(score_idx)));
            let all_row_iter = valid_rows.iter().map(|&v| (v, Unordered(0)));

            for (row_idx, Unordered(score_idx)) in unique_merging_iterator(row_iter, all_row_iter) {
                let trans_cost = score_params
                    .transition(score_idx == 0, (prior_val[row_idx] > 0) != (score_idx > 0));
                row_scores[row_idx] +=
                    trans_cost + confidence_matrix.get_sparse(score_idx, column).ln();

                // Set for the next column...
                prior_val[row_idx] = score_idx;
            }
        }

        // Compute total confidence of all entries added together for this block (in log space)...
        let total_confidence = valid_rows
            .iter()
            .map(|&row_idx| row_scores[row_idx])
            .reduce(logsumexp)
            .unwrap_or(0.0);
        let min_confidence = annotation_args.min_block_confidence.ln();

        let row_filter = |&&row_idx: &&usize| {
            (row_idx == 0)
                || (seg.ali_id != 0 && (row_scores[row_idx] - total_confidence) > min_confidence)
        };

        valid_rows.iter().filter(row_filter).for_each(|&row_idx| {
            first_segment_seen[row_idx] = if first_segment_seen[row_idx] == 0 {
                segments.len()
            } else {
                first_segment_seen[row_idx]
            };
            last_segment_seen[row_idx] = last_segment_seen[row_idx].max(segments.len());
        });

        let mut new_segment = Segment {
            start_col: seg.col_start,
            end_col: seg.col_end,
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
                        target_start: start,
                        target_end: end,
                        confidence: row_scores[row_idx],
                        can_join_up_to: s_idx,
                    }
                })
                .collect_vec(),
        };
        // Order blocks by query id, then row... This order allows for really fast intersection checks in history code...
        new_segment
            .blocks
            .sort_unstable_by_key(|b| (b.query_id, b.row_idx));

        segments.push(new_segment);
    }

    // Allow each alignment block to farthest segment it can be linked to...
    for (s_idx, seg) in segments.iter_mut().enumerate() {
        for b_idx in 0..seg.blocks.len() {
            let block = &seg.blocks[b_idx];
            if block.row_idx == 0 || s_idx != last_segment_seen[block.row_idx] {
                continue;
            }

            // Skip the skip state and tandem repeats...
            if let BlockType::TandemRepeat | BlockType::Skip = block.block_type {
                continue;
            }
            let mut best_idx = s_idx;

            for e in assembly_graph.link_graph[block.row_idx - 1].iter() {
                best_idx = best_idx.max(first_segment_seen[e.edge_to + 1]);
            }

            seg.blocks[b_idx].can_join_up_to = best_idx;
        }
    }

    segments
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
