use crate::{
    chunks::ProximityGroup, collapse::AssemblyGraph, matrix::Matrix, score_params::ScoreParams,
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

impl<T> AsRef<T> for Unordered<T> {
    fn as_ref(&self) -> &T {
        &self.0
    }
}

impl<T> AsMut<T> for Unordered<T> {
    fn as_mut(&mut self) -> &mut T {
        &mut self.0
    }
}

pub struct Block {
    pub alignment_id: usize,
    pub query_id: Option<usize>,
    pub target_start: usize,
    pub target_end: usize,
    pub confidence: f64,
    pub can_join_up_to: usize,
}

pub struct Segment {
    pub start_col: usize,
    pub end_col: usize,
    pub blocks: Vec<Block>,
}

// type SegmentedMatrix = Vec<Segment>;
pub type SegmentedMatrix = Vec<Segment>;

struct MergeIterator<I: Iterator, J: Iterator<Item = I::Item>> {
    iter1: I,
    iter2: J,
    val1: Option<I::Item>,
    val2: Option<I::Item>,
    prior_val: Option<I::Item>,
    val1_exhasted: bool,
    val2_exhasted: bool,
}

impl<I: Iterator, J: Iterator<Item = I::Item>> MergeIterator<I, J>
where
    I::Item: Copy,
{
    pub fn new(iter1: I, iter2: J) -> Self {
        Self {
            iter1,
            iter2,
            val1: None,
            val2: None,
            prior_val: None,
            val1_exhasted: false,
            val2_exhasted: false,
        }
    }
}

impl<I: Iterator, J: Iterator<Item = I::Item>> Iterator for MergeIterator<I, J>
where
    I::Item: Copy + Ord,
{
    type Item = I::Item;

    fn next(&mut self) -> Option<Self::Item> {
        let mut next_val: Option<Self::Item> = None;

        while !self.val1_exhasted || !self.val2_exhasted {
            // Fill iterators with next values...
            if self.val1.is_none() && !self.val1_exhasted {
                self.val1 = self.iter1.next();
                self.val1_exhasted = self.val1.is_none();
            }
            if self.val2.is_none() && !self.val2_exhasted {
                self.val2 = self.iter2.next();
                self.val2_exhasted = self.val2.is_none();
            }

            match (self.val1, self.val2) {
                (Some(v), None) => {
                    self.val1 = None;
                    next_val = Some(v);
                }
                (None, Some(v)) => {
                    self.val2 = None;
                    next_val = Some(v)
                }
                (Some(v1), Some(v2)) => {
                    if v1 <= v2 {
                        next_val = self.val1;
                        self.val1 = None;
                    } else {
                        next_val = self.val2;
                        self.val2 = None;
                    }
                }
                _ => {}
            }

            if self.prior_val != next_val {
                break;
            }
        }

        self.prior_val = next_val;
        next_val
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
    debug_assert!(confidence_matrix.def.num_rows == group.alignments.len() + 1);

    let matrix_definition = confidence_matrix.def;
    let mut segments: SegmentedMatrix = Vec::with_capacity(trace_segments.len());

    // Monitor alignment scores, note we'll preallocate for performance...
    let mut row_scores: Vec<f64> = vec![0.0; matrix_definition.num_rows];
    // This tracks the last segment each alignment is found in.
    let mut segment_last_seen: Vec<usize> = vec![0; matrix_definition.num_rows];
    // Tracks, for each alignment, if it existed in the prior row...
    let mut prior_val: Vec<usize> = vec![0; matrix_definition.num_rows];

    for (s_idx, seg) in trace_segments.iter().enumerate() {
        // Initialize offsets...
        for i in 0..matrix_definition.num_rows {
            row_scores[i] = 0.0;
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
            let rows = &matrix_definition.active_rows_by_col[column];
            let row_iter = rows
                .iter()
                .enumerate()
                .map(|(score_idx, &ali_idx)| (ali_idx, Unordered(score_idx)));
            let all_row_iter = valid_rows.iter().map(|&v| (v, Unordered(0)));

            for (ali_idx, Unordered(score_idx)) in unique_merging_iterator(row_iter, all_row_iter) {
                let trans_cost = score_params.transition(score_idx == 0, prior_val[ali_idx] != score_idx);
                row_scores[ali_idx] += trans_cost + confidence_matrix.data[column][score_idx];

                // Set for the next column...
                prior_val[ali_idx] = score_idx;
            }
        }

        // Compute total confidence of all entries added together for this block (in log space)...
        let total_confidence = valid_rows
            .iter()
            .map(|&ali_id| row_scores[ali_id])
            .reduce(logsumexp)
            .unwrap_or(0.0);
        let min_confidence = annotation_args.min_block_confidence.ln();

        valid_rows
            .iter()
            .filter(|&&ali_id| (row_scores[ali_id] - total_confidence) > min_confidence)
            .for_each(|&ali_id| {
                segment_last_seen[ali_id] = segments.len();
            });

        segments.push(Segment {
            start_col: seg.col_start,
            end_col: seg.col_end,
            blocks: valid_rows
                .iter()
                .filter(|&&ali_id| {
                    // Remove sections which have too low of a confidence...
                    (row_scores[ali_id] - total_confidence) > min_confidence
                })
                .map(|&ali_id| {
                    let start = seg
                        .col_start
                        .max(matrix_definition.col_range_by_logical_row[ali_id].0);
                    let end = seg
                        .col_end
                        .min(matrix_definition.col_range_by_logical_row[ali_id].1);

                    let is_alignment = ali_id > 0 && ali_id <= group.alignments.len();

                    Block {
                        alignment_id: ali_id,
                        query_id: if is_alignment {
                            Some(group.alignments[ali_id - 1].query_id)
                        } else {
                            None
                        },
                        target_start: start,
                        target_end: end,
                        confidence: row_scores[ali_id],
                        can_join_up_to: s_idx,
                    }
                })
                .collect_vec(),
        });
    }

    // Allow each alignment block to farthest segment it can be linked to...
    for (s_idx, seg) in segments.iter_mut().enumerate() {
        // 1 is to skip the skip state...
        for b_idx in 0..seg.blocks.len() {
            let block = &seg.blocks[b_idx];

            // Skip the skip state and tandem repeats...
            if block.alignment_id == 0 || block.alignment_id > group.alignments.len() {
                continue;
            }

            let alignment = &group.alignments[block.alignment_id - 1];

            let mut best_idx = s_idx;

            for &compat_al in assembly_graph.fwd_map.get(alignment).into_iter().flatten() {
                best_idx = best_idx.max(segment_last_seen[compat_al.id + 1]);
            }

            seg.blocks[b_idx].can_join_up_to = best_idx;
        }
    }

    segments
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::cmp::Ordering;

    // Type used for identifying ordering when values are the same...
    #[derive(Debug)]
    #[allow(dead_code)]
    struct LComp(u32, Unordered<u32>);

    impl PartialEq for LComp {
        fn eq(&self, other: &Self) -> bool {
            other.0.eq(&other.0)
        }
    }
    impl Eq for LComp {}
    impl PartialOrd for LComp {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            Some(self.cmp(other))
        }
    }
    impl Ord for LComp {
        fn cmp(&self, other: &Self) -> Ordering {
            self.0.cmp(&other.0)
        }
    }

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
            [lcomp(1, 30), lcomp(1, 40), lcomp(2, 4)].iter(),
            [lcomp(1, 15)].iter()
        )
        .eq([lcomp(1, 30), lcomp(2, 4)].iter()));
    }
}
