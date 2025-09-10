use itertools::Itertools;
use crate::{collapse::AssemblyGraph, matrix::Matrix, score_params::ScoreParams, viterbi::TraceSegment};


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
    fn partial_cmp(&self, _: &Self) -> Option<std::cmp::Ordering> {
        Some(std::cmp::Ordering::Equal)
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
        return Self(value);
    }
}

impl<T> AsRef<T> for Unordered<T> {
    fn as_ref(&self) -> &T {
        return &self.0;
    }
}

impl<T> AsMut<T> for Unordered<T> {
    fn as_mut(&mut self) -> &mut T {
        return &mut self.0;
    }
}


struct Block {
    pub alignment_id: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub confidence: f64,
    pub can_join_up_to: usize,
}


struct Segment {
    pub start_col: usize,
    pub end_col: usize,
    pub blocks: Vec<Block>
}

// type SegmentedMatrix = Vec<Segment>;
type SegmentedMatrix = Vec<Segment>;


struct MergeIterator<I: Iterator, J: Iterator<Item = I::Item>> {
    iter1: I,
    iter2: J,
    val1: Option<I::Item>,
    val2: Option<I::Item>,
    prior_val: Option<I::Item>,
    val1_exhasted: bool,
    val2_exhasted: bool
}


impl<I: Iterator, J: Iterator<Item = I::Item>> MergeIterator<I, J> where I::Item: Copy {
    pub fn new(iter1: I, iter2: J) -> Self {
        return Self {
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


impl<I: Iterator, J: Iterator<Item = I::Item>> Iterator for MergeIterator<I, J> where I::Item: Copy + Ord {
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
                },
                (None, Some(v)) => {
                    self.val2 = None;
                    next_val = Some(v)
                },
                (Some(v1), Some(v2)) => {
                    if v1 <= v2 {
                        next_val = self.val1;
                        self.val1 = None;
                    }
                    else {
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
        return next_val;
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let (min1, max1) = self.iter1.size_hint();
        let (min2, max2) = self.iter2.size_hint();
        return ((min1.max(min2) > 0) as usize, match (max1, max2) {
            (Some(v1), Some(v2)) => Some(v1 + v2),
            _ => None
        });
    }
}


fn unique_merging_iterator<I: Iterator, J: Iterator<Item = I::Item>>(list1: I, list2: J) -> MergeIterator<I, J> where I::Item: Copy {
    return MergeIterator::new(list1, list2);
}


pub fn segments_from_matrix_trace(
    trace_segments: &Vec<TraceSegment>,
    confidence_matrix: &Matrix<f64>,
    score_params: &ScoreParams,
    assembly_graph: &AssemblyGraph
) -> SegmentedMatrix {
    // Matrix should always have at least 1 row (for the skip state)...
    debug_assert!(confidence_matrix.def.num_rows > 0);

    let matrix_definition = confidence_matrix.def;
    let mut segments: SegmentedMatrix = Vec::with_capacity(trace_segments.len());

    // Monitor alignment scores, and start/ends... Note we'll preallocate everything for performance...
    let mut row_scores: Vec<f64> = vec![0.0; matrix_definition.num_rows];
    let mut exists_prior: Vec<bool> = vec![false; matrix_definition.num_rows];
    

    for (s_idx, seg) in trace_segments.iter().enumerate() {
        // Initialize offsets...
        for i in 0..matrix_definition.num_rows {
            row_scores[i] = 0.0;
        }

        // Identify alignments actually in this segment...
        let valid_rows = matrix_definition.col_range_by_logical_row
            .iter()
            .enumerate()
            .filter_map(|(i, (s, e))| {
                if *s <= seg.col_end && *e >= seg.col_start {Some(i)} else {None}
            }).collect_vec();
        
        // Compute scores and start/end points for all rows in this segment....
        for column in seg.col_start..=seg.col_end {
            let rows = &matrix_definition.active_rows_by_col[column];
            let row_iter = rows.iter().enumerate().map(|(score_idx, &ali_idx)| (ali_idx, Unordered(score_idx)));
            let all_row_iter = valid_rows.iter().map(|&v| (v, Unordered(0 as usize)));
            
            for (ali_idx, Unordered(score_idx)) in unique_merging_iterator(row_iter, all_row_iter) {
                let non_skip = score_idx > 0;
                let state_change = exists_prior[ali_idx] != non_skip;
                // Due to dumb rules...
                let ns = non_skip as u32 as f64;
                let sc = state_change as u32 as f64;
                let nso =  !non_skip as u32 as f64;
                let sco =  !state_change as u32 as f64;
                let jump_score = score_params.query_jump_score * ns + score_params.query_to_skip_score * nso;
                let loop_score = score_params.query_loop_score * ns + score_params.skip_loop_score * nso;

                let trans_cost = jump_score * sc + loop_score * sco;

                row_scores[ali_idx] += trans_cost + confidence_matrix.data[column][score_idx];

                // Set for the next run...
                exists_prior[ali_idx] = non_skip;
            }
        }

        segments.push(Segment { 
            start_col: seg.col_start, 
            end_col: seg.col_end, 
            blocks: valid_rows
                .iter()
                .map(|&ali_id| {
                    let start = seg.col_start.max(matrix_definition.col_range_by_logical_row[ali_id].0);
                    let end = seg.col_end.min(matrix_definition.col_range_by_logical_row[ali_id].1);

                    Block {
                        alignment_id: ali_id,
                        target_start: start,
                        target_end: end,
                        confidence: row_scores[ali_id],
                        can_join_up_to: s_idx,
                    }
                }).collect_vec()
        });
    }

    // Correct 
    /* 
    for seg in segments {
        for block in seg {
            assembly_graph.fwd_map
        }
    }
    */

    return segments;
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::cmp::Ordering;

    #[derive(Debug)]
    struct LComp(u32, Unordered<u32>);

    impl PartialEq for LComp {
        fn eq(&self, other: &Self) -> bool {
            other.0.eq(&other.0)
        }
    }
    impl Eq for LComp {}
    impl PartialOrd for LComp {
        fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
            return self.0.partial_cmp(&other.0);
        }
    }
    impl Ord for LComp {
        fn cmp(&self, other: &Self) -> Ordering { 
            return self.0.cmp(&other.0);
        }
    }
    
    fn lcomp(v1: u32, v2: u32) {
        LComp(v1, Unordered(v2));
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
            unique_merging_iterator([1, 2, 2, 3, 4, 4, 4, 10].iter(), [3, 4, 5, 8].iter()).eq([1, 2, 3, 4, 5, 8, 10].iter())
        );

        /* When two lists have the same value, values are taken from the first iterator first. */
        assert!(
            unique_merging_iterator([lcomp(1, 30), lcomp(1, 40), lcomp(2, 4)].iter(), [lcomp(1, 15)].iter()).eq([lcomp(1, 30), lcomp(2, 4)].iter())
        );
    }
}