use itertools::Itertools;

use crate::{matrix::{Matrix, MatrixDef}, score_params::ScoreParams, viterbi::TraceSegment};
use std::cmp::{Ordering};


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


struct MergeIterator<I: Iterator> where I::Item: Copy + Clone {
    iter1: I,
    iter2: I,
    val1: Option<I::Item>,
    val2: Option<I::Item>,
    prior_val: Option<I::Item>,
    val1_exhasted: bool,
    val2_exhasted: bool
}


impl<I: Iterator> MergeIterator<I> where I::Item: Copy + Clone {
    pub fn new(iter1: I, iter2: I) -> Self {
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


impl<I: Iterator> Iterator for MergeIterator<I> where I::Item: Copy + Clone + Ord {
    type Item = I::Item;
    
    fn next(&mut self) -> Option<Self::Item> {
        let mut next_val: Option<Self::Item> = None;

        // Handle exhasted cases...
        while self.val1.is_some() && self.val2.is_some() && self.prior_val == next_val {
            // Fill iterators with next values...
            if self.val1.is_none() && !self.val1_exhasted {
                self.val1 = self.iter1.next();
            }
            if self.val2.is_none() && !self.val2_exhasted {
                self.val2 = self.iter2.next();
            }

            match (self.val1, self.val2) {
                (Some(v), None) => {
                    self.val1 = None;
                    self.val2_exhasted = true;
                    next_val = Some(v);
                },
                (None, Some(v)) => {
                    self.val2 = None;
                    self.val1_exhasted = true;
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
        }

        self.prior_val = next_val;
        return next_val;
    }
}


fn unique_merging_iterator<T: Iterator>(list1: T, list2: T) {
    return MergeIterator::new(list1, list2);
}


pub fn segments_from_matrix_trace(
    trace_segments: &Vec<TraceSegment>,
    confidence_matrix: &Matrix<f64>,
    score_params: &ScoreParams,
) -> SegmentedMatrix {
    // Matrix should always have at least 1 row (for the skip state)...
    debug_assert!(confidence_matrix.def.num_rows > 0);

    let matrix_definition = confidence_matrix.def;
    let mut segments: SegmentedMatrix = Vec::with_capacity(trace_segments.len());

    // Monitor alignment scores, and start/ends... Note we'll preallocate everything for performance...
    let mut row_scores: Vec<f64> = vec![0.0; matrix_definition.num_rows];
    let mut row_starts: Vec<usize> = vec![0; matrix_definition.num_rows];
    let mut row_ends: Vec<usize> = vec![0; matrix_definition.num_rows];

    let mut found_in_segment: Vec<bool> = vec![false; matrix_definition.num_rows];
    

    for seg in trace_segments.iter().rev() {
        // Initialize offsets...
        for i in 0..matrix_definition.num_rows {
            row_starts[i] = seg.col_start;
            row_ends[i] = seg.col_start;
            row_scores[i] = 0.0;
            found_in_segment[i] = false;
        }

        // Identify alignments actually in this segment...
        let valid_rows = matrix_definition.col_range_by_logical_row
            .iter()
            .enumerate()
            .filter_map(|(i, (s, e))| {
                if *s <= seg.col_end && *e >= seg.col_start {Some(i)} else {None}
            }).collect_vec();
        
        for column in seg.col_start..=seg.col_end {
            // Add new rows to the row tracker...
            for (&row_idx, confidence) in matrix_definition.active_rows_by_col[column].iter().zip(confidence_matrix.data[column]).skip(1) {
                let is_new_row = row_tracker[row_idx] == 0;
                found_rows += is_new_row as usize;  
                row_tracker[row_idx] = if is_new_row {found_rows} else {row_tracker[row_idx]};
                // Flip found flag...
                let in_segment_index = row_tracker[row_idx];
                found_in_segment[in_segment_index] = found_in_segment[0];

                // Add score...
                row_scores[in_segment_index] += score_params.query_loop_score + confidence
            }

            // Iterate all found rows in this segment...
            for row_idx in 0..found_rows {
                let matrix_definition.active_rows_by_col[column][row_idx]
            }
        }
    }

    return segments;
}