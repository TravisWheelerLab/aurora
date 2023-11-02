use itertools::Itertools;
use serde_json::json;

use crate::{
    alignment::{Strand, VecMap},
    matrix::{Matrix, MatrixDef},
    results::Annotation,
    viterbi::Trace,
};

#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Segment {
    /// The index of this segment
    pub idx: usize,
    /// Segments that have been joined will have the same join_id
    pub join_id: usize,
    /// The start of the segment in terms of column indices
    pub start_col_idx: usize,
    /// The end of the segment in terms of column indices
    pub end_col_idx: usize,
}

#[derive(Debug, PartialEq)]
pub struct TraceFragment {
    pub query_id: usize,
    pub join_id: usize,
    pub row_idx: usize,
    pub start_col_idx: usize,
    pub end_col_idx: usize,
    pub consensus_start: usize,
    pub consensus_end: usize,
    pub avg_confidence: f64,
    pub strand: Strand,
}

impl TraceFragment {
    #[allow(dead_code)]
    pub fn overlaps(&self, other: &TraceFragment) -> bool {
        self.start_col_idx <= other.end_col_idx && self.end_col_idx >= other.start_col_idx
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Fragment {
    pub query_id: usize,
    pub row_idx: usize,
    pub start_col_idx: usize,
    pub end_col_idx: usize,
    pub consensus_start: usize,
    pub consensus_end: usize,
    pub avg_confidence: f64,
    pub strand: Strand,
    pub left_distance: usize,
    pub right_distance: usize,
}

impl Fragment {
    pub fn len(&self) -> usize {
        self.end_col_idx - self.start_col_idx + 1
    }
}

#[derive(Clone)]
pub struct Link {
    /// The query ID that is responsible for the link
    pub query_id: usize,
    /// The left segment
    pub left_segment_idx: usize,
    /// The right segment
    pub right_segment_idx: usize,
    /// The sum of the confidence of the fragments responsible for the join
    pub confidence_sum: f64,
}

#[derive(Default)]
pub struct Segments {
    pub num_segments: usize,
    pub segments: Vec<Segment>,
    pub fragments: Vec<Vec<Fragment>>,
    pub trace_fragments: Vec<Vec<TraceFragment>>,
    pub marked_for_removal: Vec<bool>,
    pub links: Vec<Link>,
    pub active_links: Vec<bool>,
}

impl Segments {
    pub fn new(trace: &Trace, confidence_matrix: &Matrix<f64>, active_cols: &[usize]) -> Self {
        // first we're going to zip through the trace
        // to identify all of the called segments
        let mut segment_start_col_idx = *active_cols
            .first()
            .expect("active_cols is empty in call to Segments::new()");
        let mut fragment_start_col_idx = segment_start_col_idx;

        let mut segment_cnt = 0usize;
        let mut segments: Vec<Segment> = vec![];
        let mut trace_fragments: Vec<Vec<TraceFragment>> = vec![vec![]];

        let mut confidence_sum: f64 = 0.0;
        active_cols
            .iter()
            .zip(active_cols.iter().skip(1))
            .map(|(&prev_col_idx, &col_idx)| {
                (prev_col_idx, &trace[prev_col_idx], col_idx, &trace[col_idx])
            })
            .for_each(|(prev_col_idx, prev_step, col_idx, step)| {
                debug_assert!(col_idx > prev_col_idx);
                let is_join = prev_step.join_id == step.join_id;
                let is_contiguous = col_idx - prev_col_idx == 1;
                let is_same_row = prev_step.row_idx == step.row_idx;

                confidence_sum +=
                    confidence_matrix.get_sparse(prev_step.sparse_row_idx, prev_col_idx);

                if !is_same_row || !is_contiguous {
                    trace_fragments[segment_cnt].push(TraceFragment {
                        query_id: prev_step.query_id,
                        join_id: prev_step.join_id,
                        row_idx: prev_step.row_idx,
                        start_col_idx: fragment_start_col_idx,
                        end_col_idx: prev_col_idx,
                        consensus_start: confidence_matrix.consensus_position_sparse(
                            trace[fragment_start_col_idx].sparse_row_idx,
                            fragment_start_col_idx,
                        ),
                        consensus_end: confidence_matrix
                            .consensus_position_sparse(prev_step.sparse_row_idx, prev_col_idx),
                        avg_confidence: confidence_sum
                            / (prev_col_idx - fragment_start_col_idx + 1) as f64,
                        strand: confidence_matrix.strand_of_row(prev_step.row_idx),
                    });
                    fragment_start_col_idx = col_idx;
                    confidence_sum = 0.0;

                    if !is_join {
                        segments.push(Segment {
                            idx: segment_cnt,
                            join_id: prev_step.join_id,
                            start_col_idx: segment_start_col_idx,
                            end_col_idx: prev_col_idx,
                        });
                        trace_fragments.push(vec![]);
                        segment_cnt += 1;
                        segment_start_col_idx = col_idx;
                    }
                }
            });

        let last_col_idx = *active_cols.last().unwrap();
        let last_step = &trace[last_col_idx];
        trace_fragments[segment_cnt].push(TraceFragment {
            query_id: last_step.query_id,
            join_id: last_step.join_id,
            row_idx: last_step.row_idx,
            start_col_idx: fragment_start_col_idx,
            end_col_idx: last_col_idx,
            consensus_start: confidence_matrix.consensus_position_sparse(
                trace[segment_start_col_idx].sparse_row_idx,
                segment_start_col_idx,
            ),
            consensus_end: confidence_matrix
                .consensus_position_sparse(last_step.sparse_row_idx, last_col_idx),
            avg_confidence: confidence_sum / (last_col_idx - segment_start_col_idx + 1) as f64,
            strand: confidence_matrix.strand_of_row(last_step.row_idx),
        });

        let end_col_idx = *active_cols.last().unwrap();
        segments.push(Segment {
            idx: segment_cnt,
            join_id: trace[end_col_idx].join_id,
            start_col_idx: segment_start_col_idx,
            end_col_idx,
        });

        // now that we have the segments identified, we can
        // figure out which fragments are in each segment
        let num_rows = confidence_matrix.num_rows();

        // this will track the rolling sum of confidence values
        let mut confidence: Vec<f64> = vec![0.0; num_rows];
        // we need to know how many valid positions
        // we have so we can take the average
        let mut position_cnts: Vec<usize> = vec![0; num_rows];
        let mut consensus_starts: Vec<usize> = confidence_matrix
            .def
            .strand_by_logical_row
            .iter()
            .map(|s| match s {
                Strand::Forward => usize::MAX,
                Strand::Reverse => 0usize,
                Strand::Unset => panic!(),
            })
            .collect();
        let mut consensus_ends: Vec<usize> = confidence_matrix
            .def
            .strand_by_logical_row
            .iter()
            .map(|s| match s {
                Strand::Forward => 0usize,
                Strand::Reverse => usize::MAX,
                Strand::Unset => panic!(),
            })
            .collect();
        let mut fragments: Vec<Vec<Fragment>> = vec![];

        segments.iter().for_each(|segment| {
            // what: we're going to zip through the segment's column range
            //       to figure out the average confidence and consensus
            //       start & end of each alignment that is inside the segment
            //
            //  why: it's *way* more efficient to do this column by column
            //       instead of row by row because of the column-oriented
            //       memory layout of the Matrix data structure
            (segment.start_col_idx..=segment.end_col_idx).for_each(|col_idx| {
                (0..confidence_matrix.col_length(col_idx))
                    .map(|row_idx| {
                        (
                            row_idx,
                            confidence_matrix.sparse_to_logical_row_idx(row_idx, col_idx),
                        )
                    })
                    .for_each(|(sparse_row_idx, logical_row_idx)| {
                        confidence[logical_row_idx] +=
                            confidence_matrix.get_sparse(sparse_row_idx, col_idx);
                        position_cnts[logical_row_idx] += 1;

                        let strand = confidence_matrix.strand_of_row(logical_row_idx);
                        let consensus_pos =
                            confidence_matrix.consensus_position_sparse(sparse_row_idx, col_idx);
                        let consensus_start = &mut consensus_starts[logical_row_idx];
                        let consensus_end = &mut consensus_ends[logical_row_idx];
                        match strand {
                            Strand::Forward => {
                                *consensus_start = (*consensus_start).min(consensus_pos);
                                *consensus_end = (*consensus_end).max(consensus_pos);
                            }
                            Strand::Reverse => {
                                *consensus_start = (*consensus_start).max(consensus_pos);
                                *consensus_end = (*consensus_end).min(consensus_pos);
                            }
                            Strand::Unset => panic!(),
                        }
                    });
            });

            confidence
                .iter_mut()
                .enumerate()
                .for_each(|(row_idx, val)| {
                    if position_cnts[row_idx] > 0 {
                        *val /= position_cnts[row_idx] as f64
                    }
                });

            let mut segment_fragments = vec![];
            (0..num_rows).for_each(|row_idx| {
                let (row_start, row_end) = confidence_matrix.col_range_by_row(row_idx);

                let segment_includes_row =
                    row_start < segment.end_col_idx && row_end > segment.start_col_idx;

                if segment_includes_row {
                    let start_col_idx = row_start.max(segment.start_col_idx);
                    let end_col_idx = row_end.min(segment.end_col_idx);

                    let query_id = confidence_matrix.query_id_of_row(row_idx);

                    debug_assert!(confidence[row_idx].is_finite());
                    let strand = confidence_matrix.strand_of_row(row_idx);
                    segment_fragments.push(Fragment {
                        query_id,
                        row_idx,
                        start_col_idx,
                        end_col_idx,
                        consensus_start: consensus_starts[row_idx],
                        consensus_end: consensus_ends[row_idx],
                        avg_confidence: confidence[row_idx],
                        strand,
                        left_distance: start_col_idx - segment.start_col_idx,
                        right_distance: segment.end_col_idx - end_col_idx,
                    });
                    confidence[row_idx] = 0.0;
                    position_cnts[row_idx] = 0;
                    match strand {
                        Strand::Forward => {
                            consensus_starts[row_idx] = usize::MAX;
                            consensus_ends[row_idx] = 0;
                        }
                        Strand::Reverse => {
                            consensus_starts[row_idx] = 0;
                            consensus_ends[row_idx] = usize::MAX;
                        }
                        Strand::Unset => panic!(),
                    }
                }
            });

            segment_fragments
                .sort_by(|a, b| a.avg_confidence.partial_cmp(&b.avg_confidence).unwrap());
            segment_fragments.reverse();
            fragments.push(segment_fragments);
        });

        let num_segments = segments.len();

        Self {
            num_segments,
            segments,
            fragments,
            trace_fragments,
            marked_for_removal: vec![true; num_segments],
            links: vec![],
            active_links: vec![],
        }
    }

    pub fn create_links(
        &mut self,
        consensus_join_distance: usize,
        target_join_distance: usize,
        num_skip_loops_eq_to_jump: usize,
        min_fragment_length: usize,
    ) {
        self.trace_fragments
            .iter()
            .enumerate()
            .filter(|(_, frags)| !frags.is_empty())
            .map(|(idx, frags)| {
                (
                    idx,
                    frags[0].query_id,
                    frags[0].strand,
                    frags[0].consensus_start,
                    frags[frags.len() - 1].consensus_end,
                    frags[0].start_col_idx,
                    frags[frags.len() - 1].end_col_idx,
                    frags[0].avg_confidence,
                    frags[frags.len() - 1].avg_confidence,
                )
            })
            .for_each(
                |(
                    segment_idx,
                    query_id,
                    strand,
                    consensus_start,
                    consensus_end,
                    start_col_idx,
                    end_col_idx,
                    left_confidence,
                    right_confidence,
                )| {
                    (0..segment_idx.saturating_sub(1))
                        .map(|idx| (idx, &self.fragments[idx]))
                        .for_each(|(other_segment_idx, fragments_to_the_left)| {
                            fragments_to_the_left
                                .iter()
                                .filter(|f| {
                                    f.query_id == query_id
                                        && f.strand == strand
                                        && f.avg_confidence >= 0.03
                                        && f.len() >= min_fragment_length
                                        && f.right_distance < num_skip_loops_eq_to_jump
                                        && f.consensus_end.abs_diff(consensus_start)
                                            <= consensus_join_distance
                                        && start_col_idx - f.end_col_idx < target_join_distance
                                })
                                .for_each(|f| {
                                    self.links.push(Link {
                                        query_id: f.query_id,
                                        left_segment_idx: other_segment_idx,
                                        right_segment_idx: segment_idx,
                                        confidence_sum: left_confidence + f.avg_confidence,
                                    });
                                })
                        });

                    (segment_idx + 2..self.num_segments)
                        .map(|idx| (idx, &self.fragments[idx]))
                        .for_each(|(other_segment_idx, fragments_to_the_right)| {
                            fragments_to_the_right
                                .iter()
                                .filter(|f| {
                                    f.query_id == query_id
                                        && f.strand == strand
                                        && f.avg_confidence >= 0.03
                                        && f.len() >= min_fragment_length
                                        && f.left_distance < num_skip_loops_eq_to_jump
                                        && f.consensus_start.abs_diff(consensus_end)
                                            <= consensus_join_distance
                                        && f.start_col_idx - end_col_idx < target_join_distance
                                })
                                .for_each(|f| {
                                    self.links.push(Link {
                                        query_id: f.query_id,
                                        left_segment_idx: segment_idx,
                                        right_segment_idx: other_segment_idx,
                                        confidence_sum: right_confidence + f.avg_confidence,
                                    });
                                })
                        });
                },
            );

        self.active_links = vec![true; self.links.len()];
    }

    pub fn process_links(&mut self) {
        fn link_map_fn(args: (usize, &Link)) -> (usize, usize, usize, f64) {
            let (link_idx, link) = args;
            (
                link_idx,
                link.left_segment_idx,
                link.right_segment_idx,
                link.confidence_sum,
            )
        }

        'outer: for (link_idx, start, end, confidence_sum) in
            self.links.iter().enumerate().map(link_map_fn)
        {
            'inner: for (other_idx, other_start, other_end, other_confidence_sum) in
                self.links.iter().enumerate().map(link_map_fn)
            {
                if !self.active_links[other_idx]
                    || (start == other_start || start == other_end || end == other_end)
                    || (start > other_start && end < other_end)
                    || (other_start > start && other_end < end)
                {
                    continue 'inner;
                };

                if start < other_end && end > other_start {
                    if confidence_sum >= other_confidence_sum {
                        self.active_links[link_idx] = false;
                        continue 'outer;
                    } else {
                        self.active_links[other_idx] = false;
                    }
                }
            }
        }

        self.links
            .iter()
            .enumerate()
            .filter(|&(idx, _)| self.active_links[idx])
            .flat_map(|(_, link)| [link.left_segment_idx, link.right_segment_idx])
            .unique()
            .for_each(|segment_idx| {
                self.marked_for_removal[segment_idx] = false;
            });
    }

    pub fn get_results_and_active_cols(
        &self,
        query_names: &VecMap<String>,
        target_name: &str,
        target_start: usize,
        region_id: usize,
    ) -> (Vec<Annotation>, Vec<usize>) {
        let mut results = vec![];
        let mut active_cols = vec![];

        (0..self.num_segments)
            // grab (its fragments, and the removal flag)
            .map(|idx| (&self.trace_fragments[idx], self.marked_for_removal[idx]))
            .for_each(|(fragments, remove_segment)| {
                fragments
                    // for all the fragments in the segment
                    .iter()
                    // the ones that called the segment have a join_id
                    .for_each(|fragment| {
                        if remove_segment {
                            // congratulations! this fragment has
                            // offically become an annotation
                            results.push(Annotation {
                                target_name: target_name.to_string(),
                                target_start: fragment.start_col_idx + target_start,
                                target_end: fragment.end_col_idx + target_start,
                                query_name: query_names.get(fragment.query_id).clone(),
                                query_start: fragment.consensus_start,
                                query_end: fragment.consensus_end,
                                strand: fragment.strand,
                                confidence: fragment.avg_confidence,
                                join_id: fragment.join_id,
                                region_id,
                            });
                        } else {
                            (fragment.start_col_idx..=fragment.end_col_idx).for_each(|col_idx| {
                                debug_assert!(!active_cols.contains(&col_idx));
                                active_cols.push(col_idx);
                            });
                        }
                    });
            });

        active_cols.sort();
        (results, active_cols)
    }
}
