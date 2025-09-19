use crate::{alignment::Strand, matrix::Matrix, score_params::ScoreParams, segments::{Block, Segment, SegmentedMatrix}};

use itertools::{izip, multizip};

pub fn viterbi_collapsed(
    confidence_matrix: &Matrix<f64>,
    viterbi_matrix: &mut Matrix<f64>,
    sources_matrix: &mut Matrix<usize>,
    active_cols: &[usize],
    score_params: &ScoreParams,
) {
    let first_col_idx = active_cols[0];

    // initialize the first viterbi column as the
    // log of the first column of confidence values
    for (sparse_row_idx, (confidence_value, viterbi_score, source_value)) in multizip((
        confidence_matrix.col_slice(first_col_idx),
        viterbi_matrix.col_slice_mut(first_col_idx),
        sources_matrix.col_slice_mut(first_col_idx),
    ))
    .enumerate()
    {
        *viterbi_score = confidence_value.ln();
        *source_value = sparse_row_idx;
    }

    for (&col_from_idx, &col_to_idx) in active_cols.iter().zip(active_cols.iter().skip(1)) {
        // TODO: we could keep track of these as they are
        //       being set in the previous loop iteration
        let (sparse_row_idx_of_max_score_in_col_from, &max_score_in_col_from) = viterbi_matrix
            .col_slice(col_from_idx)
            .iter()
            .enumerate()
            // skip the skip state
            .skip(1)
            // filter ghost skip states to prevent
            // paths from jumping to another row
            .filter(|(sparse_row_idx, _)| {
                viterbi_matrix.ali_id_sparse(*sparse_row_idx, col_from_idx) != 0usize
            })
            // find the (sparse_row_idx, score) tuple with the highest score)
            .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
            // if the column is empty, we'll return -inf, which will be ignored
            .unwrap_or((0usize, &-f64::INFINITY));

        // to get the skip state score, we only have to
        // compare the cost of looping vs the cost of
        // jumping from the best score in the previous col
        let skip_loop_score = viterbi_matrix.get_skip(col_from_idx) + score_params.skip_loop_score;
        let query_to_skip_score = max_score_in_col_from + score_params.query_to_skip_score;

        if skip_loop_score > query_to_skip_score {
            viterbi_matrix.set_skip(
                col_to_idx,
                confidence_matrix.get_skip(col_to_idx).ln() + skip_loop_score,
            );
            sources_matrix.set_skip(col_to_idx, 0);
        } else {
            viterbi_matrix.set_skip(
                col_to_idx,
                confidence_matrix.get_skip(col_to_idx).ln() + query_to_skip_score,
            );
            sources_matrix.set_skip(col_to_idx, sparse_row_idx_of_max_score_in_col_from);
        }

        let col_to_length = viterbi_matrix.col_length(col_to_idx);

        // now we need to pull the scores in
        // for every active cell in the column
        (0..col_to_length)
            // skip the skip state, since
            // we already did it
            .skip(1)
            .for_each(|sparse_row_to_idx| {
                // if the row is a ghost skip state, it will have ali_id of 0
                let row_is_ghost =
                    viterbi_matrix.ali_id_sparse(sparse_row_to_idx, col_to_idx) == 0usize;

                // if the row is a ghost, we add massive negative
                // score to all transitions except the loop
                // otherwise, we add 0.0
                let jump_modifier = (row_is_ghost as usize as f64) * -1_000_000.0;

                // now figure out what score we'd
                // get if we came from the skip state
                let skip_to_query_tuple = (
                    viterbi_matrix.get_skip(col_from_idx)
                        + score_params.query_to_skip_score
                        + jump_modifier,
                    0,
                );

                // since the sparse row indices aren't the same
                // between columns, we need to find the logical
                // row idx of the current cell we are computing
                let logical_row_to_idx =
                    viterbi_matrix.sparse_to_logical_row_idx(sparse_row_to_idx, col_to_idx);

                let query_loop_tuple =
                    if viterbi_matrix.contains_cell(logical_row_to_idx, col_from_idx) {
                        let sparse_row_from_idx = viterbi_matrix
                            .logical_to_sparse_row_idx(logical_row_to_idx, col_from_idx);
                        (
                            viterbi_matrix.get_sparse(sparse_row_from_idx, col_from_idx)
                                // NOTE: we're not adding the query loop score 
                                // here since it's currently always ZERO
                                + score_params.skip_loop_score * row_is_ghost as usize as f64,
                            sparse_row_from_idx,
                        )
                    } else {
                        (-f64::INFINITY, 0usize)
                    };

                let query_jump_tuple = (
                    max_score_in_col_from + score_params.query_jump_score + jump_modifier,
                    sparse_row_idx_of_max_score_in_col_from,
                );

                let score_and_source_tuples =
                    [query_loop_tuple, query_jump_tuple, skip_to_query_tuple];

                let (score, source) = score_and_source_tuples.iter().fold(
                    (-f64::INFINITY, 0usize),
                    |acc, &(score, source)| {
                        if score > acc.0 {
                            (score, source)
                        } else {
                            acc
                        }
                    },
                );

                viterbi_matrix.set_sparse(
                    sparse_row_to_idx,
                    col_to_idx,
                    score
                        + confidence_matrix
                            .get_sparse(sparse_row_to_idx, col_to_idx)
                            .ln(),
                );

                sources_matrix.set_sparse(sparse_row_to_idx, col_to_idx, source)
            });
    }
}

///
///
///
///
#[derive(Default, Clone, Debug)]
pub struct TraceStep {
    pub sparse_row_idx: usize,
    pub row_idx: usize,
    pub col_idx: usize,
    pub consensus_pos: usize,
    pub confidence: f64,
    pub strand: Strand,
    pub query_id: usize,
    pub ali_id: usize,
}

///
///
///
///
pub type Trace = Vec<TraceStep>;

///
///
///
///
#[derive(Clone)]
pub struct TraceSegment {
    pub query_id: usize,
    pub ali_id: usize,
    pub row_idx: usize,
    pub col_start: usize,
    pub col_end: usize,
}

///
///
///
///
pub fn trace_segments(trace: &Trace) -> Vec<TraceSegment> {
    let mut trace_segments: Vec<TraceSegment> = vec![];

    let mut start_step = &trace[0];
    trace
        .iter()
        .zip(trace.iter().skip(1))
        .for_each(|(step, next_step)| {
            if step.ali_id != next_step.ali_id || step.col_idx + 1 != next_step.col_idx {
                debug_assert_eq!(start_step.row_idx, step.row_idx);
                debug_assert_eq!(start_step.ali_id, step.ali_id);
                trace_segments.push(TraceSegment {
                    query_id: step.query_id,
                    ali_id: step.ali_id,
                    row_idx: step.row_idx,
                    col_start: start_step.col_idx,
                    col_end: step.col_idx,
                });
                start_step = next_step;
            }
        });

    let last_step = trace.last().unwrap();
    trace_segments.push(TraceSegment {
        query_id: last_step.query_id,
        ali_id: last_step.ali_id,
        row_idx: last_step.row_idx,
        col_start: start_step.col_idx,
        col_end: last_step.col_idx,
    });

    trace_segments
}

pub fn traceback(
    viterbi_matrix: &Matrix<f64>,
    confidence_matrix: &Matrix<f64>,
    sources: &Matrix<usize>,
    active_cols: &[usize],
) -> Trace {
    let col_idx = *active_cols
        .last()
        .expect("active_cols is empty in call to traceback()");

    let (_, sparse_row_idx) = viterbi_matrix.col_slice(col_idx).iter().enumerate().fold(
        (-f64::INFINITY, 0usize),
        |acc, (idx, &score)| {
            if score > acc.0 {
                (score, idx)
            } else {
                acc
            }
        },
    );

    let row_idx = viterbi_matrix.sparse_to_logical_row_idx(sparse_row_idx, col_idx);
    let query_id = viterbi_matrix.query_id_of_row(row_idx);
    let ali_id = viterbi_matrix.ali_id_sparse(sparse_row_idx, col_idx);
    let consensus_pos = viterbi_matrix.consensus_position_sparse(sparse_row_idx, col_idx);
    let confidence = confidence_matrix.get_sparse(sparse_row_idx, col_idx);
    let strand = viterbi_matrix.strand_of_cell_sparse(sparse_row_idx, col_idx);

    let mut trace = vec![TraceStep {
        sparse_row_idx,
        row_idx,
        col_idx,
        consensus_pos,
        confidence,
        strand,
        query_id,
        ali_id,
    }];

    active_cols
        .iter()
        .zip(active_cols.iter().skip(1))
        .rev()
        .for_each(|(&col_idx, &prev_col_idx)| {
            // prev_col_idx is the column we labeled in the previous iteration
            // col_idx is the column we are labeling now
            let prev_step = trace.last().unwrap();

            debug_assert_eq!(prev_step.col_idx, prev_col_idx);

            // source row of the last step
            let sparse_row_idx = sources.get_sparse(prev_step.sparse_row_idx, prev_col_idx);
            let row_idx = viterbi_matrix.sparse_to_logical_row_idx(sparse_row_idx, col_idx);
            let query_id = viterbi_matrix.query_id_of_cell_sparse(sparse_row_idx, col_idx);
            let ali_id = viterbi_matrix.ali_id_sparse(sparse_row_idx, col_idx);
            let consensus_pos = viterbi_matrix.consensus_position_sparse(sparse_row_idx, col_idx);
            let confidence = confidence_matrix.get_sparse(sparse_row_idx, col_idx);
            let strand = viterbi_matrix.strand_of_cell_sparse(sparse_row_idx, col_idx);

            trace.push(TraceStep {
                sparse_row_idx,
                row_idx,
                col_idx,
                consensus_pos,
                confidence,
                strand,
                query_id,
                ali_id,
            })
        });

    trace.reverse();

    trace
}

#[allow(dead_code)]
pub fn print_viterbi_with_sources(viterbi_matrix: &Matrix<f64>, sources_matrix: &Matrix<usize>) {
    (0..viterbi_matrix.num_rows()).for_each(|row_idx| {
        (0..viterbi_matrix.num_cols()).for_each(|col_idx| {
            //
            if viterbi_matrix.contains_cell(row_idx, col_idx) {
                print!("{:12.3} ", viterbi_matrix.get(row_idx, col_idx));
            } else {
                print!("{:>12.3} ", "x");
            }
        });
        println!();

        (0..sources_matrix.num_cols()).for_each(|col_idx| {
            if sources_matrix.contains_cell(row_idx, col_idx) {
                print!("{:12.3} ", sources_matrix.get(row_idx, col_idx));
            } else {
                print!("{:>12.3} ", "x");
            }
        });
        println!();
        println!();
    });
}


pub struct HistoryInfo {
    segment: usize,
    block: usize,
    prior_block: usize,
    prior_history: usize,
    score: f64
}


impl PartialEq for HistoryInfo {
    fn eq(&self, other: &Self) -> bool {
        self.segment == other.segment && self.block == other.block && self.prior_history == other.prior_history
    }
}

impl Eq for HistoryInfo {}

impl Ord for HistoryInfo {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        match self.segment.cmp(&other.segment) {
            core::cmp::Ordering::Equal => {}
            ord => return ord,
        }
        match self.block.cmp(&other.block) {
            core::cmp::Ordering::Equal => {}
            ord => return ord,
        }
        self.prior_history.cmp(&other.prior_history)
    }
}


impl PartialOrd for HistoryInfo {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}


#[derive(PartialEq, Eq, PartialOrd, Ord)]
pub enum HistoryEntry {
    Root,
    Join(HistoryInfo),
    Append(HistoryInfo)
}


pub struct History {
    pub segment_offsets: Vec<usize>,
    pub entries: Vec<HistoryEntry>
}


fn remove_expired_history_entries(
    history: &[HistoryEntry], 
    segments: &SegmentedMatrix, 
    current_segment: usize, 
    start_entry: usize, 
    history_depth: usize
) -> usize {
    let mut current_entry = start_entry;
    
    for _ in 0..history_depth {
        match &history[current_entry] {
            HistoryEntry::Root => {
                return current_entry;
            }
            HistoryEntry::Append(val) | HistoryEntry::Join(val) => {
                let segment = &segments[val.segment].blocks[val.block];
                if segment.can_join_up_to > current_segment {
                    return current_entry;
                }
                current_entry = val.prior_history;
            }
        }
    }

    0
}


fn history_score(entry: &HistoryEntry) -> f64 {
    match entry {
        HistoryEntry::Root => 0.0,
        HistoryEntry::Append(val) | HistoryEntry::Join(val) => val.score
    }
}


fn keep_unique_histories(histories: &mut Vec<HistoryEntry>, start_offset: usize) {
    // Sort top entries in-place...
    let h_len = histories.len();
    histories[start_offset..h_len].sort_unstable();

    let mut current_unique = start_offset;

    for next_idx in start_offset..h_len {
        if histories[next_idx] == histories[current_unique] {
            if history_score(&histories[next_idx]) > history_score(&histories[current_unique]) {
                histories.swap(next_idx, current_unique);
            }
        } else {
            current_unique += 1
        }
    }

    while histories.len() > (current_unique + 1) {
        histories.pop();
    }
}


pub fn history_viterbi_on_segments(
    segments: &SegmentedMatrix, 
    score_params: &ScoreParams,
    history_depth: usize,
) -> History {
    let block_count: usize = segments.iter().map(|s| s.blocks.len()).sum();

    let mut histories: Vec<HistoryEntry> = Vec::with_capacity(block_count + 1);
    let mut seg_offsets: Vec<usize> = Vec::with_capacity(segments.len() + 1);
    
    histories.push(HistoryEntry::Root);
    seg_offsets.push(0);
    let mut prior_step_end = histories.len();

    // For every segment...
    for segment_idx in 0..segments.len() {
        for (block_idx, current_block) in segments[segment_idx].blocks.iter().enumerate() {
            for prior_hist_idx in (*seg_offsets.last().unwrap())..prior_step_end {
                let mut last_hist = prior_hist_idx;
                // Add a join and no join history...
                let join_index = (0..history_depth)
                    .map_while(|_| {
                        match &histories[last_hist] {
                            HistoryEntry::Root => None,
                            HistoryEntry::Append(val) | HistoryEntry::Join(val) => {
                                let cur_hist = last_hist;
                                last_hist = val.prior_history;
                                let blk = &segments[val.segment].blocks[val.block];
                                if blk.query_id.is_some() && current_block.query_id.is_some() && blk.query_id == current_block.query_id && blk.can_join_up_to >= segment_idx {
                                    Some(Some(cur_hist))
                                } else {
                                    Some(None)
                                }
                            }
                        }
                    }).find(|v| v.is_some()).flatten();

                let other_index = remove_expired_history_entries(
                    &histories, 
                    segments, 
                    segment_idx, 
                    prior_hist_idx, 
                    history_depth
                );

                if let Some(mut join_index) = join_index {
                    // Clean expired history entries from the join path....
                    join_index = remove_expired_history_entries(&histories, segments, segment_idx, join_index, history_depth);
                    histories.push(HistoryEntry::Join(HistoryInfo {
                        segment: segment_idx, 
                        block: block_idx, 
                        prior_block: prior_hist_idx, 
                        prior_history: join_index, 
                        score: history_score(&histories[prior_hist_idx]) + score_params.query_loop_score + segments[segment_idx].blocks[block_idx].confidence
                    }));
                }

                // Figure out transition score...
                let trans_score = match &histories[prior_hist_idx] {
                    // First step, no cost to start in a row...
                    HistoryEntry::Root => 0.0,
                    HistoryEntry::Append(val) | HistoryEntry::Join(val) => {
                        let prior_block = &segments[val.segment].blocks[val.block];
                        score_params.transition(current_block.alignment_id == 0, prior_block.alignment_id != current_block.alignment_id)
                    }
                };

                // Add append event...
                histories.push(HistoryEntry::Append(HistoryInfo { 
                    segment: segment_idx, 
                    block: block_idx,
                    prior_block: prior_hist_idx, 
                    prior_history: other_index, 
                    score: history_score(&histories[prior_hist_idx]) + trans_score + segments[segment_idx].blocks[block_idx].confidence
                }));

            }
        }

        keep_unique_histories(&mut histories, prior_step_end);
        seg_offsets.push(prior_step_end);
        prior_step_end = histories.len();
    }

    History { segment_offsets: seg_offsets, entries: histories }
}