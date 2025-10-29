use std::cmp::Ordering;

use crate::{
    alignment::Strand,
    assembly::{AssemblyGraph, Direction, Edge, LinkType},
    matrix::Matrix,
    score_params::ScoreParams,
    segments::{Block, BlockType, SegmentedMatrix},
};

use itertools::multizip;

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

#[derive(Debug)]
pub struct HistoryInfo {
    pub segment: usize,
    pub block: usize,
    pub prior_block_history: usize,
    pub prior_history: usize,
    pub join_history: usize,
    pub score: f64,
}

impl PartialEq for HistoryInfo {
    fn eq(&self, other: &Self) -> bool {
        self.segment == other.segment
            && self.block == other.block
            && self.prior_history == other.prior_history
    }
}

impl Eq for HistoryInfo {}

impl Ord for HistoryInfo {
    fn cmp(&self, other: &Self) -> Ordering {
        match self.segment.cmp(&other.segment) {
            Ordering::Equal => {}
            ord => return ord,
        }
        match self.block.cmp(&other.block) {
            Ordering::Equal => {}
            ord => return ord,
        }
        self.prior_history.cmp(&other.prior_history)
    }
}

impl PartialOrd for HistoryInfo {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug)]
pub enum HistoryEntry {
    Root,
    Join(HistoryInfo),
    Append(HistoryInfo),
}

impl PartialEq for HistoryEntry {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Join(l0) | Self::Append(l0), Self::Join(r0) | Self::Append(r0)) => l0 == r0,
            (Self::Root, Self::Root) => true,
            _ => false,
        }
    }
}

impl Eq for HistoryEntry {}

impl Ord for HistoryEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self, other) {
            (Self::Join(l0) | Self::Append(l0), Self::Join(r0) | Self::Append(r0)) => l0.cmp(r0),
            (Self::Root, Self::Root) => Ordering::Equal,
            (Self::Root, _) => Ordering::Less,
            (_, Self::Root) => Ordering::Greater,
        }
    }
}

impl PartialOrd for HistoryEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug)]
pub struct History {
    pub segment_offsets: Vec<usize>,
    pub entries: Vec<HistoryEntry>,
}

fn remove_expired_history_entries(
    history: &[HistoryEntry],
    segments: &SegmentedMatrix,
    current_segment: usize,
    start_entry: usize,
    history_depth: usize,
) -> usize {
    let mut current_entry = start_entry;

    for _ in 0..history_depth {
        match &history[current_entry] {
            HistoryEntry::Root => {
                return current_entry;
            }
            HistoryEntry::Append(val) | HistoryEntry::Join(val) => {
                let blk = &segments[val.segment].blocks[val.block];
                if blk.can_join_up_to > current_segment {
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
        HistoryEntry::Append(val) | HistoryEntry::Join(val) => val.score,
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
            current_unique += 1;
            histories.swap(next_idx, current_unique);
        }
    }

    while histories.len() > (current_unique + 1) {
        histories.pop();
    }
}

fn check_for_forward_link(
    assembly_graph: &AssemblyGraph,
    start_block: &Block,
    later_block: &Block,
) -> bool {
    // Weight, direction, and link type are ignored for edges...
    let edge = Edge {
        edge_to: later_block.row_idx - 1,
        weight: 0.0,
        direction: Direction::Right,
        link_type: LinkType::Forward,
    };

    // If we find it in either the forward or reverse graph, check it's in front of the start alignment...
    if let Some(e1) = assembly_graph.link_graph[start_block.row_idx - 1].get(&edge) {
        e1.direction == edge.direction
    } else {
        false
    }
}

fn check_for_join(
    histories: &[HistoryEntry],
    segments: &SegmentedMatrix,
    assembly_graph: &AssemblyGraph,
    current_block_index: (usize, usize),
    start_entry: usize,
    history_depth: usize,
) -> Option<usize> {
    let mut last_hist = start_entry;
    let current_block = &segments[current_block_index.0].blocks[current_block_index.1];
    let segment_idx = current_block_index.0;

    if let Some(current_query_id) = current_block.query_id {
        for _ in 0..history_depth {
            match &histories[last_hist] {
                HistoryEntry::Root => return None,
                HistoryEntry::Append(val) | HistoryEntry::Join(val) => {
                    let cur_hist = last_hist;
                    last_hist = val.prior_history;
                    let blk = &segments[val.segment].blocks[val.block];

                    if let Some(prior_query_id) = blk.query_id {
                        if prior_query_id == current_query_id
                            && blk.can_join_up_to >= segment_idx
                            && check_for_forward_link(assembly_graph, blk, current_block)
                        {
                            return Some(cur_hist);
                        }
                    }
                }
            }
        }
    }

    None
}

fn get_owning_block(history_entry: &HistoryEntry) -> Option<usize> {
    match history_entry {
        HistoryEntry::Append(val) | HistoryEntry::Join(val) => Some(val.block),
        HistoryEntry::Root => None,
    }
}

pub fn history_viterbi_on_segments(
    segments: &SegmentedMatrix,
    score_params: &ScoreParams,
    assembly_graph: &AssemblyGraph,
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
                // Add a join and no join history...
                let join_index = check_for_join(
                    &histories,
                    segments,
                    assembly_graph,
                    (segment_idx, block_idx),
                    prior_hist_idx,
                    history_depth,
                );

                let other_index = remove_expired_history_entries(
                    &histories,
                    segments,
                    segment_idx,
                    prior_hist_idx,
                    history_depth,
                );

                if let Some(join_index) = join_index {
                    // Clean expired history entries from the join path....
                    let simplified_join_index = remove_expired_history_entries(
                        &histories,
                        segments,
                        segment_idx,
                        join_index,
                        history_depth,
                    );
                    histories.push(HistoryEntry::Join(HistoryInfo {
                        segment: segment_idx,
                        block: block_idx,
                        prior_block_history: prior_hist_idx,
                        prior_history: simplified_join_index,
                        join_history: join_index,
                        score: history_score(&histories[prior_hist_idx])
                            + score_params.query_loop_score
                            + segments[segment_idx].blocks[block_idx].confidence,
                    }));
                }

                // Figure out transition score...
                let trans_score = match &histories[prior_hist_idx] {
                    // First step, no cost to start in a row...
                    HistoryEntry::Root => 0.0,
                    HistoryEntry::Append(val) | HistoryEntry::Join(val) => {
                        let prior_block = &segments[val.segment].blocks[val.block];
                        score_params.transition(
                            current_block.row_idx == 0,
                            prior_block.row_idx != current_block.row_idx,
                        )
                    }
                };

                // Add append event...
                histories.push(HistoryEntry::Append(HistoryInfo {
                    segment: segment_idx,
                    block: block_idx,
                    prior_block_history: prior_hist_idx,
                    prior_history: other_index,
                    join_history: prior_hist_idx,
                    score: history_score(&histories[prior_hist_idx])
                        + trans_score
                        + segments[segment_idx].blocks[block_idx].confidence,
                }));
            }
        }

        keep_unique_histories(&mut histories, prior_step_end);
        seg_offsets.push(prior_step_end);
        prior_step_end = histories.len();
    }

    History {
        segment_offsets: seg_offsets,
        entries: histories,
    }
}

#[derive(Debug)]
pub struct RefinedTraceSegment {
    pub query_id: Option<usize>,
    pub row_idx: usize,
    pub col_start: usize,
    pub col_end: usize,
    pub join_index: usize,
}

fn get_max_history(history_range: &[HistoryEntry]) -> usize {
    history_range
        .iter()
        .map(history_score)
        .enumerate()
        .reduce(|(pi, pscore), (i, score)| {
            if score > pscore {
                (i, score)
            } else {
                (pi, pscore)
            }
        })
        .expect("Unable to find a max history, should not be possible!")
        .0
}

// Return is the assigned index of the new block
pub fn history_backtrace_append_block(
    refined_segments: &mut Vec<RefinedTraceSegment>,
    join_stack: &mut Vec<(usize, usize, usize)>,
    block: &Block,
    current_index: usize,
    join_index: usize,
) -> (Option<usize>, usize) {
    // Case 1: Same row index and touches start of segment in front of it, extend the segment backwards to include this...
    if let Some(ref_seg) = refined_segments.last_mut() {
        if ref_seg.row_idx == block.row_idx
            && block.target_end >= (ref_seg.col_start.saturating_sub(1))
        {
            ref_seg.col_start = block.target_start;
            return (Some(ref_seg.join_index), join_index);
        }
    }

    if let BlockType::TandemRepeat | BlockType::Alignment = block.block_type {
        // Case 2: Is part of a join, use shared join index...
        if let Some(&(check_idx, _hist_idx, group_join_idx)) = join_stack.last() {
            if current_index == check_idx {
                refined_segments.push(RefinedTraceSegment {
                    query_id: block.query_id,
                    row_idx: block.row_idx,
                    col_start: block.target_start,
                    col_end: block.target_end,
                    join_index: group_join_idx,
                });

                join_stack.pop();
                return (Some(group_join_idx), join_index);
            }
        }

        // Case 3: New segment not part of a join...
        refined_segments.push(RefinedTraceSegment {
            query_id: block.query_id,
            row_idx: block.row_idx,
            col_start: block.target_start,
            col_end: block.target_end,
            join_index,
        });

        return (Some(join_index), join_index + 1);
    }

    // Case 4: Skip state, don't add anything...
    (None, join_index)
}

pub fn backtrace_histories(
    segments: &SegmentedMatrix,
    history: &History,
) -> Vec<RefinedTraceSegment> {
    debug_assert!(segments.len() == history.segment_offsets.len() - 1);

    let mut refined_segments: Vec<RefinedTraceSegment> = Vec::new();
    let mut join_stack: Vec<(usize, usize, usize)> = Vec::new();

    let last_segment = history.segment_offsets.len() - 1;
    // Find the max in the first row....
    let mut current_idx = history.segment_offsets[last_segment]
        + get_max_history(&history.entries[history.segment_offsets[last_segment]..]);
    let mut current_entry = &history.entries[current_idx];
    let mut join_idx: usize = 0;

    while let HistoryEntry::Join(entry_info) | HistoryEntry::Append(entry_info) = current_entry {
        // Append current entry to segment stack...
        let block = &segments[entry_info.segment].blocks[entry_info.block];

        // Append block for this entry (or extend prior trace block if this is the same alignment)...
        let new_node_join_index;
        (new_node_join_index, join_idx) = history_backtrace_append_block(
            &mut refined_segments,
            &mut join_stack,
            block,
            current_idx,
            join_idx,
        );

        // If this is a join, add it so the segment it joins to can be constructed correctly later...
        if let HistoryEntry::Join(_) = current_entry {
            if let Some(new_node_join_index) = new_node_join_index {
                join_stack.push((entry_info.join_history, current_idx, new_node_join_index));
            }
        }

        // Go to the next entry in the history...
        current_idx = entry_info.prior_block_history;
        current_entry = &history.entries[current_idx];
    }

    // Reverse so trace segments go from start to end instead of end to start.
    refined_segments.reverse();
    refined_segments
}
