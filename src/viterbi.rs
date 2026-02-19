use core::f64;
use std::{cmp::Ordering, u16, usize};

use crate::{
    alignment::Strand,
    assembly::{AssemblyGraph, Direction, Edge, LinkType},
    balanced_tree::{AVLIndexSet, SetInsert},
    matrix::Matrix,
    score_params::ScoreParams,
    segments::{Block, BlockType, Segment, SegmentedMatrix},
};

use itertools::{multizip, Itertools};

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
    pub group_index: usize,
    pub prior_block_history: usize,
    pub prior_history: usize,
    pub join_history: usize,
    pub score: f64,
}

impl PartialEq for HistoryInfo {
    fn eq(&self, other: &Self) -> bool {
        self.segment == other.segment
            && self.group_index == other.group_index
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
        match self.group_index.cmp(&other.group_index) {
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
    pub segment_groups: Vec<SegmentGroups>,
    pub segment_offsets: Vec<usize>,
    pub entries: Vec<HistoryEntry>,
}

fn remove_expired_history_entries(
    history: &[HistoryEntry],
    segment_groups: &[SegmentGroups],
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
                let can_join_up_to = segment_groups[val.segment].can_join_up_to(val.group_index);
                if can_join_up_to > current_segment {
                    return current_entry;
                }
                current_entry = val.prior_history;
            }
        }
    }

    0
}

pub fn history_score(entry: &HistoryEntry) -> f64 {
    match entry {
        HistoryEntry::Root => 0.0,
        HistoryEntry::Append(val) | HistoryEntry::Join(val) => val.score,
    }
}

fn prior_history(entry: &HistoryEntry) -> usize {
    match entry {
        HistoryEntry::Root => 0,
        HistoryEntry::Append(val) | HistoryEntry::Join(val) => val.prior_history,
    }
}

fn remove_low_scoring_histories(
    histories: &mut Vec<HistoryEntry>,
    start_offset: usize,
    relative_score_bound: f64,
) {
    let h_len = histories.len();
    let limited_bound = relative_score_bound.min(0.0);
    let best_history_score = histories[start_offset..h_len]
        .iter()
        .map(history_score)
        .max_by(f64::total_cmp)
        .unwrap_or(0.0);

    let mut next_insertion_point = start_offset;

    for next_index in start_offset..h_len {
        if (history_score(&histories[next_index]) - best_history_score) >= limited_bound {
            histories.swap(next_index, next_insertion_point);
            next_insertion_point += 1
        }
    }

    histories.truncate(next_insertion_point);
}

fn limit_history_count(
    histories: &mut Vec<HistoryEntry>,
    start_offset: usize,
    max_history_count: usize,
) {
    if max_history_count == 0 {
        return;
    }

    let h_len = histories.len();

    if (h_len - start_offset) <= max_history_count {
        return;
    }

    // Sort in reverse order so best entries are at the front....
    histories[start_offset..h_len]
        .sort_unstable_by(|a, b| history_score(b).total_cmp(&history_score(a)));

    // Truncate length of histories to the limit, this removes bad histories...
    histories.truncate(start_offset + max_history_count);
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

    histories.truncate(current_unique + 1);
}

fn check_for_forward_link(
    assembly_graph: &AssemblyGraph,
    start_block: &Block,
    later_block: &Block,
) -> Option<f64> {
    // Weight, direction, and link type are ignored for edges...
    let edge = Edge {
        edge_to: later_block.row_idx - 1,
        weight: 0.0,
        direction: Direction::Right,
        link_type: LinkType::Forward,
    };

    // If we find it in either the forward or reverse graph, check it's in front of the start alignment...
    if let Some(e1) = assembly_graph.link_graph[start_block.row_idx - 1].get(&edge) {
        if e1.direction == edge.direction {
            Some(e1.weight)
        } else {
            None
        }
    } else {
        None
    }
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
enum OptionalBlock<'a> {
    Valid(&'a Block),
    After,
}

fn get_group_block_optional<'a>(
    segment: &'a Segment,
    group: &[usize],
    index: usize,
) -> OptionalBlock<'a> {
    if index < group.len() {
        OptionalBlock::Valid(&segment.blocks[group[index]])
    } else {
        OptionalBlock::After
    }
}

fn get_valid_joins_for_current_group(
    current_segment: &Segment,
    current_group: &[usize],
    current_segment_index: usize,
    prior_segment: &Segment,
    prior_group: &[usize],
    assembly_graph: &AssemblyGraph,
    epsilon: f64,
) -> (Vec<usize>, Vec<(usize, f64)>) {
    let mut values: Vec<(f64, usize)> = Vec::with_capacity(current_group.len());

    let mut current_idx = 0;
    let mut prior_idx = 0;

    while current_idx < current_group.len() || prior_idx < prior_group.len() {
        let current_block = &get_group_block_optional(current_segment, current_group, current_idx);
        let prior_block = &get_group_block_optional(prior_segment, prior_group, prior_idx);

        if let (&OptionalBlock::Valid(c_block), &OptionalBlock::Valid(p_block)) =
            (current_block, prior_block)
        {
            if let (Some(query_id1), Some(query_id2)) = (c_block.query_id, p_block.query_id) {
                if query_id1 == query_id2 {
                    if p_block.can_join_up_to >= current_segment_index {
                        // Get cost of connection...
                        if let Some(weight) =
                            check_for_forward_link(assembly_graph, p_block, c_block)
                        {
                            // If greater or equal to, add it to the list...
                            values.push((weight, current_group[current_idx]));
                        }
                    }

                    current_idx += 1;
                    //prior_idx += 1;
                    continue;
                }
            }
        }

        let current_is_smaller = current_block < prior_block;
        current_idx += current_is_smaller as usize;
        prior_idx += !current_is_smaller as usize;
    }

    // TODO: Consider replacing this with linear runtime version that possibly adds more histories, as it apears the join confidences are in almost all cases different, so this optimal algorithm is just wasting time....
    // values.sort_unstable_by(|a, b| a.0.total_cmp(&b.0));
    let last_weight = f64::NEG_INFINITY;

    let split_points = values
        .iter()
        .enumerate()
        .filter_map(|(i, v)| {
            if (last_weight - v.0).abs() > epsilon {
                Some((i, v.0))
            } else {
                None
            }
        })
        .collect_vec();

    /*
    for i in 0..split_points.len() {
        let start = split_points[i].0;
        let end = if i + 1 < split_points.len() {
            split_points[i + 1].0
        } else {
            values.len()
        };
        values[start..end].sort_unstable_by_key(|v| v.1);
    }*/

    (values.iter().map(|v| v.1).collect_vec(), split_points)
}

fn get_valid_appends_for_current_group(
    current_segment: &Segment,
    current_group: &[usize],
    prior_segment: &Segment,
    prior_group: &[usize],
) -> (Vec<usize>, usize) {
    let mut current_values = vec![0; current_group.len()];
    let mut split_point = 0;
    let mut different_insert_point = current_values.len() - 1;

    let mut current_idx = 0;
    let mut prior_idx = 0;

    while current_idx < current_group.len() || prior_idx < prior_group.len() {
        let current_block = &get_group_block_optional(current_segment, current_group, current_idx);
        let prior_block = &get_group_block_optional(prior_segment, prior_group, prior_idx);

        if let (&OptionalBlock::Valid(c_block), &OptionalBlock::Valid(p_block)) =
            (current_block, prior_block)
        {
            if c_block.row_idx == p_block.row_idx {
                current_values[split_point] = current_group[current_idx];
                split_point += 1;

                current_idx += 1;
                //prior_idx += 1;
                continue;
            }
        }

        let is_current_smaller = current_block < prior_block;
        if is_current_smaller {
            current_values[different_insert_point] = current_group[current_idx];
            different_insert_point = different_insert_point.saturating_sub(1);
        }
        current_idx += is_current_smaller as usize;
        prior_idx += !is_current_smaller as usize;
    }

    current_values[split_point..].reverse();

    (current_values, split_point)
}

fn check_for_join(
    histories: &[HistoryEntry],
    segments: &SegmentedMatrix,
    segment_groups: &[SegmentGroups],
    assembly_graph: &AssemblyGraph,
    current_group_reference: (usize, usize),
    start_entry: usize,
    history_depth: usize,
    epsilon: f64,
) -> Option<(usize, Vec<usize>, Vec<(usize, f64)>)> {
    let mut last_hist = start_entry;
    let segment_idx = current_group_reference.0;
    let group_idx = current_group_reference.1;
    let current_group_indexes = segment_groups[segment_idx].get_group(group_idx);

    let has_joinable_blocks = current_group_indexes.iter().any(|&v| {
        let block = &segments[segment_idx].blocks[v];
        block.query_id.is_some()
    });

    if has_joinable_blocks {
        for _ in 0..history_depth {
            match &histories[last_hist] {
                HistoryEntry::Root => return None,
                HistoryEntry::Append(val) | HistoryEntry::Join(val) => {
                    let cur_hist = last_hist;
                    last_hist = val.prior_history;

                    let (valid_blocks, valid_group_splits) = get_valid_joins_for_current_group(
                        &segments[segment_idx],
                        current_group_indexes,
                        segment_idx,
                        &segments[val.segment],
                        segment_groups[val.segment].get_group(val.group_index),
                        assembly_graph,
                        epsilon,
                    );

                    if !valid_blocks.is_empty() {
                        return Some((cur_hist, valid_blocks, valid_group_splits));
                    }
                }
            }
        }
    }

    None
}

#[derive(Debug)]
pub struct SegmentGroups {
    // GROUP VECTORS: All contain number of elements matching total number of groups.
    // Search index, this is are used and enforcing uniqueness of groups. log(n) search, linear insertion time. May need to be update, b-tree would get faster insertion time but may have higher constant factor (and also more memory).
    pub groups_ordered: AVLIndexSet<u16>,
    // Stores the starting offset for each group. The next index is the end of this group, exclusive...
    pub group_offsets: Vec<usize>,
    // For each group, what it can join to...
    pub can_join_to: Vec<usize>,
    // INDEX VECTORS: All contain values for every block in the segment.
    // Groups are not ordered. Indexes within a group ARE ordered, by query id, and then by row index if the query id is the same or non-existent. This order comes from the segmented matrix logic, which sorts each segment in this way.
    pub indexes: Vec<usize>,
}

fn get_offset_range_from_vector(
    offset_vector: &[usize],
    array_length: usize,
    index: usize,
) -> (usize, usize) {
    (
        offset_vector[index],
        if index + 1 < offset_vector.len() {
            offset_vector[index + 1]
        } else {
            array_length
        },
    )
}

impl SegmentGroups {
    pub fn from_segment(segment: &Segment, delta_threshold: f64) -> Self {
        // Sort the blocks by score, collect the indexes for that...
        let mut score_ordered = (0..segment.blocks.len())
            .sorted_by(|&a, &b| {
                let type_cmp = segment.blocks[a]
                    .block_type
                    .cmp(&segment.blocks[b].block_type);
                if matches!(type_cmp, Ordering::Equal) {
                    f64::total_cmp(&segment.blocks[a].confidence, &segment.blocks[b].confidence)
                } else {
                    type_cmp
                }
            })
            .collect_vec();

        // Merge segments if they are close enough in score...
        let mut last_index = 0;
        let multi_segment_offsets = (0..score_ordered.len())
            .filter_map(|next_index| {
                if last_index == next_index // First iteration, just add a group...
                    || segment.blocks[next_index].block_type != segment.blocks[last_index].block_type
                    || (segment.blocks[score_ordered[next_index]].confidence
                        - segment.blocks[score_ordered[last_index]].confidence)
                        .abs()
                        > delta_threshold
                {
                    last_index = next_index;
                    Some(last_index)
                } else {
                    None
                }
            })
            .collect_vec();

        let mut can_join_to_vec = Vec::with_capacity(multi_segment_offsets.len());

        // Sort each group so indexes run in increasing order...
        for index in 0..multi_segment_offsets.len() {
            let (start, end) =
                get_offset_range_from_vector(&multi_segment_offsets, score_ordered.len(), index);
            score_ordered[start..end].sort_unstable();
            let max_join = score_ordered[start..end]
                .iter()
                .map(|&v| segment.blocks[v].can_join_up_to)
                .max()
                .unwrap_or(0);
            can_join_to_vec.push(max_join);
        }

        let mut ordered_group_indexes = AVLIndexSet::new();
        for index in 0..multi_segment_offsets.len() {
            ordered_group_indexes
                .add(|v| {
                    let v_slice = get_offset_range_from_vector(
                        &multi_segment_offsets,
                        score_ordered.len(),
                        v,
                    );
                    let idx_slice = get_offset_range_from_vector(
                        &multi_segment_offsets,
                        score_ordered.len(),
                        index,
                    );

                    v_slice.cmp(&idx_slice)
                })
                .expect("Error adding block group, hit tree capacity...");
        }

        Self {
            can_join_to: can_join_to_vec,
            groups_ordered: ordered_group_indexes,
            group_offsets: multi_segment_offsets,
            indexes: score_ordered,
        }
    }

    pub fn get_first_block<'a>(&self, segment: &'a Segment, group_idx: usize) -> &'a Block {
        &segment.blocks[*self
            .get_group(group_idx)
            .first()
            .expect("Empty group, should be impossible.")]
    }

    pub fn get_range(&self, group_idx: usize) -> (usize, usize) {
        get_offset_range_from_vector(&self.group_offsets, self.index_count(), group_idx)
    }

    pub fn get_group(&self, group_idx: usize) -> &[usize] {
        let range = self.get_range(group_idx);
        &self.indexes[range.0..range.1]
    }

    pub fn can_join_up_to(&self, group_idx: usize) -> usize {
        self.can_join_to[group_idx]
    }

    #[allow(dead_code)]
    pub fn iter_group_ranges(&self) -> impl Iterator<Item = (usize, usize)> + use<'_> {
        (0..self.group_offsets.len()).map(|v| self.get_range(v))
    }

    #[allow(dead_code)]
    pub fn iter_groups(&self) -> impl Iterator<Item = &[usize]> {
        (0..self.group_offsets.len()).map(|v| self.get_group(v))
    }

    pub fn add_group(&mut self, segment: &Segment, new_segment: &[usize]) -> usize {
        let group_offsets = &self.group_offsets;
        let indexes = &self.indexes;

        let insert = self
            .groups_ordered
            .add(|probe_idx| {
                let slice = get_offset_range_from_vector(group_offsets, indexes.len(), probe_idx);
                indexes[slice.0..slice.1].cmp(new_segment)
            })
            .expect("Failed to add a new entry! Ran out of space in the block groups tree!");

        match insert {
            SetInsert::New(idx) => {
                let new_offset = self.index_count();
                self.group_offsets.push(new_offset);
                self.indexes.extend_from_slice(new_segment);
                let max_join = new_segment
                    .iter()
                    .map(|&i| segment.blocks[i].can_join_up_to)
                    .max()
                    .unwrap_or(0);
                self.can_join_to.push(max_join);
                idx
            }
            SetInsert::Found(idx) => idx,
        }
    }

    pub fn group_count(&self) -> usize {
        self.group_offsets.len()
    }

    pub fn index_count(&self) -> usize {
        self.indexes.len()
    }
}

pub fn try_add_history_entry(
    history_entries: &mut Vec<HistoryEntry>,
    segment: &Segment,
    entry: HistoryEntry,
) -> bool {
    if let HistoryEntry::Append(info) | HistoryEntry::Join(info) = &entry {
        if info.score >= segment.absolute_score_bound {
            history_entries.push(entry);
            return true;
        }
    }

    false
}

pub fn history_viterbi_on_segments(
    segments: &SegmentedMatrix,
    score_params: &ScoreParams,
    assembly_graph: &AssemblyGraph,
    history_depth: usize,
    max_history_count: usize,
    min_rel_history_score: f64,
) -> History {
    let corrected_min_history_score = if min_rel_history_score >= 0.0 {
        f64::NEG_INFINITY
    } else {
        min_rel_history_score
    };

    let block_count: usize = segments.iter().map(|s| s.blocks.len()).sum();

    let mut histories: Vec<HistoryEntry> = Vec::with_capacity(block_count + 1);
    let mut seg_offsets: Vec<usize> = Vec::with_capacity(segments.len() + 1);
    let mut segment_groups: Vec<SegmentGroups> = Vec::with_capacity(segments.len());

    histories.push(HistoryEntry::Root);
    seg_offsets.push(0);
    let mut prior_step_end = histories.len();

    // For every segment...
    for segment_idx in 0..segments.len() {
        // Group the blocks in the next segment so ones with basically identical score combine into a single history...
        // TOOD: Parameterize the score threshold...
        let sg = SegmentGroups::from_segment(&segments[segment_idx], 0.001);
        segment_groups.push(sg);

        for group_idx in 0..segment_groups[segment_idx].group_count() {
            for prior_hist_idx in (*seg_offsets.last().unwrap())..prior_step_end {
                // Add a join and no join history...
                let possible_join = check_for_join(
                    &histories,
                    segments,
                    &segment_groups,
                    assembly_graph,
                    (segment_idx, group_idx),
                    prior_hist_idx,
                    history_depth,
                    1e-2,
                );

                let other_index = remove_expired_history_entries(
                    &histories,
                    &segment_groups,
                    segment_idx,
                    prior_hist_idx,
                    history_depth,
                );

                // JOIN HISTORY...
                if let Some((join_index, join_blocks, join_group_starts)) = possible_join {
                    // Clean expired history entries from the join path....
                    let simplified_join_index = remove_expired_history_entries(
                        &histories,
                        &segment_groups,
                        segment_idx,
                        prior_history(&histories[join_index]),
                        history_depth,
                    );

                    for i in 0..join_group_starts.len() {
                        let (group_start, group_transition_cost) = join_group_starts[i];
                        let group_end = if i + 1 < join_group_starts.len() {
                            join_group_starts[i + 1].0
                        } else {
                            join_blocks.len()
                        };

                        let new_group_index = segment_groups[segment_idx].add_group(
                            &segments[segment_idx],
                            &join_blocks[group_start..group_end],
                        );

                        let new_score = history_score(&histories[prior_hist_idx])
                            + group_transition_cost
                            + segment_groups[segment_idx]
                                .get_first_block(&segments[segment_idx], new_group_index)
                                .confidence;

                        try_add_history_entry(
                            &mut histories,
                            &segments[segment_idx],
                            HistoryEntry::Join(HistoryInfo {
                                segment: segment_idx,
                                group_index: new_group_index,
                                prior_block_history: prior_hist_idx,
                                prior_history: simplified_join_index,
                                join_history: join_index,
                                score: new_score,
                            }),
                        );
                    }
                }

                match &histories[prior_hist_idx] {
                    // First step, no cost to start in a row...
                    HistoryEntry::Root => {
                        // Add append event with 0 transition score since were coming from the root...
                        let new_score = history_score(&histories[prior_hist_idx])
                            + segment_groups[segment_idx]
                                .get_first_block(&segments[segment_idx], group_idx)
                                .confidence;

                        try_add_history_entry(
                            &mut histories,
                            &segments[segment_idx],
                            HistoryEntry::Append(HistoryInfo {
                                segment: segment_idx,
                                group_index: group_idx,
                                prior_block_history: prior_hist_idx,
                                prior_history: other_index,
                                join_history: prior_hist_idx,
                                score: new_score,
                            }),
                        );
                    }
                    HistoryEntry::Append(val) | HistoryEntry::Join(val) => {
                        // Can add up to two append events for blocks with multiple alignments:
                        let current_rep_block = &segment_groups[segment_idx]
                            .get_first_block(&segments[segment_idx], group_idx);
                        let prior_rep_block = &segment_groups[val.segment]
                            .get_first_block(&segments[val.segment], val.group_index);
                        let is_skip =
                            current_rep_block.row_idx == 0 || prior_rep_block.row_idx == 0;

                        let (current_blocks, split_point) = get_valid_appends_for_current_group(
                            &segments[segment_idx],
                            segment_groups[segment_idx].get_group(group_idx),
                            &segments[val.segment],
                            segment_groups[val.segment].get_group(val.group_index),
                        );

                        let has_matching = split_point > 0;
                        let new_group = if has_matching {
                            &current_blocks[..split_point]
                        } else {
                            &current_blocks[split_point..]
                        };

                        if !new_group.is_empty() {
                            let new_group_idx = segment_groups[segment_idx]
                                .add_group(&segments[segment_idx], new_group);
                            let new_score = history_score(&histories[prior_hist_idx])
                                + score_params.transition(is_skip, !has_matching)
                                + current_rep_block.confidence;

                            // Add append event for matching blocks, this will have no transition penalty...
                            try_add_history_entry(
                                &mut histories,
                                &segments[segment_idx],
                                HistoryEntry::Append(HistoryInfo {
                                    segment: segment_idx,
                                    group_index: new_group_idx,
                                    prior_block_history: prior_hist_idx,
                                    prior_history: other_index,
                                    join_history: prior_hist_idx,
                                    score: new_score,
                                }),
                            );
                        }
                    }
                };
            }
        }

        remove_low_scoring_histories(
            &mut histories,
            prior_step_end,
            segments[segment_idx]
                .relative_score_bound
                .max(corrected_min_history_score),
        );
        limit_history_count(&mut histories, prior_step_end, max_history_count);
        keep_unique_histories(&mut histories, prior_step_end);

        seg_offsets.push(prior_step_end);
        prior_step_end = histories.len();
    }

    histories.shrink_to_fit();

    History {
        segment_groups,
        segment_offsets: seg_offsets,
        entries: histories,
    }
}

#[derive(Debug, Clone)]
pub struct AnnotatedRange {
    pub query_id: Option<usize>,
    pub row_idx: usize,
    pub col_start: usize,
    pub col_end: usize,
}

#[derive(Debug)]
pub struct RefinedTraceSegment {
    pub annotated: Vec<AnnotatedRange>,
    pub join_index: usize,
    pub score: f64,
    pub segment: usize,
}

impl RefinedTraceSegment {
    pub fn max_bounds(&self) -> (usize, usize) {
        (
            self.annotated
                .iter()
                .map(|v| v.col_start)
                .min()
                .expect("No start column!"),
            self.annotated
                .iter()
                .map(|v| v.col_end)
                .max()
                .expect("No end column!"),
        )
    }
}

fn to_comparable(annot: Option<&AnnotatedRange>) -> Option<(Option<usize>, usize)> {
    if let Some(inner_annot) = annot {
        return Some((inner_annot.query_id, inner_annot.row_idx));
    }
    None
}

fn get_max_history(history_range: &[HistoryEntry], region_idx: usize) -> usize {
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
        .unwrap_or_else(|| {
            panic!(
                "Unable to find a max history, should not be possible! Region: {}",
                region_idx
            )
        })
        .0
}

fn get_matching_blocks<'a>(
    new_blocks: impl Iterator<Item = &'a Block>,
    prior_entries: impl Iterator<Item = &'a AnnotatedRange>,
    exact_match: bool,
) -> impl Iterator<Item = (&'a Block, &'a AnnotatedRange)> {
    let mut start = true;
    let mut prior_iter = prior_entries.fuse();
    let mut prior_val: Option<&AnnotatedRange> = None;

    new_blocks.filter_map(move |b| {
        let new_val = b.to_comparable();

        while start || (prior_val.is_some() && to_comparable(prior_val) < Some(new_val)) {
            if let Some(next_val) = prior_iter.next() {
                prior_val = Some(next_val);
            } else {
                break;
            }
            start = false;
        }

        if let Some(annot_range) = prior_val {
            let is_match = if exact_match || new_val.0.is_none() {
                new_val == (annot_range.query_id, annot_range.row_idx)
            } else {
                new_val.0 == annot_range.query_id
            };

            if is_match {
                Some((b, annot_range))
            } else {
                None
            }
        } else {
            None
        }
    })
}

fn get_possible_extensions<'a>(
    new_blocks: impl Iterator<Item = &'a Block>,
    prior_entries: impl Iterator<Item = &'a AnnotatedRange>,
) -> Vec<AnnotatedRange> {
    get_matching_blocks(new_blocks, prior_entries, true)
        .filter_map(|(b, a)| {
            let mut annot = a.clone();

            if b.target_end >= annot.col_start.saturating_sub(1) {
                annot.col_start = b.target_start;
                return Some(annot);
            }
            None
        })
        .collect_vec()
}

fn get_joinable_extensions<'a>(
    new_blocks: impl Iterator<Item = &'a Block>,
    prior_entries: impl Iterator<Item = &'a AnnotatedRange>,
) -> Vec<AnnotatedRange> {
    get_matching_blocks(new_blocks, prior_entries, false)
        .map(|(b, _a)| AnnotatedRange {
            query_id: b.query_id,
            row_idx: b.row_idx,
            col_start: b.target_start,
            col_end: b.target_end,
        })
        .collect_vec()
}

// Return is the assigned index of the new block, and the new join index to use for the following block...
pub fn history_backtrace_append_block(
    refined_segments: &mut Vec<RefinedTraceSegment>,
    join_stack: &mut Vec<(usize, usize, usize, usize)>,
    blocks: &[&Block],
    current_index: usize,
    join_index: usize,
    score: f64,
    segment: usize,
) -> (Option<usize>, usize) {
    // Case 1: Same row index and touches start of segment in front of it, extend the segment backwards to include this...
    if let Some(ref_seg) = refined_segments.last_mut() {
        let direct_extensions =
            get_possible_extensions(blocks.iter().copied(), ref_seg.annotated.iter());

        if !direct_extensions.is_empty() {
            ref_seg.annotated = direct_extensions;
            return (Some(ref_seg.join_index), join_index);
        }
    }

    if blocks
        .iter()
        .any(|&b| matches!(b.block_type, BlockType::Alignment))
    {
        // Case 2: Is part of a join, use shared join index...
        if let Some(&(check_idx, _hist_idx, stack_idx, group_join_idx)) = join_stack.last() {
            if current_index == check_idx {
                let joins = get_joinable_extensions(
                    blocks.iter().copied(),
                    refined_segments[stack_idx].annotated.iter(),
                );

                // Should not be possible assuming a join was allowed in the first place...
                if joins.is_empty() {
                    panic!(
                        "Annotation from join made with 0 elements! This should not be possible!"
                    );
                }

                refined_segments.push(RefinedTraceSegment {
                    annotated: joins,
                    join_index: group_join_idx,
                    score,
                    segment,
                });

                join_stack.pop();
                return (Some(group_join_idx), join_index);
            }
        }

        // Case 3: New segment not part of a join...
        refined_segments.push(RefinedTraceSegment {
            annotated: blocks
                .iter()
                .map(|&b| AnnotatedRange {
                    query_id: b.query_id,
                    row_idx: b.row_idx,
                    col_start: b.target_start,
                    col_end: b.target_end,
                })
                .collect_vec(),
            join_index,
            score,
            segment,
        });

        return (Some(join_index), join_index + 1);
    }

    // Case 4: Skip state, don't add anything...
    (None, join_index)
}

pub fn backtrace_histories(
    segments: &SegmentedMatrix,
    history: &History,
    region_idx: usize,
) -> Vec<RefinedTraceSegment> {
    debug_assert!(segments.len() == history.segment_offsets.len() - 1);

    let mut refined_segments: Vec<RefinedTraceSegment> = Vec::new();
    let mut join_stack: Vec<(usize, usize, usize, usize)> = Vec::new();

    let last_segment = history.segment_offsets.len() - 1;
    // Find the max in the first row....
    let mut current_idx = history.segment_offsets[last_segment]
        + get_max_history(
            &history.entries[history.segment_offsets[last_segment]..],
            region_idx,
        );
    let mut current_entry = &history.entries[current_idx];
    let mut join_idx: usize = 0;

    while let HistoryEntry::Join(entry_info) | HistoryEntry::Append(entry_info) = current_entry {
        // Append current entry to segment stack...
        let blocks = history.segment_groups[entry_info.segment]
            .get_group(entry_info.group_index)
            .iter()
            .map(|&i| &segments[entry_info.segment].blocks[i])
            .collect_vec();

        // Append block for this entry (or extend prior trace block if this is the same alignment)...
        let new_node_join_index;
        (new_node_join_index, join_idx) = history_backtrace_append_block(
            &mut refined_segments,
            &mut join_stack,
            &blocks,
            current_idx,
            join_idx,
            entry_info.score,
            entry_info.segment,
        );

        // If this is a join, add it so the segment it joins to can be constructed correctly later...
        if let HistoryEntry::Join(_) = current_entry {
            if let Some(new_node_join_index) = new_node_join_index {
                join_stack.push((
                    entry_info.join_history,
                    current_idx,
                    refined_segments.len() - 1,
                    new_node_join_index,
                ));
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
