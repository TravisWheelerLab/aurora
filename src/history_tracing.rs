use crate::{
    assembly::{Edge, LinkType, SegmentAssemblyGraph},
    score_params::ScoreParams,
    segment_groups::SegmentGroups,
    segments::{Block, BlockType, Segment, SegmentedMatrix},
};
use itertools::Itertools;
use std::cmp::Ordering;

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
        .unwrap_or(f64::NEG_INFINITY);

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

    let max_score = histories[start_offset..h_len]
        .iter()
        .map(history_score)
        .max_by(f64::total_cmp)
        .unwrap_or(f64::NEG_INFINITY);

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

    assert!(histories[start_offset..current_unique + 1]
        .iter()
        .any(|v| history_score(v) >= max_score));

    histories.truncate(current_unique + 1);
}

fn check_for_forward_link(
    assembly_graph: &SegmentAssemblyGraph,
    start_segment: usize,
    later_segment: usize,
    start_block: &Block,
    later_block: &Block,
) -> Option<f64> {
    if let Some(e1) = assembly_graph.link_graph.get(&(
        (start_segment, start_block.row_idx),
        (later_segment, later_block.row_idx),
    )) {
        Some(e1.weight)
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
    prior_segment_index: usize,
    assembly_graph: &SegmentAssemblyGraph,
    epsilon: f64,
) -> (Vec<usize>, Vec<(usize, f64)>, Vec<usize>) {
    let mut values: Vec<(f64, usize)> = Vec::with_capacity(current_group.len());
    let mut remaining_values = Vec::with_capacity(current_group.len());

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
                        if let Some(weight) = check_for_forward_link(
                            assembly_graph,
                            prior_segment_index,
                            current_segment_index,
                            p_block,
                            c_block,
                        ) {
                            // If greater or equal to, add it to the list...
                            values.push((weight, current_group[current_idx]));
                            current_idx += 1;
                            continue;
                        }
                    }

                    remaining_values.push(current_group[current_idx]);
                    current_idx += 1;
                    continue;
                }
            }
        }

        let current_is_smaller = current_block < prior_block;
        if current_is_smaller && current_idx < current_group.len() {
            remaining_values.push(current_group[current_idx]);
        }
        current_idx += current_is_smaller as usize;
        prior_idx += !current_is_smaller as usize;
    }

    // Comments below are for implementation that minimizes the number of newly created groups, but at O(n log(n)) cost...
    // Since in most cases a new group has to be made for all, currently using faster linear method...
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

    (
        values.iter().map(|v| v.1).collect_vec(),
        split_points,
        remaining_values,
    )
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

struct JoinCheckArgs<'a> {
    histories: &'a [HistoryEntry],
    segments: &'a SegmentedMatrix,
    segment_groups: &'a [SegmentGroups],
    assembly_graph: &'a SegmentAssemblyGraph,
    current_group_reference: (usize, usize),
    start_entry: usize,
    history_depth: usize,
    epsilon: f64,
}

type JoinHistoryIndex = usize;
type JoinBlockIndexes = Vec<usize>;
type JoinSegmentOffsetsAndConfidences = Vec<(usize, f64)>;
type PossibleJoins = Vec<(
    JoinHistoryIndex,
    JoinBlockIndexes,
    JoinSegmentOffsetsAndConfidences,
)>;

fn check_for_joins(args: JoinCheckArgs) -> PossibleJoins {
    let JoinCheckArgs {
        histories,
        segments,
        segment_groups,
        assembly_graph,
        current_group_reference,
        start_entry,
        history_depth,
        epsilon,
    } = args;

    let mut last_hist = start_entry;
    let segment_idx = current_group_reference.0;
    let group_idx = current_group_reference.1;
    let mut current_group_indexes = segment_groups[segment_idx].get_group(group_idx).to_vec();

    let has_joinable_blocks = current_group_indexes.iter().any(|&v| {
        let block = &segments[segment_idx].blocks[v];
        block.query_id.is_some()
    });

    let mut possible_joins = Vec::new();

    if has_joinable_blocks {
        for _ in 0..history_depth {
            match &histories[last_hist] {
                HistoryEntry::Root => break,
                HistoryEntry::Append(val) | HistoryEntry::Join(val) => {
                    if current_group_indexes.is_empty() {
                        break;
                    }

                    let cur_hist = last_hist;
                    last_hist = val.prior_history;

                    let (valid_blocks, valid_group_splits, remaining_block_indexes) =
                        get_valid_joins_for_current_group(
                            &segments[segment_idx],
                            &current_group_indexes,
                            segment_idx,
                            &segments[val.segment],
                            segment_groups[val.segment].get_group(val.group_index),
                            val.segment,
                            assembly_graph,
                            epsilon,
                        );
                    current_group_indexes = remaining_block_indexes;

                    if !valid_blocks.is_empty() {
                        possible_joins.push((cur_hist, valid_blocks, valid_group_splits));
                    }
                }
            }
        }
    }

    possible_joins
}

pub fn history_viterbi_on_segments(
    segments: &SegmentedMatrix,
    score_params: &ScoreParams,
    assembly_graph: &SegmentAssemblyGraph,
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
        let sg = SegmentGroups::from_segment(&segments[segment_idx], 0.001);
        segment_groups.push(sg);

        for group_idx in 0..segment_groups[segment_idx].group_count() {
            for prior_hist_idx in (*seg_offsets.last().unwrap())..prior_step_end {
                // Add a join and no join history...
                let possible_joins = check_for_joins(JoinCheckArgs {
                    histories: &histories,
                    segments,
                    segment_groups: &segment_groups,
                    assembly_graph,
                    current_group_reference: (segment_idx, group_idx),
                    start_entry: prior_hist_idx,
                    history_depth,
                    epsilon: 1e-2,
                });

                let other_index = remove_expired_history_entries(
                    &histories,
                    &segment_groups,
                    segment_idx,
                    prior_hist_idx,
                    history_depth,
                );

                // JOIN HISTORIES...
                for (join_index, join_blocks, join_group_starts) in possible_joins.iter() {
                    // Clean expired history entries from the join path....
                    let simplified_join_index = remove_expired_history_entries(
                        &histories,
                        &segment_groups,
                        segment_idx,
                        prior_history(&histories[*join_index]),
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
                                .alignment_score;

                        histories.push(HistoryEntry::Join(HistoryInfo {
                            segment: segment_idx,
                            group_index: new_group_index,
                            prior_block_history: prior_hist_idx,
                            prior_history: simplified_join_index,
                            join_history: *join_index,
                            score: new_score,
                        }));
                    }
                }

                match &histories[prior_hist_idx] {
                    // First step, no cost to start in a row...
                    HistoryEntry::Root => {
                        // Add append event with 0 transition score since were coming from the root...
                        let new_score = history_score(&histories[prior_hist_idx])
                            + segment_groups[segment_idx]
                                .get_first_block(&segments[segment_idx], group_idx)
                                .alignment_score;

                        histories.push(HistoryEntry::Append(HistoryInfo {
                            segment: segment_idx,
                            group_index: group_idx,
                            prior_block_history: prior_hist_idx,
                            prior_history: other_index,
                            join_history: prior_hist_idx,
                            score: new_score,
                        }));
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

                        let new_matching_group = &current_blocks[..split_point];
                        let new_mismatching_group = &current_blocks[split_point..];

                        if !new_matching_group.is_empty() {
                            let new_group_idx = segment_groups[segment_idx]
                                .add_group(&segments[segment_idx], new_matching_group);
                            let new_score = history_score(&histories[prior_hist_idx])
                                + score_params.transition(is_skip, false)
                                + current_rep_block.alignment_score;

                            // Add append event for matching blocks, this will have no transition penalty...
                            histories.push(HistoryEntry::Append(HistoryInfo {
                                segment: segment_idx,
                                group_index: new_group_idx,
                                prior_block_history: prior_hist_idx,
                                prior_history: other_index,
                                join_history: prior_hist_idx,
                                score: new_score,
                            }));
                        }

                        if !new_mismatching_group.is_empty() {
                            let new_group_idx = segment_groups[segment_idx]
                                .add_group(&segments[segment_idx], new_mismatching_group);
                            let new_score = history_score(&histories[prior_hist_idx])
                                + score_params.transition(is_skip, true)
                                + current_rep_block.alignment_score;

                            // Add append event for mismatching blocks, this will have a transition penalty...
                            histories.push(HistoryEntry::Append(HistoryInfo {
                                segment: segment_idx,
                                group_index: new_group_idx,
                                prior_block_history: prior_hist_idx,
                                prior_history: other_index,
                                join_history: prior_hist_idx,
                                score: new_score,
                            }));
                        }
                    }
                };
            }
        }

        if histories.len() == prior_step_end {
            panic!("No histories added by absolute threshold!!!!")
        }

        remove_low_scoring_histories(
            &mut histories,
            prior_step_end,
            segments[segment_idx]
                .relative_score_bound
                .max(corrected_min_history_score),
        );

        if histories.len() == prior_step_end {
            panic!("All histories deleted by relative thresholding!!!!")
        }

        limit_history_count(&mut histories, prior_step_end, max_history_count);
        keep_unique_histories(&mut histories, prior_step_end);

        seg_offsets.push(prior_step_end);
        if histories.len() == prior_step_end {
            panic!("No histories added!!!!")
        }
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
    pub query_start: usize,
    pub query_end: usize,
    pub avg_confidence: f64,
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

            if b.col_end >= annot.col_start.saturating_sub(1) {
                let w = (b.col_end - b.col_start + 1) as f64 / (a.col_end - b.col_start + 1) as f64;

                annot.col_start = b.col_start;
                annot.query_start = b.query_start;
                annot.avg_confidence = a.avg_confidence * (1.0 - w) + b.avg_confidence * w;
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
            col_start: b.col_start,
            col_end: b.col_end,
            query_start: b.query_start,
            query_end: b.query_end,
            avg_confidence: b.avg_confidence,
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
        .any(|&b| matches!(b.block_type, BlockType::Alignment | BlockType::TandemRepeat))
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
                    col_start: b.col_start,
                    col_end: b.col_end,
                    query_start: b.query_start,
                    query_end: b.query_end,
                    avg_confidence: b.avg_confidence,
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
