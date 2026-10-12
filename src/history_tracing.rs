use crate::{
    assembly::{LinkType, SegmentAssemblyGraph, Side},
    score_params::ScoreParams,
    segment_groups::SegmentGroups,
    segments::{Block, BlockType, Segment, SegmentedMatrix},
};
use itertools::Itertools;
use std::cmp::Ordering;

/**
 * Represents a joinable side of a history if n
 */
#[derive(Debug, Clone)]
pub struct JoinSide {
    /// This the history that this side change was caused by. Note this history can already be linked to other blocks durring history building, making this not match the link history.
    /// For non-join events, this is always just equal to the index of the history itself.
    pub caused_by_history: usize,
    /// Stores the segment and block that can be joined to, it's always the farthest for the given side...
    pub link_history: usize,
    pub side: Side,
}

#[derive(Debug)]
pub struct HistoryInfo {
    pub segment: usize,
    pub group_index: usize,
    pub prior_block_history: usize,
    pub prior_history: usize,
    pub join_left_block: JoinSide,
    pub join_right_block: JoinSide,
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

    // If we reach max history depth, simply return the index of the root history...
    0
}

pub fn history_score(entry: &HistoryEntry) -> f64 {
    match entry {
        HistoryEntry::Root => 0.0,
        HistoryEntry::Append(val) | HistoryEntry::Join(val) => val.score,
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

    debug_assert!(histories[start_offset..current_unique + 1]
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
) -> Option<(f64, LinkType)> {
    assembly_graph
        .link_graph
        .get(&(
            (start_segment, start_block.row_idx),
            (later_segment, later_block.row_idx),
        ))
        .map(|e1| (e1.weight, e1.link_type))
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

struct SegmentGroupInfo<'a> {
    segment: &'a Segment,
    group: &'a [usize],
    segment_index: usize,
}

/// Represents a possible join link. Tuple of the originating history, the history actually linked to, the side linked to on that history, and the score of completing that join.
#[derive(Debug, Clone)]
struct JoinLink {
    origin_history: usize,
    linked_history: usize,
    link_side: Side,
    probability: f64,
}

impl JoinLink {
    fn match_within_epsilon(&self, other: &Self, epsilon: f64) -> bool {
        self.origin_history == other.origin_history
            && self.linked_history == other.linked_history
            && self.link_side == other.link_side
            && ((self.probability - other.probability).abs() < epsilon)
    }
}

impl PartialEq for JoinLink {
    fn eq(&self, other: &Self) -> bool {
        matches!(self.cmp(other), Ordering::Equal)
    }
}

impl Eq for JoinLink {}

impl PartialOrd for JoinLink {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for JoinLink {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.origin_history
            .cmp(&other.origin_history)
            .then_with(|| self.linked_history.cmp(&other.linked_history))
            .then_with(|| self.link_side.cmp(&other.link_side))
            .then_with(|| self.probability.total_cmp(&other.probability))
    }
}

fn get_valid_joins_for_current_group(
    current_group_info: SegmentGroupInfo,
    prior_group_info: SegmentGroupInfo,
    sides_solved: &mut [[Option<JoinLink>; 2]],
    prior_origin_history: usize,
    prior_link_history: usize,
    prior_linkable_side: Side,
    assembly_graph: &SegmentAssemblyGraph,
) -> bool {
    let SegmentGroupInfo {
        segment: current_segment,
        group: current_group,
        segment_index: current_segment_index,
    } = current_group_info;

    let SegmentGroupInfo {
        segment: prior_segment,
        group: prior_group,
        segment_index: prior_segment_index,
    } = prior_group_info;

    let mut current_idx: usize = 0;
    let mut prior_idx: usize = 0;
    let mut solved_current: usize = 0;

    while current_idx < current_group.len() || prior_idx < prior_group.len() {
        let current_block = &get_group_block_optional(current_segment, current_group, current_idx);
        let prior_block = &get_group_block_optional(prior_segment, prior_group, prior_idx);

        if let (&OptionalBlock::Valid(c_block), &OptionalBlock::Valid(p_block)) =
            (current_block, prior_block)
        {
            // If both sides have been linked, skip this block...
            if sides_solved[current_idx].iter().all(|v| v.is_some()) {
                solved_current += 1;
                current_idx += 1;
                continue;
            }

            if let (Some(query_id1), Some(query_id2)) = (c_block.query_id, p_block.query_id) {
                if query_id1 == query_id2 {
                    if p_block.can_join_up_to >= current_segment_index {
                        // Get cost of connection...
                        if let Some((weight, link_type)) = check_for_forward_link(
                            assembly_graph,
                            prior_segment_index,
                            current_segment_index,
                            p_block,
                            c_block,
                        ) {
                            let (proposed_prior_side, current_side) = link_type.get_linked_sides();

                            // Check if the sides are actually available to link (not already taken by another join)...
                            if &proposed_prior_side == &prior_linkable_side
                                && sides_solved[current_idx][current_side.to_index()].is_none()
                            {
                                sides_solved[current_idx][current_side.to_index()] =
                                    Some(JoinLink {
                                        origin_history: prior_origin_history,
                                        linked_history: prior_link_history,
                                        link_side: prior_linkable_side,
                                        probability: weight,
                                    });

                                solved_current +=
                                    sides_solved[current_idx][current_side.flip().to_index()]
                                        .is_some() as usize;
                            }
                        }
                    }

                    current_idx += 1;
                    continue;
                }
            }
        }

        let current_is_smaller = current_block < prior_block;
        current_idx += current_is_smaller as usize;
        prior_idx += !current_is_smaller as usize;
    }

    solved_current == current_group.len()
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

type PossibleJoinLinks = Vec<[Option<JoinLink>; 2]>;
type JoinBlockIndexes = Vec<usize>;
type SplitPoints = Vec<usize>;

fn check_for_joins(args: JoinCheckArgs) -> (PossibleJoinLinks, JoinBlockIndexes, SplitPoints) {
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
    let current_group_indexes = segment_groups[segment_idx].get_group(group_idx).to_vec();

    let mut sides_solved: Vec<[Option<JoinLink>; 2]> =
        vec![[None, None]; current_group_indexes.len()];

    let has_joinable_blocks = current_group_indexes.iter().any(|&v| {
        let block = &segments[segment_idx].blocks[v];
        block.query_id.is_some()
    });

    if has_joinable_blocks {
        for _ in 0..history_depth {
            match &histories[last_hist] {
                HistoryEntry::Root => break,
                HistoryEntry::Append(val) | HistoryEntry::Join(val) => {
                    let cur_hist = last_hist;
                    last_hist = val.prior_history;

                    // Prefer linking to end of prior history first...
                    let all_resolved = get_valid_joins_for_current_group(
                        SegmentGroupInfo {
                            segment: &segments[segment_idx],
                            group: &current_group_indexes,
                            segment_index: segment_idx,
                        },
                        SegmentGroupInfo {
                            segment: &segments[val.segment],
                            group: segment_groups[val.segment].get_group(val.group_index),
                            segment_index: val.segment,
                        },
                        &mut sides_solved,
                        cur_hist,
                        val.join_right_block.link_history,
                        val.join_right_block.side,
                        assembly_graph,
                    );

                    if all_resolved {
                        break;
                    }

                    // Now check if we can link to the start...
                    let all_resolved = get_valid_joins_for_current_group(
                        SegmentGroupInfo {
                            segment: &segments[segment_idx],
                            group: &current_group_indexes,
                            segment_index: segment_idx,
                        },
                        SegmentGroupInfo {
                            segment: &segments[val.segment],
                            group: segment_groups[val.segment].get_group(val.group_index),
                            segment_index: val.segment,
                        },
                        &mut sides_solved,
                        cur_hist,
                        val.join_left_block.link_history,
                        val.join_left_block.side,
                        assembly_graph,
                    );

                    if all_resolved {
                        break;
                    }
                }
            }
        }
    }

    // Resolve into possible joins, merging joins with the same score...
    let (block_indexes, join_links): (Vec<_>, Vec<_>) = sides_solved
        .iter()
        .enumerate()
        .sorted_by_key(|v| v.1)
        .filter_map(|(i, v)| {
            if v[0].is_none() && v[1].is_none() {
                None
            } else {
                Some((current_group_indexes[i], v.clone()))
            }
        })
        .unzip();

    let mut prior_val: Option<&[Option<JoinLink>; 2]> = None;

    let split_indexes = join_links
        .iter()
        .enumerate()
        .filter_map(|(i, v)| {
            if let Some(val) = prior_val {
                if v.iter().zip(val.iter()).any(|(v1, v2)| match (v1, v2) {
                    (Some(jl1), Some(jl2)) => !jl1.match_within_epsilon(jl2, epsilon),
                    (None, None) => false,
                    _ => true,
                }) {
                    prior_val = Some(v);
                    Some(i)
                } else {
                    None
                }
            } else {
                prior_val = Some(v);
                Some(i)
            }
        })
        .collect_vec();

    (join_links, block_indexes, split_indexes)
}

const JOIN_SIDE_SELF_PLACEHOLDER: usize = 0;

fn finalize_history_indexes_for_join_side(side: &mut JoinSide, hist_idx: usize) {
    if side.caused_by_history == JOIN_SIDE_SELF_PLACEHOLDER {
        side.caused_by_history = hist_idx;
    }
    if side.link_history == JOIN_SIDE_SELF_PLACEHOLDER {
        side.link_history = hist_idx;
    }
}

fn finalize_join_sides(histories: &mut [HistoryEntry], last_segment_start: usize) {
    for hist_idx in last_segment_start..histories.len() {
        if let HistoryEntry::Append(info) | HistoryEntry::Join(info) = &mut histories[hist_idx] {
            finalize_history_indexes_for_join_side(&mut info.join_left_block, hist_idx);
            finalize_history_indexes_for_join_side(&mut info.join_right_block, hist_idx);
        }
    }
}

fn iterate_from_split_points<'a>(
    splits: &'a [usize],
    iterator_length: usize,
) -> impl Iterator<Item = (usize, usize)> + use<'a> {
    (0..splits.len()).map(move |i| {
        (
            splits[i],
            if i + 1 < splits.len() {
                splits[i + 1]
            } else {
                iterator_length
            },
        )
    })
}

fn get_history_link(entry: &HistoryEntry, entry_index: usize, side: Side) -> JoinSide {
    match entry {
        HistoryEntry::Root => JoinSide {
            caused_by_history: entry_index,
            link_history: entry_index,
            side: side,
        },
        HistoryEntry::Append(val) | HistoryEntry::Join(val) => match side {
            Side::Left => val.join_left_block.clone(),
            Side::Right => val.join_right_block.clone(),
        },
    }
}

fn get_join_endpoints_from_links(
    histories: &[HistoryEntry],
    left_join: Option<&JoinLink>,
    right_join: Option<&JoinLink>,
) -> (JoinSide, JoinSide) {
    let left_side = if let Some(join) = left_join {
        let mut link = get_history_link(
            &histories[join.linked_history],
            join.linked_history,
            join.link_side.flip(),
        );
        link.caused_by_history = join.origin_history;

        link
    } else {
        JoinSide {
            caused_by_history: JOIN_SIDE_SELF_PLACEHOLDER,
            link_history: JOIN_SIDE_SELF_PLACEHOLDER,
            side: Side::Left,
        }
    };

    let right_side = if let Some(join) = right_join {
        let mut link = get_history_link(
            &histories[join.linked_history],
            join.linked_history,
            join.link_side.flip(),
        );
        link.caused_by_history = join.origin_history;

        link
    } else {
        JoinSide {
            caused_by_history: JOIN_SIDE_SELF_PLACEHOLDER,
            link_history: JOIN_SIDE_SELF_PLACEHOLDER,
            side: Side::Right,
        }
    };

    (left_side, right_side)
}

fn prior_history_index(entry: &HistoryEntry) -> usize {
    if let HistoryEntry::Append(info) | HistoryEntry::Join(info) = entry {
        info.prior_history
    } else {
        // 0 for the root of the histories...
        0
    }
}

fn add_single_join(
    histories: &mut Vec<HistoryEntry>,
    segments: &[Segment],
    segment_groups: &mut [SegmentGroups],
    prior_hist_idx: usize,
    segment_idx: usize,
    join_blocks: &[usize],
    left_join_link: Option<&JoinLink>,
    right_join_link: Option<&JoinLink>,
    score_params: &ScoreParams,
    history_depth: usize,
) {
    if right_join_link.is_none() && left_join_link.is_none() {
        return;
    }

    // Create a new group for the blocks in the join...
    let new_group_index =
        segment_groups[segment_idx].add_group(&segments[segment_idx], join_blocks);

    let join_prior_index = prior_history_index(
        &histories[match (left_join_link, right_join_link) {
            (Some(lv), Some(rv)) => lv.origin_history.min(rv.origin_history),
            (Some(v), None) | (None, Some(v)) => v.origin_history,
            _ => panic!("Unreachable branch here, something went really wrong..."),
        }],
    );

    // Clean expired history entries from the join path....
    let simplified_join_index = remove_expired_history_entries(
        &histories,
        &segment_groups,
        segment_idx,
        join_prior_index,
        history_depth,
    );

    let prior_is_skip = match &histories[prior_hist_idx] {
        HistoryEntry::Append(val) | HistoryEntry::Join(val) => val.group_index == 0,
        HistoryEntry::Root => true,
    };

    let right_score = right_join_link
        .as_ref()
        .map(|v| score_params.join_transition(prior_is_skip, v.probability))
        .unwrap_or(0.0);
    let left_score = left_join_link
        .as_ref()
        .map(|v| score_params.join_transition(prior_is_skip, v.probability))
        .unwrap_or(0.0);
    // If two joins, we incurred a expensive query-to-query jump in the past, so now we undo that cost...
    let bonus = if left_join_link.is_some() && right_join_link.is_some() {
        -score_params.query_jump_score
    } else {
        0.0
    };
    let transition_cost = right_score + left_score + bonus;

    let new_score = history_score(&histories[prior_hist_idx])
        + transition_cost
        + segment_groups[segment_idx]
            .get_first_block(&segments[segment_idx], new_group_index)
            .alignment_score;

    let (left_join_side, right_join_side) =
        get_join_endpoints_from_links(histories, left_join_link, right_join_link);

    histories.push(HistoryEntry::Join(HistoryInfo {
        segment: segment_idx,
        group_index: new_group_index,
        prior_block_history: prior_hist_idx,
        prior_history: simplified_join_index,
        join_left_block: left_join_side,
        join_right_block: right_join_side,
        score: new_score,
    }));
}

fn add_joins_to_history(
    histories: &mut Vec<HistoryEntry>,
    segments: &[Segment],
    segment_groups: &mut [SegmentGroups],
    prior_hist_idx: usize,
    segment_idx: usize,
    join_links: &[[Option<JoinLink>; 2]],
    join_block_indexes: &[usize],
    join_split_points: &[usize],
    score_params: &ScoreParams,
    history_depth: usize,
) {
    // JOIN HISTORIES...
    for (join_start, join_end) in iterate_from_split_points(join_split_points, join_links.len()) {
        // All joins in group should match, so just use the first one...
        let [left_join, right_join] = &join_links[join_start];

        if left_join.is_some() && right_join.is_some() {
            add_single_join(
                histories,
                segments,
                segment_groups,
                prior_hist_idx,
                segment_idx,
                &join_block_indexes[join_start..join_end],
                left_join.as_ref(),
                right_join.as_ref(),
                score_params,
                history_depth,
            );
        }

        if left_join.is_some() {
            add_single_join(
                histories,
                segments,
                segment_groups,
                prior_hist_idx,
                segment_idx,
                &join_block_indexes[join_start..join_end],
                left_join.as_ref(),
                None,
                score_params,
                history_depth,
            );
        }

        if right_join.is_some() {
            add_single_join(
                histories,
                segments,
                segment_groups,
                prior_hist_idx,
                segment_idx,
                &join_block_indexes[join_start..join_end],
                None,
                right_join.as_ref(),
                score_params,
                history_depth,
            );
        }
    }
}

fn add_appends_to_history(
    histories: &mut Vec<HistoryEntry>,
    segments: &[Segment],
    segment_groups: &mut [SegmentGroups],
    segment_idx: usize,
    prior_hist_idx: usize,
    group_idx: usize,
    history_depth: usize,
    score_params: &ScoreParams,
) {
    let other_index = remove_expired_history_entries(
        &histories,
        &segment_groups,
        segment_idx,
        prior_hist_idx,
        history_depth,
    );

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
                prior_history: other_index, // We use 0 as a place holder (can't join to root, it's impossible), we replace this later once histories indexes are locked in for this step...
                join_left_block: JoinSide {
                    caused_by_history: JOIN_SIDE_SELF_PLACEHOLDER,
                    link_history: JOIN_SIDE_SELF_PLACEHOLDER,
                    side: Side::Left,
                },
                join_right_block: JoinSide {
                    caused_by_history: JOIN_SIDE_SELF_PLACEHOLDER,
                    link_history: JOIN_SIDE_SELF_PLACEHOLDER,
                    side: Side::Right,
                },
                score: new_score,
            }));
        }
        HistoryEntry::Append(val) | HistoryEntry::Join(val) => {
            // Can add up to two append events for blocks with multiple alignments:
            let current_rep_block =
                &segment_groups[segment_idx].get_first_block(&segments[segment_idx], group_idx);
            let prior_rep_block = &segment_groups[val.segment]
                .get_first_block(&segments[val.segment], val.group_index);
            let is_skip = current_rep_block.row_idx == 0 || prior_rep_block.row_idx == 0;

            let can_append = val.join_right_block.caused_by_history == prior_hist_idx
                && val.join_right_block.side == Side::Right;
            let current_group = segment_groups[segment_idx].get_group(group_idx);

            let (current_blocks, split_point) = if can_append {
                get_valid_appends_for_current_group(
                    &segments[segment_idx],
                    current_group,
                    &segments[val.segment],
                    segment_groups[val.segment].get_group(val.group_index),
                )
            } else {
                (current_group.to_vec(), 0)
            };

            let new_matching_group = &current_blocks[..split_point];
            let new_mismatching_group = &current_blocks[split_point..];

            // Notice that we don't add a history for adjacent blocks if it is an alignment.
            // This is because they are handled (correctly) by the join system...
            // So this is only needed for the skip state and tandem repeats.
            if !new_matching_group.is_empty()
                && !matches!(current_rep_block.block_type, BlockType::Alignment)
            {
                let new_group_idx = segment_groups[segment_idx]
                    .add_group(&segments[segment_idx], new_matching_group);
                let new_score = history_score(&histories[prior_hist_idx])
                    + score_params.transition(is_skip, false)
                    + current_rep_block.alignment_score;

                let (left_side, right_side) = get_join_endpoints_from_links(histories, None, None);

                // Add append event for mismatching blocks, this will have a transition penalty...
                histories.push(HistoryEntry::Append(HistoryInfo {
                    segment: segment_idx,
                    group_index: new_group_idx,
                    prior_block_history: prior_hist_idx,
                    prior_history: other_index,
                    join_left_block: left_side,
                    join_right_block: right_side,
                    score: new_score,
                }));
            }

            if !new_mismatching_group.is_empty() {
                let new_group_idx = segment_groups[segment_idx]
                    .add_group(&segments[segment_idx], new_mismatching_group);
                let new_score = history_score(&histories[prior_hist_idx])
                    + score_params.transition(is_skip, true)
                    + current_rep_block.alignment_score;

                let (left_side, right_side) = get_join_endpoints_from_links(histories, None, None);

                // Add append event for mismatching blocks, this will have a transition penalty...
                histories.push(HistoryEntry::Append(HistoryInfo {
                    segment: segment_idx,
                    group_index: new_group_idx,
                    prior_block_history: prior_hist_idx,
                    prior_history: other_index,
                    join_left_block: left_side,
                    join_right_block: right_side,
                    score: new_score,
                }));
            }
        }
    }
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
                let (join_links, join_block_indexes, join_split_points) =
                    check_for_joins(JoinCheckArgs {
                        histories: &histories,
                        segments,
                        segment_groups: &segment_groups,
                        assembly_graph,
                        current_group_reference: (segment_idx, group_idx),
                        start_entry: prior_hist_idx,
                        history_depth,
                        epsilon: 1e-2,
                    });

                add_joins_to_history(
                    &mut histories,
                    segments,
                    &mut segment_groups,
                    prior_hist_idx,
                    segment_idx,
                    &join_links,
                    &join_block_indexes,
                    &join_split_points,
                    score_params,
                    history_depth,
                );

                add_appends_to_history(
                    &mut histories,
                    segments,
                    &mut segment_groups,
                    segment_idx,
                    prior_hist_idx,
                    group_idx,
                    history_depth,
                    score_params,
                );
            }
        }

        if histories.len() == prior_step_end {
            panic!("No histories added by step logic!!!!")
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
        // Replace placeholder self-referencing history indexes (0) with actual index of the new histories...
        finalize_join_sides(&mut histories, prior_step_end);

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

#[derive(Debug)]
struct JoinStackEntry {
    joined_history_offset: usize,
    trace_segment_offset: usize,
    join_index: usize,
}

struct AddedBlockInfo {
    history_index: usize,
    trace_stack_index: usize,
    join_index: usize,
}

#[derive(Debug)]
struct JoinStack {
    pub stack: Vec<JoinStackEntry>,
    pub next_join_index: usize,
}

enum BlockAction {
    // join index, trace stack index...
    Join(usize, usize),
    // join index...
    Add(usize),
}

impl JoinStack {
    fn new() -> Self {
        JoinStack {
            stack: Vec::new(),
            next_join_index: 0,
        }
    }

    /// Try adding one or two joins to the join stack if this block is a join and has linked edges.
    fn try_push(&mut self, entry: &HistoryEntry, added_block_info: Option<&AddedBlockInfo>) {
        if let (HistoryEntry::Join(val), Some(info)) = (entry, added_block_info) {
            let top_offset = self.stack.len();

            if val.join_left_block.caused_by_history != info.history_index {
                self.stack.push(JoinStackEntry {
                    joined_history_offset: val.join_left_block.caused_by_history,
                    trace_segment_offset: info.trace_stack_index,
                    join_index: info.join_index,
                });
            }

            if val.join_right_block.caused_by_history != info.history_index {
                self.stack.push(JoinStackEntry {
                    joined_history_offset: val.join_right_block.caused_by_history,
                    trace_segment_offset: info.trace_stack_index,
                    join_index: info.join_index,
                });
            }

            self.stack[top_offset..].sort_unstable_by_key(|v| v.joined_history_offset);
        }
    }

    /// Attempt to join top-most join on the join stack if history index matches...
    fn check_for_join(&mut self, history_index: usize) -> BlockAction {
        if let Some(join_entry) = self.stack.last() {
            if join_entry.joined_history_offset == history_index {
                let block =
                    BlockAction::Join(join_entry.join_index, join_entry.trace_segment_offset);
                self.stack.pop();
                return block;
            }
        }

        let result = BlockAction::Add(self.next_join_index);
        self.next_join_index += 1;
        result
    }

    /// Same as add-block, but does not return information, simply just consumes the value on the top of the stack.
    fn try_pop(&mut self, history_index: usize) {
        if let Some(join_entry) = self.stack.last() {
            if join_entry.joined_history_offset == history_index {
                self.stack.pop();
            }
        }
    }
}

// Return is the assigned index of the new block, and the new join index to use for the following block...
fn history_backtrace_append_block(
    refined_segments: &mut Vec<RefinedTraceSegment>,
    joiner: &mut JoinStack,
    blocks: &[&Block],
    score: f64,
    current_history_index: usize,
    segment: usize,
    region_idx: usize,
) -> Option<AddedBlockInfo> {
    // Case 1: Same row index and touches start of segment in front of it, extend the segment backwards to include this...
    let stack_size = refined_segments.len();

    if let Some(ref_seg) = refined_segments.last_mut() {
        let direct_extensions =
            get_possible_extensions(blocks.iter().copied(), ref_seg.annotated.iter());

        if !direct_extensions.is_empty() {
            // If marked as join (happens for regular blocks), consume the value added to the join stack...
            joiner.try_pop(current_history_index);
            ref_seg.annotated = direct_extensions;
            return Some(AddedBlockInfo {
                history_index: current_history_index,
                trace_stack_index: stack_size - 1,
                join_index: ref_seg.join_index,
            });
        }
    }

    if blocks
        .iter()
        .any(|&b| matches!(b.block_type, BlockType::Alignment | BlockType::TandemRepeat))
    {
        match joiner.check_for_join(current_history_index) {
            // Case 2: Involved in a join, add new block.
            BlockAction::Join(join_index, stack_pos) => {
                let joins = get_joinable_extensions(
                    blocks.iter().copied(),
                    refined_segments[stack_pos].annotated.iter(),
                );

                // Should not be possible assuming a join was allowed in the first place...
                if joins.is_empty() {
                    panic!("Annotation from join made with 0 elements! This should not be possible! (Region: {}, Segment: {})", region_idx, segment);
                }

                refined_segments.push(RefinedTraceSegment {
                    annotated: joins,
                    join_index,
                    score,
                    segment,
                });

                Some(AddedBlockInfo {
                    history_index: current_history_index,
                    trace_stack_index: refined_segments.len() - 1,
                    join_index,
                })
            }
            BlockAction::Add(join_index) => {
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

                Some(AddedBlockInfo {
                    history_index: current_history_index,
                    trace_stack_index: refined_segments.len() - 1,
                    join_index,
                })
            }
        }
    } else {
        // Case 4: Skip state, don't add anything...=
        None
    }
}

pub fn backtrace_histories(
    segments: &SegmentedMatrix,
    history: &History,
    region_idx: usize,
) -> Vec<RefinedTraceSegment> {
    debug_assert!(segments.len() == history.segment_offsets.len() - 1);

    let mut refined_segments: Vec<RefinedTraceSegment> = Vec::new();
    // Both structures below used for tracking joins...
    let mut joiner: JoinStack = JoinStack::new();

    let last_segment = history.segment_offsets.len() - 1;
    // Find the max in the first row....
    let mut current_idx = history.segment_offsets[last_segment]
        + get_max_history(
            &history.entries[history.segment_offsets[last_segment]..],
            region_idx,
        );
    let mut current_entry = &history.entries[current_idx];

    while let HistoryEntry::Join(entry_info) | HistoryEntry::Append(entry_info) = current_entry {
        // Append current entry to segment stack...
        let blocks = history.segment_groups[entry_info.segment]
            .get_group(entry_info.group_index)
            .iter()
            .map(|&i| &segments[entry_info.segment].blocks[i])
            .collect_vec();

        // Append block for this entry (or extend prior trace block if this is the same alignment)...
        let added_block = history_backtrace_append_block(
            &mut refined_segments,
            &mut joiner,
            &blocks,
            entry_info.score,
            current_idx,
            entry_info.segment,
            region_idx,
        );

        // If this is a join, link segments it joins to so they can be constructed correctly later...
        joiner.try_push(current_entry, added_block.as_ref());

        // Go to the next entry in the history...
        current_idx = entry_info.prior_block_history;
        current_entry = &history.entries[current_idx];
    }

    if joiner.stack.len() != 0 {
        panic!(
            "Backtrace not done properly, there are {} leftover values on the join stack!",
            joiner.stack.len()
        );
    }

    // Reverse so trace segments go from start to end instead of end to start.
    refined_segments.reverse();
    refined_segments
}
