use crate::balanced_tree::{AVLIndexSet, SetInsert};
use crate::segments::{Block, Segment};
use itertools::Itertools;
use std::cmp::Ordering;

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
                    f64::total_cmp(
                        &segment.blocks[a].alignment_score,
                        &segment.blocks[b].alignment_score,
                    )
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
                    || (segment.blocks[score_ordered[next_index]].alignment_score
                        - segment.blocks[score_ordered[last_index]].alignment_score)
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
