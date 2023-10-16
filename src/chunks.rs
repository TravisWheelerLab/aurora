use crate::alignment::{Alignment, AlignmentData};

#[derive(Copy, Clone)]
struct Interval {
    target_id: usize,
    target_start: usize,
    target_end: usize,
}

impl Interval {
    pub fn overlaps(&self, other: &Self) -> bool {
        self.target_start <= other.target_end && self.target_end >= other.target_start
    }
}

/// A group of alignments that are in
/// proximity to each other. The alignments
/// are implicitly in the same TargetGroup
pub struct ProximityGroup<'a> {
    pub target_id: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub alignments: &'a [Alignment],
}

impl<'a> std::fmt::Debug for ProximityGroup<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}-{}: {}",
            self.target_start,
            self.target_end,
            self.alignments.len()
        )
    }
}

impl<'a> ProximityGroup<'a> {
    pub fn from_alignment_data(
        alignment_data: &AlignmentData,
        join_distance: usize,
    ) -> Vec<ProximityGroup> {
        alignment_data.target_groups.iter().for_each(|g| {
            debug_assert!(g
                .alignments
                .windows(2)
                .all(|w| w[0].target_start <= w[1].target_start));
        });

        alignment_data
            .target_groups
            .iter()
            .flat_map(|target_group| {
                let query_ids = 0..alignment_data.query_name_map.size();
                let mut target_intervals: Vec<Interval> = query_ids
                    .skip(1)
                    .flat_map(|query_id| {
                        // grab all alignments to this query
                        let alignments: Vec<&Alignment> = target_group
                            .alignments
                            .iter()
                            .filter(|a| a.query_id == query_id)
                            .collect();

                        (0..alignments.len()).map(move |idx| {
                            let first_alignment = alignments[idx];
                            let other_alignments = &alignments[(idx + 1)..];

                            let mut interval = Interval {
                                target_id: target_group.target_id,
                                target_start: first_alignment.target_start,
                                target_end: first_alignment.target_end,
                            };

                            for next_alignment in other_alignments {
                                let distance = next_alignment
                                    .target_start
                                    .saturating_sub(first_alignment.target_end);

                                if distance > join_distance {
                                    break;
                                }

                                if first_alignment.strand == next_alignment.strand {
                                    interval.target_end =
                                        interval.target_end.max(next_alignment.target_end);
                                }
                            }

                            interval
                        })
                    })
                    .collect();

                target_intervals.sort_by_key(|interval| interval.target_start);

                let merged_intervals = target_intervals.iter().skip(1).fold(
                    vec![target_intervals[0]],
                    |mut acc, interval| {
                        let last = acc.last_mut().unwrap();
                        if last.overlaps(interval) {
                            last.target_end = last.target_end.max(interval.target_end)
                        } else {
                            acc.push(*interval);
                        }
                        acc
                    },
                );

                let mut start_idx = 0usize;
                merged_intervals.into_iter().map(move |interval| {
                    let end_idx = target_group
                        .alignments
                        .iter()
                        .enumerate()
                        .skip(start_idx)
                        .find(|(_, a)| a.target_start > interval.target_end)
                        .unwrap_or((target_group.alignments.len(), &Alignment::default()))
                        .0;

                    let alignments = &target_group.alignments[start_idx..end_idx];

                    start_idx = end_idx;
                    ProximityGroup {
                        target_id: interval.target_id,
                        target_start: interval.target_start,
                        target_end: interval.target_end,
                        alignments,
                    }
                })
            })
            .collect()
    }
}

pub fn validate_groups(groups: &[ProximityGroup], join_distance: usize) -> bool {
    for (group_idx, group) in groups.iter().enumerate() {
        for ali in group.alignments {
            // check if any alignments are outside of the chunk boundaries
            if ali.target_start < group.target_start || ali.target_end > group.target_end {
                return false;
            }

            // now check if we violate the join conditions

            // first check all chunks to the left
            for other_group_idx in 0..group_idx {
                let other_group = &groups[other_group_idx];

                // if the other chunk is farther than the join distance
                // away from this alignment, then there's nothing to check
                if (ali.target_start - other_group.target_end + 1) > join_distance {
                    continue;
                }

                for other_ali in other_group.alignments {
                    if other_ali.query_id != ali.query_id {
                        continue;
                    }

                    if ali.strand == other_ali.strand
                        && (ali.target_start - other_ali.target_end + 1) < join_distance
                    {
                        return false;
                    }
                }
            }

            // then check all chunks to the right
            for other_group_idx in (group_idx..groups.len()).skip(1) {
                let other_group = &groups[other_group_idx];

                if (other_group.target_start - ali.target_end + 1) > join_distance {
                    continue;
                }

                for other_ali in other_group.alignments {
                    if other_ali.query_id != ali.query_id {
                        continue;
                    }

                    if other_ali.strand == ali.strand
                        && (other_ali.target_start - ali.target_end + 1) < join_distance
                    {
                        return false;
                    }
                }
            }
        }
    }
    true
}
