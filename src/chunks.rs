use crate::alignment::{Alignment, AlignmentData, VecMap};

#[derive(Copy, Clone)]
struct Interval {
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
    pub target_start: usize,
    pub target_end: usize,
    pub alignments: &'a [Alignment],
    pub row_map: VecMap<usize>,
}

impl<'a> ProximityGroup<'a> {
    pub fn from_alignment_data(
        alignment_data: &mut AlignmentData,
        join_distance: usize,
    ) -> Vec<ProximityGroup> {
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
                            let ali = alignments[idx];
                            let other_ali = &alignments[(idx + 1)..];
                            let mut interval = Interval {
                                target_start: alignments[0].target_start,
                                target_end: alignments[0].target_end,
                            };
                            for next_ali in other_ali {
                                let distance = next_ali.target_start.saturating_sub(ali.target_end);

                                if distance > join_distance {
                                    break;
                                }

                                if ali.strand == next_ali.strand {
                                    interval.target_end = next_ali.target_end;
                                }
                            }

                            interval
                        })
                    })
                    .collect();

                target_intervals.sort_by_key(|interval| interval.target_start);

                let merged_intervals = target_intervals.iter().skip(1).fold(
                    vec![target_intervals[0].clone()],
                    |mut acc, interval| {
                        let last = acc.last_mut().unwrap();
                        if last.overlaps(&interval) {
                            last.target_end = interval.target_end
                        } else {
                            acc.push(interval.clone())
                        }
                        acc
                    },
                );

                merged_intervals.into_iter().map(|interval| {
                    let start_idx = target_group
                        .alignments
                        .iter()
                        .enumerate()
                        .find(|(_, a)| a.target_start == interval.target_start)
                        .expect("failed to match interval start")
                        .0;

                    let end_idx = target_group
                        .alignments
                        .iter()
                        .enumerate()
                        .rev()
                        .find(|(_, a)| a.target_end == interval.target_end)
                        .expect("failed to match interval end")
                        .0;

                    let alignments = &target_group.alignments[start_idx..=end_idx];
                    let mut row_map = VecMap::new();
                    alignments.iter().for_each(|a| {
                        row_map.insert(a.query_id);
                    });

                    ProximityGroup {
                        target_start: interval.target_start,
                        target_end: interval.target_end,
                        alignments,
                        row_map,
                    }
                })
            })
            .collect()
    }
}

//fn validate_chunks(chunks: &[Chunk], alignments: &[Alignment], join_distance: usize) -> bool {
//    for (chunk_idx, chunk) in chunks.iter().enumerate() {
//        for ali in &alignments[chunk.start_idx..=chunk.end_idx] {
//            // check if any alignments are outside of the chunk boundaries
//            if ali.target_start < chunk.target_start || ali.target_end > chunk.target_end {
//                return false;
//            }

//            // now check if we violate the join conditions

//            // first check all chunks to the left
//            for other_chunk_idx in 0..chunk_idx {
//                let other_chunk = &chunks[other_chunk_idx];

//                // if the other chunk is farther than the join distance
//                // away from this alignment, then there's nothing to check
//                if (ali.target_start - other_chunk.target_end + 1) > join_distance {
//                    continue;
//                }

//                for other_ali in &alignments[other_chunk.start_idx..=other_chunk.end_idx] {
//                    if other_ali.query_id != ali.query_id {
//                        continue;
//                    }

//                    if (ali.target_start - other_ali.target_end + 1) < join_distance {
//                        return false;
//                    }
//                }
//            }

//            // then check all chunks to the right
//            for other_chunk_idx in (chunk_idx..chunks.len()).skip(1) {
//                let other_chunk = &chunks[other_chunk_idx];

//                if (other_chunk.target_start - ali.target_end + 1) > join_distance {
//                    continue;
//                }

//                for other_ali in &alignments[other_chunk.start_idx..=other_chunk.end_idx] {
//                    if other_ali.query_id != ali.query_id {
//                        continue;
//                    }

//                    if (other_ali.target_start - ali.target_end + 1) < join_distance {
//                        return false;
//                    }
//                }
//            }
//        }
//    }
//    true
//}

//#[allow(dead_code)]
//pub fn chunk_memory_estimate(alignments: &[Alignment], chunks: &[Chunk]) {
//    let num_bytes_per_cell = 32usize;
//    let mut mem = vec![];

//    for chunk in chunks.iter() {
//        let mut target_start = usize::MAX;
//        let mut target_end = 0usize;
//        let mut num_cells_sparse = 0usize;

//        alignments[chunk.start_idx..=chunk.end_idx]
//            .iter()
//            .for_each(|a| {
//                target_start = target_start.min(a.target_start);
//                target_end = target_end.max(a.target_end);
//                num_cells_sparse += a.target_end - a.target_start + 1;
//            });

//        let target_width = target_end - target_start + 1;
//        let num_ali = chunk.end_idx - chunk.start_idx + 1;
//        num_cells_sparse += target_width;
//        let num_cells_full = target_width * num_ali + target_width;

//        mem.push((num_ali, num_cells_full, num_cells_sparse));
//    }

//    mem.sort();
//    mem.reverse();

//    mem.iter().for_each(|(a, f, s)| {
//        //
//        println!(
//            "{:10}, {:07.2}, {:0.2}",
//            a,
//            (f * num_bytes_per_cell) as f64 / 1e9,
//            (s * num_bytes_per_cell) as f64 / 1e9
//        )
//    });
//}
