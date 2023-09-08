use crate::alignment::Alignment;

#[derive(Copy, Clone)]
struct TargetSequenceChunk {
    target_start: usize,
    target_end: usize,
}

#[derive(Copy, Clone)]
struct AlignmentIndexChunk {
    start_idx: usize,
    end_idx: usize,
}

pub struct Chunk {
    /// The index of the first alignment in the chunk
    pub start_idx: usize,
    /// The index of the last alignment in the chunk
    pub end_idx: usize,
    /// The start position of the chunk on the target sequence.
    ///
    /// Note: this is the same as the `target_start` position of the first alignment in the chunk
    /// because we know there are no alignments that start earlier than the first one
    pub target_start: usize,
    /// The end position of the chunk on the target sequence
    ///
    /// Note: this is NOT the same as the `target_end` position of the last alignment in the chunk
    /// because there may be an earlier alignment that extends past the end of the last one
    pub target_end: usize,
}

///
///
///
///
///
pub fn chunks_by_query_distance(
    alignments: &mut [Alignment],
    num_queries: usize,
    join_distance: usize,
) -> Vec<Chunk> {
    let num_ali = alignments.len();

    alignments.sort_by(|a, b| {
        a.query_id
            .cmp(&b.query_id)
            .then_with(|| a.target_start.cmp(&b.target_start))
    });

    // first we'll chunk the alignments up by the query
    let mut chunks_by_query_id: Vec<AlignmentIndexChunk> = vec![
        AlignmentIndexChunk {
            start_idx: 0,
            end_idx: 0,
        };
        num_queries
    ];

    let mut prev_query_id = 0usize;
    let mut current_ali_start_idx = 0usize;
    for (idx, ali) in alignments.iter().enumerate() {
        let query_id = ali.query_id;
        if query_id != prev_query_id {
            chunks_by_query_id[query_id - 1].start_idx = current_ali_start_idx;
            chunks_by_query_id[query_id - 1].end_idx = idx.saturating_sub(1);
            current_ali_start_idx = idx;
            prev_query_id = query_id;
        }
    }

    let last_chunk_by_query_id = chunks_by_query_id
        .last_mut()
        .expect("chunks_by_query_id is empty");

    last_chunk_by_query_id.start_idx = current_ali_start_idx;
    last_chunk_by_query_id.end_idx = num_ali - 1;

    // now that we have a nice grouping of alignments with the
    // same query, we can process each slice of alignments to
    // figure out which groups are close enough to join
    let mut chunks_by_target_distance: Vec<TargetSequenceChunk> = vec![];

    for index_chunk in chunks_by_query_id.iter() {
        let ali_slice = &alignments[index_chunk.start_idx..=index_chunk.end_idx];

        let first_alignment = ali_slice.first().expect("alignments are empty");

        let mut current_target_start = first_alignment.target_start;
        let mut current_target_end = first_alignment.target_end;

        for ali in ali_slice.iter() {
            let distance = ali.target_start.saturating_sub(current_target_end + 1);
            if distance > join_distance {
                chunks_by_target_distance.push(TargetSequenceChunk {
                    target_start: current_target_start,
                    target_end: current_target_end,
                });

                current_target_start = ali.target_start;
                current_target_end = ali.target_end;
            } else {
                current_target_end = current_target_end.max(ali.target_end);
            }
        }

        chunks_by_target_distance.push(TargetSequenceChunk {
            target_start: current_target_start,
            target_end: current_target_end,
        });
    }

    // now we need to look at all alignments, in order, regardless
    // of their query, to identify anything that overlaps at all
    alignments.sort_by_key(|a| a.target_start);

    let first_alignment = alignments.first().expect("alignments are empty");

    let mut current_target_start = first_alignment.target_start;
    let mut current_target_end = first_alignment.target_end;

    for ali in alignments.iter() {
        if current_target_end >= ali.target_start {
            current_target_end = ali.target_end;
        } else {
            chunks_by_target_distance.push(TargetSequenceChunk {
                target_start: current_target_start,
                target_end: current_target_end,
            });
            current_target_start = ali.target_start;
            current_target_end = ali.target_end;
        }
    }
    chunks_by_target_distance.push(TargetSequenceChunk {
        target_start: current_target_start,
        target_end: current_target_end,
    });

    // now we are going to merge all of the overlapping chunks
    chunks_by_target_distance.sort_by_key(|c| c.target_start);
    let mut merged_chunks = vec![];

    let mut current_chunk = chunks_by_target_distance[0];
    for chunk in chunks_by_target_distance.iter().skip(1) {
        if current_chunk.target_end >= chunk.target_start {
            current_chunk.target_end = current_chunk.target_end.max(chunk.target_end);
        } else {
            merged_chunks.push(current_chunk);
            current_chunk = *chunk;
        }
    }
    merged_chunks.push(current_chunk);

    // figure out which alignments are in the chunks
    let mut final_chunks: Vec<Chunk> = vec![];

    let mut chunk_iter = merged_chunks.iter();
    let mut current_chunk = chunk_iter.next().expect("merged_chunks is empty");
    let mut current_ali_start_idx = 0usize;

    for (idx, ali) in alignments.iter().enumerate() {
        if ali.target_start > current_chunk.target_end {
            final_chunks.push(Chunk {
                start_idx: current_ali_start_idx,
                end_idx: idx - 1,
                target_start: current_chunk.target_start,
                target_end: current_chunk.target_end,
            });
            current_ali_start_idx = idx;
            current_chunk = chunk_iter.next().unwrap();
        }
    }

    final_chunks.push(Chunk {
        start_idx: current_ali_start_idx,
        end_idx: num_ali - 1,
        target_start: current_chunk.target_start,
        target_end: current_chunk.target_end,
    });

    debug_assert!(validate_chunks(&final_chunks, alignments, join_distance));
    final_chunks
}

///
///
///
///
///
fn validate_chunks(chunks: &[Chunk], alignments: &[Alignment], join_distance: usize) -> bool {
    for (chunk_idx, chunk) in chunks.iter().enumerate() {
        for ali in &alignments[chunk.start_idx..=chunk.end_idx] {
            // check if any alignments are outside of the chunk boundaries
            if ali.target_start < chunk.target_start || ali.target_end > chunk.target_end {
                return false;
            }

            // now check if we violate the join conditions

            // first check all chunks to the left
            for other_chunk_idx in 0..chunk_idx {
                let other_chunk = &chunks[other_chunk_idx];

                // if the other chunk is farther than the join distance
                // away from this alignment, then there's nothing to check
                if (ali.target_start - other_chunk.target_end + 1) > join_distance {
                    continue;
                }

                for other_ali in &alignments[other_chunk.start_idx..=other_chunk.end_idx] {
                    if other_ali.query_id != ali.query_id {
                        continue;
                    }

                    if (ali.target_start - other_ali.target_end + 1) < join_distance {
                        return false;
                    }
                }
            }

            // then check all chunks to the right
            for other_chunk_idx in (chunk_idx..chunks.len()).skip(1) {
                let other_chunk = &chunks[other_chunk_idx];

                if (other_chunk.target_start - ali.target_end + 1) > join_distance {
                    continue;
                }

                for other_ali in &alignments[other_chunk.start_idx..=other_chunk.end_idx] {
                    if other_ali.query_id != ali.query_id {
                        continue;
                    }

                    if (other_ali.target_start - ali.target_end + 1) < join_distance {
                        return false;
                    }
                }
            }
        }
    }
    true
}

#[allow(dead_code)]
pub fn chunk_memory_estimate(alignments: &[Alignment], chunks: &[Chunk]) {
    let num_bytes_per_cell = 32usize;
    let mut mem = vec![];

    for chunk in chunks.iter() {
        let mut target_start = usize::MAX;
        let mut target_end = 0usize;
        let mut num_cells_sparse = 0usize;

        alignments[chunk.start_idx..=chunk.end_idx]
            .iter()
            .for_each(|a| {
                target_start = target_start.min(a.target_start);
                target_end = target_end.max(a.target_end);
                num_cells_sparse += a.target_end - a.target_start + 1;
            });

        let target_width = target_end - target_start + 1;
        let num_ali = chunk.end_idx - chunk.start_idx + 1;
        num_cells_sparse += target_width;
        let num_cells_full = target_width * num_ali + target_width;

        mem.push((num_ali, num_cells_full, num_cells_sparse));
    }

    mem.sort();
    mem.reverse();

    mem.iter().for_each(|(a, f, s)| {
        //
        println!(
            "{:10}, {:07.2}, {:0.2}",
            a,
            (f * num_bytes_per_cell) as f64 / 1e9,
            (s * num_bytes_per_cell) as f64 / 1e9
        )
    });
}
