use crate::{
    history_tracing::{history_score, History, RefinedTraceSegment},
    segments::SegmentedMatrix,
};
use itertools::{izip, Itertools};
use std::{
    fs::File,
    io::{self, Write},
    path::Path,
};

pub fn dump_final_trace_statistics(
    history: &History,
    segments: &SegmentedMatrix,
    trace_segments: &[RefinedTraceSegment],
    path: impl AsRef<Path>,
) -> io::Result<()> {
    let mut file = File::create(path)?;

    writeln!(
        &mut file,
        "Segment, History Count, Trace Score Rank, Relative Trace Score, Absolute Trace Score, Computed Absolute Bound, Computed Relative Bound"
    )?;

    for seg in trace_segments.iter() {
        let i = seg.segment;
        let start_off = history.segment_offsets[i + 1];
        let end_off = if i + 2 >= history.segment_offsets.len() {
            history.entries.len()
        } else {
            history.segment_offsets[i + 2]
        };

        let mut sorted_scores = history.entries[start_off..end_off]
            .iter()
            .map(history_score)
            .collect_vec();
        sorted_scores.sort_by(f64::total_cmp);
        let index = match sorted_scores.binary_search_by(|v| v.total_cmp(&seg.score)) {
            Result::Ok(index) | Result::Err(index) => index,
        };
        let rank = sorted_scores.len() - index;
        let score_below_best = seg.score - sorted_scores.last().unwrap_or(&0.0);

        writeln!(
            &mut file,
            "{}, {}, {}, {}, {}, {}, {}",
            i,
            end_off - start_off,
            rank,
            score_below_best,
            seg.score,
            segments[i].absolute_score_bound,
            segments[i].relative_score_bound
        )?;
    }

    Ok(())
}

pub fn dump_debug_history_info(
    history: &History,
    segments: &SegmentedMatrix,
    target_start: usize,
    history_lengths: &[usize],
    path: impl AsRef<Path>,
) -> io::Result<()> {
    let mut file = File::create(path)?;

    let segment_lengths = segments.iter().map(|s| s.blocks.len());

    let segment_ranges = segments
        .iter()
        .map(|s| (target_start + s.start_col, target_start + s.end_col));
    let num_groups = history.segment_groups.iter().map(|s| s.group_count());
    let group_sizes = history.segment_groups.iter().map(|s| s.index_count());
    let segment_score_bounds = segments
        .iter()
        .map(|v| (v.absolute_score_bound, v.relative_score_bound));

    let max_history = (0..history.segment_offsets.len()).map(|i| {
        let start = history.segment_offsets[i];
        let end = if i + 1 < history.segment_offsets.len() {
            history.segment_offsets[i + 1]
        } else {
            history.entries.len()
        };

        history.entries[start..end]
            .iter()
            .map(history_score)
            .max_by(f64::total_cmp)
            .unwrap_or(f64::NEG_INFINITY)
    });

    writeln!(
        &mut file,
        "Segment, History Count, Group Count, Index Count, Block Count, Target Start, Target End, Computed Absolute Bound, Computed Relative Bound, Max History"
    )?;

    izip!(
        (0..segments.len()),
        history_lengths.iter(),
        num_groups,
        group_sizes,
        segment_lengths,
        segment_ranges,
        segment_score_bounds,
        max_history
    )
    .try_for_each(|v| {
        writeln!(
            &mut file,
            "{}, {}, {}, {}, {}, {}, {}, {}, {}, {}",
            v.0, v.1, v.2, v.3, v.4, v.5 .0, v.5 .1, v.6 .0, v.6 .1, v.7
        )
    })?;

    Ok(())
}
