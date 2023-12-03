use std::collections::HashMap;

use itertools::Itertools;

use crate::{
    collapse::{Assembly, AssemblyGroup},
    viterbi::TraceSegment,
    Args,
};

#[derive(Debug)]
struct MatrixRange {
    pub col_start: usize,
    pub col_end: usize,
}

struct AssemblyStuff<'a> {
    assembly: &'a Assembly<'a>,
    col_start: usize,
    col_end: usize,
    row_idx: usize,
    alignment_matrix_ranges: Vec<MatrixRange>,
    gap_matrix_ranges: Vec<MatrixRange>,
}

pub fn split_trace(
    trace_segments: Vec<TraceSegment>,
    assembly_group: &AssemblyGroup,
    active_cols: &Vec<usize>,
    confidence_avg_by_id: &HashMap<usize, f64>,
    args: &Args,
) -> (
    Vec<TraceSegment>,
    Vec<TraceSegment>,
    Vec<usize>,
    Vec<usize>,
    Vec<usize>,
) {
    let mut rows_in_trace = trace_segments
        .iter()
        .map(|s| s.row_idx)
        // we don't want to include
        // the skip state trace rows
        .filter(|&r| r != 0)
        .unique()
        .collect_vec();

    // figure out which assemblies have a
    // section that was hit by the trace
    let assemblies_in_trace = assembly_group
        .assemblies
        .iter()
        .enumerate()
        // +1 to the assembly index to get the row index
        .map(|(idx, assembly)| (idx + 1, assembly))
        // include assemblies that were hit by the trace
        .filter(|(row_idx, _)| rows_in_trace.contains(row_idx))
        // include assemblies with more than one alignment
        .filter(|(_, assembly)| assembly.alignments.len() > 1)
        .collect_vec();

    let stuff: Vec<AssemblyStuff> = assemblies_in_trace
        .iter()
        .map(|(row_idx, assembly)| {
            let assembly_col_start = assembly.target_start - assembly_group.target_start;
            let assembly_col_end = assembly.target_end - assembly_group.target_start;

            AssemblyStuff {
                assembly,
                col_start: assembly_col_start,
                col_end: assembly_col_end,
                row_idx: *row_idx,
                alignment_matrix_ranges: assembly
                    .alignment_ranges
                    .iter()
                    .map(|a| MatrixRange {
                        col_start: a.assembly_col_start + assembly_col_start,
                        col_end: a.assembly_col_end + assembly_col_start,
                    })
                    .collect_vec(),
                gap_matrix_ranges: assembly
                    .alignment_ranges
                    .iter()
                    .zip(assembly.alignment_ranges.iter().skip(1))
                    .map(|(a, b)| MatrixRange {
                        col_start: a.assembly_col_end + assembly_col_start,
                        col_end: b.assembly_col_start + assembly_col_start,
                    })
                    .collect_vec(),
            }
        })
        .collect_vec();

    // determine assembly conflicts
    let mut conflicts = vec![];
    stuff.iter().enumerate().for_each(|(idx_a, stuff_a)| {
        stuff
            .iter()
            .enumerate()
            .skip(idx_a + 1)
            .for_each(|(idx_b, stuff_b)| {
                // first check if B happens to be properly contained in a gap of A
                let b_in_a = stuff_a.gap_matrix_ranges.iter().any(|gap| {
                    stuff_b.col_start >= gap.col_start - args.fudge_distance
                        && stuff_b.col_end <= gap.col_end + args.fudge_distance
                });

                if b_in_a {
                    return;
                }

                // now check if A happens to be properly contained in a gap of B
                let a_in_b = stuff_b.gap_matrix_ranges.iter().any(|gap| {
                    stuff_a.col_start + args.fudge_distance >= gap.col_start - args.fudge_distance
                        && stuff_a.col_end - args.fudge_distance
                            <= gap.col_end + args.fudge_distance
                });

                if a_in_b {
                    return;
                }

                // now check if any fragment of B is inside a gap in A
                let conflict_a = stuff_a.gap_matrix_ranges.iter().any(|gap| {
                    stuff_b.alignment_matrix_ranges.iter().any(|frag| {
                        frag.col_start < gap.col_end - args.fudge_distance
                            && frag.col_end > gap.col_start + args.fudge_distance
                    })
                });

                if conflict_a {
                    conflicts.push((idx_a, idx_b));
                    return;
                }

                // now check if any fragment of A is inside a gap in B
                let conflict_b = stuff_b.gap_matrix_ranges.iter().any(|gap| {
                    stuff_a.alignment_matrix_ranges.iter().any(|frag| {
                        frag.col_start < gap.col_end - args.fudge_distance
                            && frag.col_end > gap.col_start + args.fudge_distance
                    })
                });

                if conflict_b {
                    conflicts.push((idx_a, idx_b));
                }
            });
    });

    let mut remaining_conflicts: HashMap<usize, bool> = HashMap::new();
    conflicts.iter().for_each(|&(a, b)| {
        remaining_conflicts.insert(a, true);
        remaining_conflicts.insert(b, true);
    });

    let mut competed_assembly_rows: Vec<usize> = vec![];

    conflicts.iter().for_each(|(a, b)| {
        if *remaining_conflicts.get(a).unwrap() && *remaining_conflicts.get(b).unwrap() {
            let ali_a = &stuff[*a].assembly.alignments;
            let ali_b = &stuff[*b].assembly.alignments;

            let conf_a: f64 = ali_a
                .iter()
                .map(|a| confidence_avg_by_id.get(&a.id).unwrap())
                .sum::<f64>()
                / ali_a.len() as f64;

            let conf_b: f64 = ali_b
                .iter()
                .map(|a| confidence_avg_by_id.get(&a.id).unwrap())
                .sum::<f64>()
                / ali_b.len() as f64;

            if conf_a >= conf_b {
                competed_assembly_rows.push(stuff[*b].row_idx);
                remaining_conflicts.insert(*b, false);
            } else {
                competed_assembly_rows.push(stuff[*a].row_idx);
                remaining_conflicts.insert(*a, false);
            }
        }
    });

    let winner_assemblies = stuff
        .into_iter()
        .filter(|s| !competed_assembly_rows.contains(&s.row_idx))
        .collect_vec();

    // convert the active columns into inactive ranges so
    // they are easier to compare against assembly gaps
    let mut inactive_col_ranges: Vec<MatrixRange> = vec![];
    active_cols
        .iter()
        .zip(active_cols.iter().skip(1))
        .for_each(|(&a, &b)| {
            if b - 1 != a {
                inactive_col_ranges.push(MatrixRange {
                    col_start: a + 1,
                    col_end: b - 1,
                });
            }
        });

    let (resolved_assemblies, unresolved_assemblies): (Vec<AssemblyStuff>, Vec<AssemblyStuff>) =
        winner_assemblies.into_iter().partition(|stuff| {
            // if for ALL of the gaps in the assembly
            stuff.gap_matrix_ranges.iter().all(|assembly_gap| {
                // if the gap is smaller than the fudge distance
                if assembly_gap.col_end - assembly_gap.col_start <= args.fudge_distance {
                    return true;
                }
                // or if the gap is contained in ANY inactive column range
                inactive_col_ranges.iter().any(|range| {
                    assembly_gap.col_start >= range.col_start - args.fudge_distance
                        && assembly_gap.col_end <= range.col_end + args.fudge_distance
                })
            })
        });

    let resolved_assembly_rows = resolved_assemblies.iter().map(|s| s.row_idx).collect_vec();

    let unresolved_assembly_rows = unresolved_assemblies
        .iter()
        .map(|s| s.row_idx)
        .collect_vec();

    let mut ambiguous = vec![];
    let mut conclusive = vec![];

    trace_segments
        .into_iter()
        .filter(|s| s.row_idx != 0)
        .for_each(|segment| {
            if unresolved_assembly_rows.contains(&segment.row_idx) {
                ambiguous.push(segment);
                return;
            }

            let overlapping_unresolved_assembly_parts = unresolved_assemblies
                .iter()
                .filter(|a| a.row_idx != segment.row_idx)
                .flat_map(|a| &a.alignment_matrix_ranges)
                .filter(|f| segment.col_start < f.col_end && segment.col_end > f.col_start)
                .collect_vec();

            if overlapping_unresolved_assembly_parts.is_empty() {
                conclusive.push(segment);
            } else {
                let min_overlap = overlapping_unresolved_assembly_parts
                    .iter()
                    .map(|f| f.col_start)
                    .min()
                    .unwrap_or(segment.col_start);

                let max_overlap = overlapping_unresolved_assembly_parts
                    .iter()
                    .map(|f| f.col_end)
                    .max()
                    .unwrap_or(segment.col_end);

                let slice_start = min_overlap.max(segment.col_start);
                let slice_end = max_overlap.min(segment.col_end);

                if slice_start > segment.col_start {
                    conclusive.push(TraceSegment {
                        query_id: segment.query_id,
                        row_idx: segment.row_idx,
                        col_start: segment.col_start,
                        col_end: slice_start - 1,
                        confidence: segment.confidence,
                    })
                }

                if slice_end < segment.col_end {
                    conclusive.push(TraceSegment {
                        query_id: segment.query_id,
                        row_idx: segment.row_idx,
                        col_start: slice_end + 1,
                        col_end: segment.col_end,
                        confidence: segment.confidence,
                    })
                }

                ambiguous.push(TraceSegment {
                    query_id: segment.query_id,
                    row_idx: segment.row_idx,
                    col_start: slice_start,
                    col_end: slice_end,
                    confidence: segment.confidence,
                });
            }
        });

    (
        ambiguous,
        conclusive,
        resolved_assembly_rows,
        unresolved_assembly_rows,
        competed_assembly_rows,
    )
}
