use std::{collections::HashMap, fs};

use itertools::Itertools;

use crate::{
    alignment::AlignmentData,
    alphabet::UTF8_TO_DIGITAL_NUCLEOTIDE,
    chunks::ProximityGroup,
    collapse::{Assembly, AssemblyGroup},
    confidence::confidence,
    matrix::{Matrix, MatrixDef},
    results::Annotation,
    score_params::ScoreParams,
    segments::Segments,
    support::windowed_confidence_slow,
    viterbi::{trace_segments, traceback, traceback2, viterbi, viterbi_collapsed, TraceSegment},
    viz::{write_soda_html, AdjudicationSodaData, AdjudicationSodaData2},
    windowed_scores::windowed_score,
    Args, BACKGROUND_WINDOW_SIZE, SCORE_WINDOW_SIZE,
};

pub fn run_assembly_pipeline(
    group: &ProximityGroup,
    alignment_data: &AlignmentData,
    region_idx: usize,
    mut args: Args,
) {
    if args.viz {
        args.viz_output_path.push(format!(
            "{}-{}-{}-{}",
            group.target_id, region_idx, group.target_start, group.target_end
        ));

        fs::create_dir_all(&args.viz_output_path).unwrap();
    }

    let score_params = ScoreParams::new(
        group.alignments.len(),
        args.query_jump_probability,
        args.num_skip_loops_eq_to_jump,
    );

    let matrix_def = MatrixDef::new(group);
    let mut confidence_matrix = Matrix::<f64>::new(&matrix_def);

    windowed_score(
        &mut confidence_matrix,
        group.alignments,
        &alignment_data.substitution_matrices,
        SCORE_WINDOW_SIZE,
        BACKGROUND_WINDOW_SIZE,
    );

    // TODO: probably remove this skip confidence_by_col thingy
    let _ = confidence(&mut confidence_matrix);

    let (confidence_avg_by_id, confidence_by_id) = windowed_confidence_slow(&mut confidence_matrix);

    // adjust the skip state to include skip-loop penalty
    let skip_adjust = args
        .query_jump_probability
        .powf(1.0 / args.num_skip_loops_eq_to_jump as f64);

    (0..confidence_matrix.num_cols()).for_each(|col_idx| {
        confidence_matrix.set_skip(col_idx, confidence_matrix.get_skip(col_idx) * skip_adjust);
    });

    let assembly_group = AssemblyGroup::new(group, &confidence_avg_by_id, &confidence_by_id, &args);

    let collapsed_matrix_def = MatrixDef::from_assembly_group(&assembly_group);

    let mut collapsed_confidence_matrix = Matrix::<f64>::new(&collapsed_matrix_def);
    let mut viterbi_matrix = Matrix::<f64>::new(&collapsed_matrix_def);
    let mut sources_matrix = Matrix::<usize>::new(&collapsed_matrix_def);

    collapsed_confidence_matrix.copy_fill(&confidence_matrix);

    // the initial active cols just removes
    // the dead space between alignments
    let mut active_cols = collapsed_confidence_matrix.initial_active_cols();
    let mut join_id_cnt = 0usize;
    let mut relevant_rows: Vec<Vec<usize>> = vec![];
    let mut trace_conclusive: Vec<Vec<TraceSegment>> = vec![];
    let mut trace_ambiguous: Vec<Vec<TraceSegment>> = vec![];

    let mut STOP = false;

    while !active_cols.is_empty() {
        viterbi_collapsed(
            &collapsed_confidence_matrix,
            &mut viterbi_matrix,
            &mut sources_matrix,
            &active_cols,
            &score_params,
        );

        let trace = traceback2(
            &viterbi_matrix,
            &collapsed_confidence_matrix,
            &sources_matrix,
            &active_cols,
            join_id_cnt,
        );

        debug_assert_eq!(trace.len(), active_cols.len());

        let trace_segments = trace_segments(&trace);

        join_id_cnt = trace
            .iter()
            .map(|t| t.join_id)
            .max()
            .expect("trace is empty");

        // -----------
        // row culling
        // -----------
        // note: this is currently just for visualization purposes.

        let mut relevant = vec![];

        assembly_group
            .assemblies
            .iter()
            .enumerate()
            .filter(|(_, a)| a.alignments.len() > 1)
            .for_each(|(idx, _)| relevant.push(idx + 1));

        trace_segments.iter().for_each(|s| relevant.push(s.row_idx));
        relevant = relevant.into_iter().unique().collect();
        relevant_rows.push(relevant);

        // ---------------
        // trace splitting
        // ---------------

        let mut rows_in_trace = trace_segments
            .iter()
            .map(|s| s.row_idx)
            .unique()
            .collect_vec();

        // we don't want to include the skip state trace rows
        rows_in_trace.retain(|&r| r != 0);

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
            fragments: Vec<MatrixRange>,
            gaps: Vec<MatrixRange>,
        }

        let stuff: Vec<AssemblyStuff> = assemblies_in_trace
            .iter()
            .map(|(row_idx, assembly)| {
                let assembly_col_start = assembly.target_start - group.target_start;
                let assembly_col_end = assembly.target_end - group.target_start;

                AssemblyStuff {
                    assembly,
                    col_start: assembly_col_start,
                    col_end: assembly_col_end,
                    row_idx: *row_idx,
                    fragments: assembly
                        .alignment_ranges
                        .iter()
                        .map(|a| MatrixRange {
                            col_start: a.assembly_col_start + assembly_col_start,
                            col_end: a.assembly_col_end + assembly_col_start,
                        })
                        .collect_vec(),
                    gaps: assembly
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
                    let b_in_a = stuff_a.gaps.iter().any(|gap| {
                        stuff_b.col_start >= gap.col_start - args.fudge_distance
                            && stuff_b.col_end <= gap.col_end + args.fudge_distance
                    });

                    if b_in_a {
                        return;
                    }

                    // now check if A happens to be properly contained in a gap of B
                    let a_in_b = stuff_b.gaps.iter().any(|gap| {
                        stuff_a.col_start + args.fudge_distance
                            >= gap.col_start - args.fudge_distance
                            && stuff_a.col_end - args.fudge_distance
                                <= gap.col_end + args.fudge_distance
                    });

                    if a_in_b {
                        return;
                    }

                    // now check if any fragment of B is inside a gap in A
                    let conflict_a = stuff_a.gaps.iter().any(|gap| {
                        stuff_b.fragments.iter().any(|frag| {
                            frag.col_start < gap.col_end - args.fudge_distance
                                && frag.col_end > gap.col_start + args.fudge_distance
                        })
                    });

                    if conflict_a {
                        conflicts.push((idx_a, idx_b));
                        return;
                    }

                    // now check if any fragment of A is inside a gap in B
                    let conflict_b = stuff_b.gaps.iter().any(|gap| {
                        stuff_a.fragments.iter().any(|frag| {
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

        let mut unassembled_rows: Vec<usize> = vec![];

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
                    unassembled_rows.push(stuff[*b].row_idx);
                    remaining_conflicts.insert(*b, false);
                } else {
                    unassembled_rows.push(stuff[*a].row_idx);
                    remaining_conflicts.insert(*a, false);
                }
            }
        });

        let assembled_stuff = stuff
            .into_iter()
            .filter(|s| !unassembled_rows.contains(&s.row_idx))
            .collect_vec();

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
            assembled_stuff.into_iter().partition(|stuff| {
                // if for ALL of the gaps in the assembly
                stuff.gaps.iter().all(|assembly_gap| {
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

        let rows = rows_in_trace
            .into_iter()
            .filter(|row_idx| {
                !unresolved_assemblies
                    .iter()
                    .map(|s| s.row_idx)
                    .contains(row_idx)
            })
            .collect_vec();

        let mut ambiguous = vec![];
        let mut conclusive = vec![];

        trace_segments
            .into_iter()
            .filter(|s| s.row_idx != 0)
            .for_each(|segment| {
                if !rows.contains(&segment.row_idx) {
                    ambiguous.push(segment);
                    return;
                }

                let overlapping_fragments = unresolved_assemblies
                    .iter()
                    .filter(|a| a.row_idx != segment.row_idx)
                    .flat_map(|a| &a.fragments)
                    .filter(|f| segment.col_start < f.col_end && segment.col_end > f.col_start)
                    .collect_vec();

                if overlapping_fragments.is_empty() {
                    conclusive.push(segment);
                } else {
                    let min_overlap = overlapping_fragments
                        .iter()
                        .map(|f| f.col_start)
                        .min()
                        .unwrap_or(segment.col_start);

                    let max_overlap = overlapping_fragments
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

        // ------------------
        // update active cols
        // ------------------

        let new_active_cols = ambiguous
            .iter()
            .flat_map(|s| s.col_start..=s.col_end)
            .collect_vec();

        debug_assert!(new_active_cols.len() <= active_cols.len());

        if new_active_cols.len() == active_cols.len() {
            break;
        }

        active_cols = new_active_cols;

        trace_conclusive.push(conclusive);
        trace_ambiguous.push(ambiguous);
    }

    let annotations: Vec<Annotation> = trace_conclusive
        .iter()
        .flat_map(|iter_segments| iter_segments.iter().map(Annotation::new))
        .collect();

    if args.viz {
        let s1 = trace_conclusive[0]
            .iter()
            .map(|seg| {
                format!(
                    "{},{},{},{}",
                    seg.col_start, seg.col_end, seg.query_id, seg.row_idx
                )
            })
            .join("|");

        let s2 = trace_ambiguous[0]
            .iter()
            .map(|seg| {
                format!(
                    "{},{},{},{}",
                    seg.col_start, seg.col_end, seg.query_id, seg.row_idx
                )
            })
            .join("|");

        let data = AdjudicationSodaData2::new(
            &assembly_group,
            &alignment_data.query_name_map,
            &annotations,
            &vec![],
            vec![s1, s2],
            &args,
        );

        let out_path = args.viz_output_path.join(format!(
            "{}-{}-{}.html",
            alignment_data.target_name_map.get(group.target_id),
            matrix_def.target_start,
            matrix_def.target_start + matrix_def.num_cols
        ));

        write_soda_html(
            &data,
            "./fixtures/soda/annotations.html",
            "./fixtures/soda/annotations.js",
            out_path,
        );
    }

    if STOP {
        panic!();
    }
}

pub fn run_pipeline(
    group: &ProximityGroup,
    alignment_data: &AlignmentData,
    region_idx: usize,
    args: &Args,
) {
    let score_params = ScoreParams::new(
        group.alignments.len(),
        args.query_jump_probability,
        args.num_skip_loops_eq_to_jump,
    );

    let matrix_def = MatrixDef::new(group);
    let mut confidence_matrix = Matrix::<f64>::new(&matrix_def);
    let mut viterbi_matrix = Matrix::<f64>::new(&matrix_def);
    let mut sources_matrix = Matrix::<usize>::new(&matrix_def);

    windowed_score(
        &mut confidence_matrix,
        group.alignments,
        &alignment_data.substitution_matrices,
        SCORE_WINDOW_SIZE,
        BACKGROUND_WINDOW_SIZE,
    );

    confidence(&mut confidence_matrix);

    windowed_confidence_slow(&mut confidence_matrix);

    // this removes any initial dead space between alignments

    // let mut active_cols = confidence_matrix.initial_active_cols();
    let mut active_cols = (0..confidence_matrix.num_cols()).collect::<Vec<_>>();

    let mut last_num_cols = active_cols.len();
    let mut results: Vec<Annotation> = vec![];

    let mut trace_strings = vec![];
    let mut fragment_strings = vec![];
    let mut segment_strings = vec![];

    let mut join_id_cnt = 0usize;

    while !active_cols.is_empty() {
        viterbi(
            &confidence_matrix,
            &mut viterbi_matrix,
            &mut sources_matrix,
            &active_cols,
            &score_params,
            alignment_data.query_name_map.size(),
            args.consensus_join_distance,
        );

        // let data = hyper_traceback(&sources_matrix, group.alignments, alignment_data);
        // write_soda_html(
        //     &data,
        //     "./fixtures/soda/trace-template.html",
        //     "./fixtures/soda/trace.js",
        //     "./trace.html",
        // )n;

        let trace = traceback(&viterbi_matrix, &sources_matrix, &active_cols, join_id_cnt);
        join_id_cnt = trace
            .iter()
            .map(|t| t.join_id)
            .max()
            .expect("trace is empty");

        let mut segments = Segments::new(&trace, &confidence_matrix, &active_cols);
        segments.create_links(
            args.consensus_join_distance,
            args.target_join_distance,
            args.num_skip_loops_eq_to_jump,
            args.min_fragment_length,
        );
        segments.process_links();

        let (mut iteration_results, cols) = segments.get_results_and_active_cols(
            &alignment_data.query_name_map,
            alignment_data.target_name_map.get(group.target_id),
            group.target_start,
            region_idx,
        );

        results.append(&mut iteration_results);
        active_cols = cols;

        if args.viz {
            trace_strings.push(segments.trace_soda_string(&alignment_data.query_name_map));
            fragment_strings.push(segments.fragments_soda_string(&alignment_data.query_name_map));
            segment_strings.push(segments.segments_soda_string());
        }

        if active_cols.len() == last_num_cols {
            panic!();
        }

        last_num_cols = active_cols.len();
    }

    results.sort_by_key(|r| (r.join_id, r.target_start));
    results.sort_by_key(|r| r.target_start);
    results.retain(|r| r.query_name != "skip");

    Annotation::write(&results, &mut std::io::stdout());

    if args.viz {
        let data = AdjudicationSodaData::new(
            group,
            &alignment_data.query_name_map,
            &results,
            trace_strings,
            fragment_strings,
            segment_strings,
            args,
        );

        let out_path = args.viz_output_path.join(format!(
            "{}-{}-{}.html",
            alignment_data.target_name_map.get(group.target_id),
            matrix_def.target_start,
            matrix_def.target_start + matrix_def.num_cols
        ));

        write_soda_html(
            &data,
            "./fixtures/soda/template.html",
            "./fixtures/soda/viz.js",
            out_path,
        );
    }
}

pub trait StrSliceExt {
    fn to_digital_nucleotides(self) -> Vec<u8>;
}

impl StrSliceExt for &str {
    fn to_digital_nucleotides(self) -> Vec<u8> {
        self.as_bytes()
            .iter()
            .map(|byte| {
                *UTF8_TO_DIGITAL_NUCLEOTIDE
                    .get(byte)
                    .unwrap_or_else(|| panic!("unknown byte: {byte}"))
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::alignment::{Alignment, Strand};

    use super::{run_pipeline, StrSliceExt};

    #[test]
    fn test_pipeline() {
        let ali = [
            Alignment {
                query_id: 1,
                target_seq: "AAAAAAAAAAGGGGGGGGGG".to_digital_nucleotides(),
                query_seq: "AAAAAAAAAAGGGGGGGGGG".to_digital_nucleotides(),
                target_start: 1,
                target_end: 20,
                query_start: 1,
                query_end: 20,
                strand: Strand::Forward,
                substitution_matrix_id: 0,
            },
            Alignment {
                query_id: 2,
                target_seq: "AAAAAAAAAAGGGGGGGGGG".to_digital_nucleotides(),
                query_seq: "AAAAAAAAAAGGGGGGGGGG".to_digital_nucleotides(),
                target_start: 21,
                target_end: 40,
                query_start: 1,
                query_end: 20,
                strand: Strand::Forward,
                substitution_matrix_id: 0,
            },
        ];
    }
}
