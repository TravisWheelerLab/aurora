use std::{
    fs::{self, File},
    io::{self, Write},
    path::Path,
};

use itertools::{izip, Itertools};

use crate::{
    alignment::AlignmentData,
    annotation::{AmbiguousAnnotation, SimpleAnnotation},
    assembly::AssemblyGraph,
    chunks::ProximityGroup,
    confidence::confidence,
    matrix::{Matrix, MatrixDef},
    score_params::{approximate_ideal_skip_state_score, ScoreParams},
    segments::{segments_from_matrix_trace, SegmentedMatrix},
    support::windowed_confidence,
    viterbi::{
        backtrace_histories, history_score, history_viterbi_on_segments, trace_segments, traceback,
        viterbi_collapsed, AnnotatedRange, History, HistoryEntry, RefinedTraceSegment,
    },
    viz::AdjudicationSodaData,
    windowed_scores::{build_target_seq_from_alignments, windowed_score, Background},
    AuroraArgs,
};

pub fn to_annotations(
    proximity_group: &ProximityGroup,
    alignment_data: &AlignmentData,
    confidence_matrix: &Matrix<f64>,
    trace_segments: &[RefinedTraceSegment],
    region_idx: usize,
) -> Vec<AmbiguousAnnotation> {
    trace_segments
        .iter()
        .filter(|v| v.annotated.iter().any(|a| a.row_idx != 0))
        .map(|s| {
            let confidence = s
                .annotated
                .iter()
                .map(|a| {
                    (a.col_start..=a.col_end)
                        .map(|col_idx| confidence_matrix.get(a.row_idx, col_idx))
                        .sum::<f64>()
                        / (a.col_end - a.col_start + 1) as f64
                })
                .sum::<f64>()
                / (s.annotated.len().max(1) as f64);

            AmbiguousAnnotation {
                target_name: alignment_data
                    .target_name_map
                    .get(proximity_group.target_id)
                    .clone(),
                annotations: s
                    .annotated
                    .iter()
                    .map(|a| {
                        SimpleAnnotation {
                            target_start: a.col_start + proximity_group.target_start,
                            target_end: a.col_end + proximity_group.target_start,
                            query_id: a.query_id.unwrap_or(0),
                            query_name: match a.row_idx {
                                // 0 is the skip state row
                                // then 1..=(num_assemblies) are alignment rows
                                // so anything >(num_assemblies) is a tandem repeat
                                r if r > proximity_group.alignments.len() => {
                                    //
                                    let tandem_repeat_idx =
                                        a.row_idx - proximity_group.alignments.len() - 1;
                                    let repeat = &proximity_group.tandem_repeats[tandem_repeat_idx];
                                    format!(
                                        "({}:{})#tandem repeat",
                                        repeat.period, repeat.consensus_pattern,
                                    )
                                }
                                _ => alignment_data
                                    .query_name_map
                                    .get(a.query_id.expect("Annotation has no query id!"))
                                    .clone(),
                            },
                            query_start: confidence_matrix
                                .consensus_position(a.row_idx, a.col_start),
                            query_end: confidence_matrix.consensus_position(a.row_idx, a.col_end),
                            strand: confidence_matrix.strand_of_row(a.row_idx),
                        }
                    })
                    .collect_vec(),
                confidence,
                join_id: s.join_index,
                region_id: region_idx,
            }
        })
        .collect_vec()
}

fn get_history_lengths(history: &History) -> Vec<usize> {
    history
        .segment_offsets
        .iter()
        .skip(1)
        .zip(
            history
                .segment_offsets
                .iter()
                .skip(2)
                .chain([history.entries.len()].iter()),
        )
        .map(|(a, b)| b - a)
        .collect_vec()
}

fn dump_history_scores(history: &History, path: impl AsRef<Path>) -> io::Result<()> {
    let mut file = File::create(path)?;

    history
        .segment_offsets
        .iter()
        .skip(1)
        .zip(
            history
                .segment_offsets
                .iter()
                .skip(2)
                .chain([history.entries.len()].iter()),
        )
        .try_for_each(|(&start, &end)| {
            for entry in history.entries[start..end].iter() {
                if let HistoryEntry::Append(val) | HistoryEntry::Join(val) = entry {
                    write!(&mut file, "{:e} ", val.score)?;
                }
            }
            writeln!(&mut file)
        })
}

fn dump_final_trace_statistics(
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

fn dump_debug_history_info(
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

    writeln!(
        &mut file,
        "Segment, History Count, Group Count, Index Count, Block Count, Target Start, Target End"
    )?;

    izip!(
        (0..segments.len()),
        history_lengths.iter(),
        num_groups,
        group_sizes,
        segment_lengths,
        segment_ranges
    )
    .try_for_each(|v| {
        writeln!(
            &mut file,
            "{}, {}, {}, {}, {}, {}, {}",
            v.0, v.1, v.2, v.3, v.4, v.5 .0, v.5 .1
        )
    })?;

    Ok(())
}

pub fn run_pipeline(
    proximity_group: &ProximityGroup,
    alignment_data: &AlignmentData,
    region_idx: usize,
    mut args: AuroraArgs,
    output_file: &mut impl Write,
    ambiguity_file: Option<&mut impl Write>,
) {
    let annot_args = &args.annotation_args;

    if args.visualization_args.viz {
        args.visualization_args
            .viz_output_path
            .push(format!("{}", region_idx));
        fs::create_dir_all(&args.visualization_args.viz_output_path).unwrap();
    }

    let vis_args = &args.visualization_args;

    let score_params = ScoreParams::new(
        proximity_group.alignments.len(),
        annot_args.query_jump_penalty,
        annot_args.num_skip_loops_eq_to_jump,
    );

    let matrix_def = MatrixDef::from_proximity_group(proximity_group);

    let mut confidence_matrix = Matrix::<f64>::new(&matrix_def);

    let target_start = proximity_group.target_start;
    let target_end = proximity_group.target_end;
    let target_length = target_end - target_start + 1;

    let target_seq =
        build_target_seq_from_alignments(proximity_group.alignments, target_start, target_end);

    let background = Background::new(
        &target_seq,
        target_start,
        target_length,
        args.annotation_args.background_window_size,
    );

    let skip_state_score = approximate_ideal_skip_state_score(
        annot_args.num_skip_loops_eq_to_jump as f64,
        annot_args.query_jump_penalty,
        annot_args.skip_state_score_shift,
    );

    windowed_score(
        &mut confidence_matrix,
        proximity_group.alignments,
        proximity_group.tandem_repeats,
        &alignment_data.substitution_matrices,
        &background,
        args.annotation_args.score_window_size,
        skip_state_score,
    )
    .unwrap();

    confidence(&mut confidence_matrix);
    let (_confidence_avg_by_id, _confidence_by_id) = windowed_confidence(&mut confidence_matrix);

    let assembly_graph = AssemblyGraph::new(proximity_group, &score_params, &args.annotation_args);
    let segments;
    let simple_trace;

    // In a new block so initial viterbi matricies/sources are freed right after being used...
    {
        let mut sources_matrix = Matrix::<usize>::new(&matrix_def);
        let mut viterbi_matrix = Matrix::<f64>::new(&matrix_def);
        // the initial active cols just removes the dead space between alignments
        let active_cols = confidence_matrix.initial_active_cols();

        viterbi_collapsed(
            &confidence_matrix,
            &mut viterbi_matrix,
            &mut sources_matrix,
            &active_cols,
            &score_params,
        );

        let trace = traceback(
            &viterbi_matrix,
            &confidence_matrix,
            &sources_matrix,
            &active_cols,
        );

        simple_trace = trace_segments(&trace);

        segments = segments_from_matrix_trace(
            proximity_group,
            &simple_trace,
            &confidence_matrix,
            &score_params,
            &assembly_graph,
            &args.annotation_args,
        );
    }

    let history_lengths;
    let refined_trace_segments;

    if !vis_args.disable_tracing {
        let history = history_viterbi_on_segments(
            &segments,
            &score_params,
            &assembly_graph,
            args.annotation_args.max_history_depth,
            args.annotation_args.max_histories_per_segment,
            args.annotation_args.min_relative_history_score,
        );

        history_lengths = get_history_lengths(&history);

        if args.visualization_args.debug {
            dump_debug_history_info(
                &history,
                &segments,
                proximity_group.target_start,
                &history_lengths,
                vis_args.viz_output_path.join("history_info.csv"),
            )
            .map_err(|_| eprintln!("Unable to save debug history info!"))
            .ok();
            dump_history_scores(
                &history,
                vis_args.viz_output_path.join("history_scores.txt"),
            )
            .map_err(|_| eprintln!("Unable to save history scores!"))
            .ok();
        }

        refined_trace_segments = backtrace_histories(&segments, &history, region_idx);

        if args.visualization_args.debug {
            dump_final_trace_statistics(
                &history,
                &segments,
                &refined_trace_segments,
                vis_args.viz_output_path.join("final_trace_stats.csv"),
            )
            .map_err(|_| eprintln!("Unable to save final trace statistics!"))
            .ok();
        }
    } else {
        refined_trace_segments = simple_trace
            .iter()
            .enumerate()
            .map(|(i, v)| RefinedTraceSegment {
                annotated: vec![AnnotatedRange {
                    query_id: Some(v.query_id),
                    row_idx: v.row_idx,
                    col_start: v.col_start,
                    col_end: v.col_end,
                }],
                join_index: i,
                score: 0.0,
                segment: i,
            })
            .collect_vec();
        history_lengths = vec![0; segments.len()];
    }

    // if we're going to produce visualizations, this will
    // keep track of all of the data needed to do so
    let mut soda_data = AdjudicationSodaData::new(
        proximity_group,
        &confidence_matrix,
        alignment_data,
        &target_seq,
        &refined_trace_segments,
        &segments,
        &history_lengths,
        &assembly_graph,
        args.visualization_args.viz_enable_scores,
        &args,
    );

    // Grab the annotations...
    let mut annotations: Vec<AmbiguousAnnotation> = to_annotations(
        proximity_group,
        alignment_data,
        &confidence_matrix,
        &refined_trace_segments,
        region_idx,
    );

    if vis_args.viz {
        // TODO: this is kind of awkward
        soda_data.set_annotations(annotations.clone());

        let out_path = vis_args.viz_output_path.join("index.html");
        soda_data.write(out_path);
    }

    if !vis_args.viz_constraints.is_empty() {
        let target_name = alignment_data
            .target_name_map
            .get(proximity_group.target_id)
            .clone();

        let target_start = proximity_group.target_start;
        let target_end = proximity_group.target_end;

        let constraints = vis_args
            .viz_constraints
            .iter()
            .filter(|c| c.target_name == target_name)
            .filter(|c| c.target_start < target_end && c.target_end > target_start)
            .collect_vec();

        constraints.iter().for_each(|constraint| {
            let out_path = vis_args.viz_output_path.parent().unwrap().join(format!(
                "{}-{}-{}.html",
                constraint.target_name, constraint.target_start, constraint.target_end
            ));
            soda_data.constrain(constraint);
            soda_data.write(out_path);
        });
    }

    annotations.sort_by_key(|r| r.annotations.iter().map(|a| a.target_start).min());
    annotations.retain(|r| r.annotations.iter().any(|a| a.query_name != "skip"));

    AmbiguousAnnotation::write(&annotations, output_file, true);
    if let Some(ambig_file) = ambiguity_file {
        AmbiguousAnnotation::write(&annotations, ambig_file, false);
    }
}
