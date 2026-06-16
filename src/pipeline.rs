use std::fmt::Display;

use itertools::Itertools;

use crate::{
    alignment::{AlignmentData, Strand},
    annotation::{AmbiguousAnnotation, SimpleAnnotation},
    assembly::gather_join_statistics,
    chunks::ProximityGroup,
    confidence::confidence,
    history_tracing::{
        backtrace_histories, history_viterbi_on_segments, History, RefinedTraceSegment,
    },
    join_estimation::{JoinEstimator, JoinStatisticsCollector},
    matrix::{Matrix, MatrixDef},
    score_params::{approximate_ideal_skip_state_score, ScoreParams},
    segments::{assemble_and_link_segments, segments_from_matrix_trace, InitialSegments},
    support::windowed_confidence,
    trace_statistics::TraceStatistics,
    viterbi::{trace_segments, traceback, viterbi_collapsed, TraceSegment},
    viz::{
        debug::{dump_debug_history_info, dump_final_trace_statistics},
        AdjudicationSodaDataArgs, AdjudicationSodaWriter,
    },
    windowed_scores::{build_target_seq_from_alignments, windowed_score, Background},
    AuroraArgs,
};

pub fn to_annotations(
    proximity_group: &ProximityGroup,
    alignment_data: &AlignmentData,
    trace_segments: &[RefinedTraceSegment],
    region_idx: usize,
) -> Vec<AmbiguousAnnotation> {
    trace_segments
        .iter()
        .filter(|v| v.annotated.iter().any(|a| a.row_idx != 0))
        .map(|s| {
            let confidence = s.annotated.iter().map(|a| a.avg_confidence).sum::<f64>()
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
                                        "({}:{})#tandem-repeat",
                                        repeat.period, repeat.consensus_pattern,
                                    )
                                }
                                _ => alignment_data
                                    .query_name_map
                                    .get(a.query_id.expect("Annotation has no query id!"))
                                    .clone(),
                            },
                            query_start: a.query_start,
                            query_end: a.query_end,
                            strand: if a.row_idx > 0
                                && a.row_idx <= proximity_group.alignments.len()
                            {
                                proximity_group.alignments[a.row_idx - 1].strand
                            } else {
                                Strand::Forward
                            },
                            kimura80: if a.row_idx > 0
                                && a.row_idx <= proximity_group.alignments.len()
                            {
                                proximity_group.alignments[a.row_idx - 1]
                                    .kimura80(a.query_start, a.query_end)
                            } else {
                                0.0
                            },
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

fn get_active_columns<T: Copy + Default + Display>(matrix: &Matrix<T>) -> Vec<(usize, usize)> {
    let mut active_cols = Vec::new();
    let mut prior_start = None;
    let mut prior_end: usize = 0;
    let cols = matrix.initial_active_cols();

    for &col in cols.iter() {
        if let Some(val) = prior_start {
            if prior_end + 1 != col {
                active_cols.push((val, prior_end));
                prior_start = Some(col);
            }

            prior_end = col;
        } else {
            prior_start = Some(col);
            prior_end = col;
        }
    }

    if let Some(val) = prior_start {
        active_cols.push((val, prior_end));
    }

    active_cols
}

pub struct NaiveTraceResults<T: JoinStatisticsCollector> {
    pub target_start: usize,
    pub target_end: usize,
    pub trace_segments: Vec<TraceSegment>,
    pub segments: InitialSegments,
    pub score_params: ScoreParams,
    pub alignment_confidences: Vec<f64>,
    pub active_columns: Vec<(usize, usize)>,
    pub query_join_statistics: Vec<(usize, T)>,
    pub viz_writer: AdjudicationSodaWriter,
    pub region_index: usize,
}

pub fn run_naive_trace<T: JoinStatisticsCollector>(
    proximity_group: &ProximityGroup,
    alignment_data: &AlignmentData,
    region_idx: usize,
    args: &AuroraArgs,
) -> NaiveTraceResults<T> {
    let annot_args = &args.annotation_args;

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
        build_target_seq_from_alignments(proximity_group.alignments, target_start, target_length);

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
    let confidence_by_row = windowed_confidence(
        &mut confidence_matrix,
        args.annotation_args.score_window_size,
    );

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

        let trace = traceback(&viterbi_matrix, &sources_matrix, &active_cols);

        simple_trace = trace_segments(&trace);

        segments = segments_from_matrix_trace(
            proximity_group,
            &simple_trace,
            &confidence_matrix,
            &score_params,
            &args.annotation_args,
        );
    }

    let mut viz_writer = AdjudicationSodaWriter::new(
        proximity_group,
        alignment_data,
        &args.visualization_args.viz_output_path,
        region_idx,
        &args.visualization_args.viz_constraints,
    );

    if args.visualization_args.viz && args.visualization_args.viz_enable_scores {
        viz_writer
            .write_confidences(&confidence_matrix)
            .expect("Unable to write confidences!!!");
    }

    let query_join_statistics = gather_join_statistics(
        proximity_group,
        &segments,
        &alignment_data.query_lengths,
        &args.annotation_args,
    );

    NaiveTraceResults {
        target_start: proximity_group.target_start,
        target_end: proximity_group.target_end,
        trace_segments: simple_trace,
        segments,
        score_params,
        alignment_confidences: confidence_by_row,
        active_columns: get_active_columns(&confidence_matrix),
        query_join_statistics,
        viz_writer,
        region_index: region_idx,
    }
}

pub fn run_history_trace<T: JoinEstimator, S: JoinStatisticsCollector>(
    proximity_group: &ProximityGroup,
    alignment_data: &AlignmentData,
    trace_statistics: &TraceStatistics<T>,
    naive_trace: &mut NaiveTraceResults<S>,
    args: &AuroraArgs,
) -> Vec<AmbiguousAnnotation> {
    let vis_args = &args.visualization_args;

    let (segments, assembly_graph) = assemble_and_link_segments(
        proximity_group,
        &mut naive_trace.segments,
        &naive_trace.trace_segments,
        &trace_statistics.region_statistics[naive_trace.region_index],
        &trace_statistics.query_statistics,
        &naive_trace.score_params,
        &args.annotation_args,
        &alignment_data.query_lengths,
    );

    let history = history_viterbi_on_segments(
        segments,
        &naive_trace.score_params,
        &assembly_graph,
        args.annotation_args.max_history_depth,
        args.annotation_args.max_histories_per_segment,
        args.annotation_args.min_relative_history_score,
    );

    let history_lengths = get_history_lengths(&history);

    if args.visualization_args.debug {
        dump_debug_history_info(
            &history,
            segments,
            proximity_group.target_start,
            &history_lengths,
            vis_args.viz_output_path.join("history_info.csv"),
        )
        .map_err(|_| eprintln!("Unable to save debug history info!"))
        .ok();
    }

    let refined_trace_segments = backtrace_histories(segments, &history, naive_trace.region_index);

    if args.visualization_args.debug {
        dump_final_trace_statistics(
            &history,
            segments,
            &refined_trace_segments,
            vis_args.viz_output_path.join("final_trace_stats.csv"),
        )
        .map_err(|_| eprintln!("Unable to save final trace statistics!"))
        .ok();
    }

    // Grab the annotations...
    let annotations: Vec<AmbiguousAnnotation> = to_annotations(
        proximity_group,
        alignment_data,
        &refined_trace_segments,
        naive_trace.region_index,
    );

    if vis_args.viz {
        naive_trace
            .viz_writer
            .write(AdjudicationSodaDataArgs {
                group: proximity_group,
                alignment_confidences: &naive_trace.alignment_confidences,
                active_columns: &naive_trace.active_columns,
                alignment_data,
                annotations: &annotations,
                target_seq: &build_target_seq_from_alignments(
                    proximity_group.alignments,
                    proximity_group.target_start,
                    proximity_group.target_end - proximity_group.target_start + 1,
                ),
                trace: &refined_trace_segments,
                segments,
                history_counts: &get_history_lengths(&history),
                links: &assembly_graph,
                viz_args: vis_args,
            })
            .expect("Unable to write visualization!");
    }

    annotations
}
