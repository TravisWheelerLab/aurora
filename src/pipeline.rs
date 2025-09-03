use std::fs;

use itertools::Itertools;

use crate::{
    alignment::AlignmentData,
    annotation::Annotation,
    chunks::ProximityGroup,
    collapse::AssemblyGraph,
    confidence::confidence,
    matrix::{Matrix, MatrixDef},
    score_params::{approximate_ideal_skip_state_score, ScoreParams},
    support::windowed_confidence,
    viterbi::{trace_segments, traceback, viterbi_collapsed, TraceSegment},
    viz::AdjudicationSodaData,
    windowed_scores::{build_target_seq_from_alignments, windowed_score, Background},
    AuroraArgs,
};

pub fn run_pipeline(
    proximity_group: &ProximityGroup,
    alignment_data: &AlignmentData,
    region_idx: usize,
    mut args: AuroraArgs,
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
    let (confidence_avg_by_id, _confidence_by_id) = windowed_confidence(&mut confidence_matrix);
    let _assembly_graph = AssemblyGraph::new(
        &proximity_group,
        &confidence_avg_by_id,
        &args,
        alignment_data,
    );

    // let mut collapsed_confidence_matrix = Matrix::<f64>::new(&collapsed_matrix_def);
    let mut viterbi_matrix = Matrix::<f64>::new(&matrix_def);
    let mut sources_matrix = Matrix::<usize>::new(&matrix_def);

    // the initial active cols just removes the dead space between alignments
    let active_cols = confidence_matrix.initial_active_cols();
    let mut trace_conclusive: Vec<Vec<TraceSegment>> = vec![];

    // if we're going to produce visualizations, this will
    // keep track of all of the data needed to do so
    let mut soda_data = AdjudicationSodaData::new(
        &proximity_group,
        &confidence_matrix,
        alignment_data,
        &target_seq,
        &args,
    );

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

    let trace_segments = trace_segments(&trace);
    trace_conclusive.push(trace_segments);

    // TODO: function for this
    let mut annotations: Vec<Annotation> = trace_conclusive
        .iter()
        .flat_map(|iter_segments| {
            iter_segments
                .iter()
                .filter(|s| s.ali_id != 0)
                .map(|s| Annotation {
                    target_name: alignment_data
                        .target_name_map
                        .get(proximity_group.target_id)
                        .clone(),
                    target_start: s.col_start + proximity_group.target_start,
                    target_end: s.col_end + proximity_group.target_start,
                    query_id: s.query_id,
                    query_name: match s.row_idx {
                        // 0 is the skip state row
                        // then 1..=(num_assemblies) are alignment rows
                        // so anything >(num_assemblies) is a tandem repeat
                        r if r > proximity_group.alignments.len() => {
                            //
                            let tandem_repeat_idx =
                                s.row_idx - proximity_group.alignments.len() - 1;
                            let repeat = &proximity_group.tandem_repeats[tandem_repeat_idx];
                            format!(
                                "({}:{})#tandem repeat",
                                repeat.period, repeat.consensus_pattern,
                            )
                        }
                        _ => alignment_data.query_name_map.get(s.query_id).clone(),
                    },
                    query_start: viterbi_matrix.consensus_position(s.row_idx, s.col_start),
                    query_end: viterbi_matrix.consensus_position(s.row_idx, s.col_end),
                    strand: viterbi_matrix.strand_of_row(s.row_idx),
                    confidence: (s.col_start..=s.col_end)
                        .map(|col_idx| confidence_matrix.get(s.row_idx, col_idx))
                        .sum::<f64>()
                        / (s.col_end - s.col_start + 1) as f64,
                    join_id: s.row_idx,
                    region_id: region_idx,
                })
        })
        .collect_vec();

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

    annotations.sort_by_key(|r| r.target_start);
    annotations.retain(|r| r.query_name != "skip");

    Annotation::write(&annotations, &mut std::io::stdout());
}
