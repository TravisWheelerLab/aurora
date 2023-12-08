use std::fs;

use itertools::Itertools;

use crate::{
    alignment::AlignmentData,
    alphabet::UTF8_TO_DIGITAL_NUCLEOTIDE,
    annotation::Annotation,
    chunks::ProximityGroup,
    collapse::AssemblyGroup,
    confidence::confidence,
    matrix::{Matrix, MatrixDef},
    score_params::ScoreParams,
    split::{split_trace, MatrixRange},
    support::windowed_confidence,
    viterbi::{trace_segments, traceback, viterbi_collapsed, TraceSegment},
    viz::{write_soda_html, AdjudicationSodaData},
    windowed_scores::windowed_score,
    Args, BACKGROUND_WINDOW_SIZE, SCORE_WINDOW_SIZE,
};

pub fn run_pipeline(
    proximity_group: &ProximityGroup,
    alignment_data: &AlignmentData,
    region_idx: usize,
    mut args: Args,
) {
    if args.viz {
        args.viz_output_path.push(format!(
            "{}-{}-{}-{}",
            proximity_group.target_id,
            region_idx,
            proximity_group.target_start,
            proximity_group.target_end
        ));

        fs::create_dir_all(&args.viz_output_path).unwrap();
    }

    let score_params = ScoreParams::new(
        proximity_group.alignments.len(),
        args.query_jump_probability,
        args.num_skip_loops_eq_to_jump,
    );

    let matrix_def = MatrixDef::from_proximity_group(proximity_group);
    let mut confidence_matrix = Matrix::<f64>::new(&matrix_def);

    windowed_score(
        &mut confidence_matrix,
        proximity_group.alignments,
        proximity_group.tandem_repeats,
        &alignment_data.substitution_matrices,
        SCORE_WINDOW_SIZE,
        BACKGROUND_WINDOW_SIZE,
    );

    confidence(&mut confidence_matrix);

    let (confidence_avg_by_id, confidence_by_id) = windowed_confidence(&mut confidence_matrix);

    // adjust the skip state to include skip-loop penalty
    let skip_adjust = args
        .query_jump_probability
        .powf(1.0 / args.num_skip_loops_eq_to_jump as f64);

    (0..confidence_matrix.num_cols()).for_each(|col_idx| {
        confidence_matrix.set_skip(col_idx, confidence_matrix.get_skip(col_idx) * skip_adjust);
    });

    // convert the ProximityGroup into an AssemblyGroup
    let assembly_group = AssemblyGroup::new(
        proximity_group,
        &confidence_avg_by_id,
        &confidence_by_id,
        &args,
    );

    // initialize the collpased DP matrices
    let collapsed_matrix_def = MatrixDef::from_assembly_group(&assembly_group);

    let mut collapsed_confidence_matrix = Matrix::<f64>::new(&collapsed_matrix_def);
    let mut viterbi_matrix = Matrix::<f64>::new(&collapsed_matrix_def);
    let mut sources_matrix = Matrix::<usize>::new(&collapsed_matrix_def);

    collapsed_confidence_matrix.copy_fill(&confidence_matrix);

    // the initial active cols just removes the dead space between alignments
    let mut active_cols = collapsed_confidence_matrix.initial_active_cols();
    let mut inactive_segments: Vec<Vec<MatrixRange>> = vec![];
    let mut trace_conclusive: Vec<Vec<TraceSegment>> = vec![];
    let mut trace_ambiguous: Vec<Vec<TraceSegment>> = vec![];
    let mut resolved_assembly_rows: Vec<Vec<usize>> = vec![];
    let mut unresolved_assembly_rows: Vec<Vec<usize>> = vec![];
    let mut competed_assembly_rows: Vec<Vec<usize>> = vec![];

    while !active_cols.is_empty() {
        viterbi_collapsed(
            &collapsed_confidence_matrix,
            &mut viterbi_matrix,
            &mut sources_matrix,
            &active_cols,
            &score_params,
        );

        let trace = traceback(
            &viterbi_matrix,
            &collapsed_confidence_matrix,
            &sources_matrix,
            &active_cols,
        );

        // we should always have one trace step for every active column
        debug_assert_eq!(trace.len(), active_cols.len());

        let trace_segments = trace_segments(&trace);

        let (ambiguous, conclusive, resolved, unresolved, competed, inactive) = split_trace(
            trace_segments,
            &assembly_group,
            &active_cols,
            &confidence_avg_by_id,
            &args,
        );

        let new_active_cols = ambiguous
            .iter()
            .flat_map(|s| s.col_start..=s.col_end)
            .collect_vec();

        debug_assert!(new_active_cols.len() <= active_cols.len());

        inactive_segments.push(inactive);
        trace_conclusive.push(conclusive);
        trace_ambiguous.push(ambiguous);
        resolved_assembly_rows.push(resolved);
        unresolved_assembly_rows.push(unresolved);
        competed_assembly_rows.push(competed);

        if new_active_cols.len() == active_cols.len() {
            break;
        }

        active_cols = new_active_cols;
    }

    let final_ambigous = trace_ambiguous.last().expect("no ambiguous trace").to_vec();
    trace_conclusive.push(final_ambigous);

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
                        r if r > assembly_group.assemblies.len() => {
                            //
                            let tandem_repeat_idx = s.row_idx - assembly_group.assemblies.len() - 1;
                            let repeat = &assembly_group.tandem_repeats[tandem_repeat_idx];
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
                        .map(|col_idx| collapsed_confidence_matrix.get(s.row_idx, col_idx))
                        .sum::<f64>()
                        / (s.col_end - s.col_start + 1) as f64,
                    join_id: s.row_idx,
                    region_id: region_idx,
                })
        })
        .collect_vec();

    annotations.sort_by_key(|r| r.target_start);
    annotations.retain(|r| r.query_name != "skip");

    Annotation::write(&annotations, &mut std::io::stdout());

    if args.viz {
        let data = AdjudicationSodaData::new(
            &assembly_group,
            &collapsed_confidence_matrix,
            alignment_data,
            &annotations,
            trace_conclusive,
            trace_ambiguous,
            resolved_assembly_rows,
            unresolved_assembly_rows,
            competed_assembly_rows,
            inactive_segments,
            &args,
        );

        let out_path = args.viz_output_path.join(format!(
            "{}-{}-{}.html",
            alignment_data
                .target_name_map
                .get(proximity_group.target_id),
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
