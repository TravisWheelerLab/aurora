use std::{collections::HashMap, fs};

use itertools::Itertools;

use crate::{
    alignment::AlignmentData,
    alphabet::UTF8_TO_DIGITAL_NUCLEOTIDE,
    chunks::ProximityGroup,
    collapse::AssemblyGroup,
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

    let mut active_cols = collapsed_confidence_matrix.initial_active_cols();
    let mut join_id_cnt = 0usize;
    let mut relevant_rows: Vec<Vec<usize>> = vec![];
    let mut trace_conclusive: Vec<Vec<TraceSegment>> = vec![];
    let mut trace_ambiguous: Vec<Vec<TraceSegment>> = vec![];

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

        let trace_segments = trace_segments(&trace);

        join_id_cnt = trace
            .iter()
            .map(|t| t.join_id)
            .max()
            .expect("trace is empty");

        // -----------
        // row culling
        // -----------

        let mut relevant = vec![];

        assembly_group
            .assemblies
            .iter()
            .enumerate()
            // assemblies with more than one alignment stay relevant
            // **TODO: need to filter for assemblies that have been resolved
            // **TODO: maybe filter based on trace continuity
            .filter(|(_, a)| a.alignments.len() > 1)
            .for_each(|(idx, _)| relevant.push(idx + 1));

        trace_segments.iter().for_each(|s| relevant.push(s.row_idx));

        relevant = relevant.into_iter().unique().collect();

        relevant_rows.push(relevant);

        // ---------------
        // trace splitting
        // ---------------

        let assembly_fragment_ranges = assembly_group
            .assemblies
            .iter()
            .filter(|g| g.alignments.len() > 1)
            .flat_map(|a| &a.alignments)
            .map(|a| {
                (
                    a.target_start - group.target_start,
                    a.target_end - group.target_start,
                )
            })
            .collect_vec();

        let (ambiguous, conclusive): (Vec<TraceSegment>, Vec<TraceSegment>) =
            trace_segments.into_iter().partition(|s| {
                let seg_start = s.col_start + 10;
                let seg_end = s.col_end - 10;

                for &(s, e) in assembly_fragment_ranges.iter() {
                    if seg_start < e && seg_end > s {
                        return true;
                    }
                }
                false
            });

        // ------------------
        // update active cols
        // ------------------

        active_cols = ambiguous
            .iter()
            .flat_map(|s| s.col_start..=s.col_end)
            .collect();

        trace_conclusive.push(conclusive);
        trace_ambiguous.push(ambiguous);
    }

    let annotations: Vec<Annotation> = trace_conclusive
        .iter()
        .flat_map(|iter_segments| iter_segments.iter().map(Annotation::new))
        .collect();

    if args.viz {
        // let trace_strings = trace_clear
        //     .iter()
        //     .map(|seg| {
        //         format!(
        //             "{},{},{},{}",
        //             seg.col_start, seg.col_end, seg.query_id, seg.row_idx
        //         )
        //     })
        //     .join("|");

        let data = AdjudicationSodaData2::new(
            &assembly_group,
            &alignment_data.query_name_map,
            &annotations,
            &vec![],
            vec![],
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
