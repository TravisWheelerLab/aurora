use std::{fs, path::PathBuf};

use crate::{
    alignment::AlignmentData,
    alphabet::UTF8_TO_DIGITAL_NUCLEOTIDE,
    chunks::ProximityGroup,
    collapse::collapse,
    confidence::confidence,
    matrix::{Matrix, MatrixDef},
    results::Annotation,
    score_params::ScoreParams,
    segments::Segments,
    support::windowed_confidence_slow,
    viterbi::{traceback, viterbi},
    viz::{write_soda_html, AuroraAdjudicationSodaData},
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

    let skip_confidence_by_col = confidence(&mut confidence_matrix);

    let confidence_avg_by_row = windowed_confidence_slow(&mut confidence_matrix);

    let assemblies = collapse(group, &confidence_avg_by_row, &args);

    let (assembled, single): (Vec<_>, Vec<_>) = assemblies.into_iter().partition(|a| a.len() > 1);

    // thoughts:
    //  - find out how often do assemblies overlap in the target
    //  - should be able to filter reasonably well by comparing groups of rows that are
    //    significantly long (50+?) and reasonably aligned (+/- 10 on either side?, maybe a
    //    percentage of length))
    //
    // vague steps:
    //  - probably don't filter grouped
    //  - probably try to filter the singles
    //  - rebuild matrix
    //  - gotta take care of target overlapping assemblies (mini-viterbi)
    //  - run DP on reduced matrix
    //    - try same-row consensus-based loop transitions
    //  - remove all columns called by singles that are either clearly
    //    contained, or don't compete with an assembly
    //  - run again
    //  - remove resolved assemblies (all pieces had the chance to join)
    //  - iterate
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
        // );

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
        let data = AuroraAdjudicationSodaData::new(
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
