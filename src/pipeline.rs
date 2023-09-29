use std::fs;

use crate::{
    alignment::Alignment,
    alphabet::UTF8_TO_DIGITAL_NUCLEOTIDE,
    confidence::confidence,
    matrix::{Matrix, MatrixDef},
    results::Annotation,
    score_params::ScoreParams,
    segments::Segments,
    substitution_matrix::SubstitutionMatrix,
    support::windowed_confidence_slow,
    viterbi::{traceback, viterbi},
    viz::AuroraSodaData,
    windowed_scores::windowed_score,
    Args, BACKGROUND_WINDOW_SIZE, SCORE_WINDOW_SIZE,
};

pub fn run_pipeline(
    alignments: &[Alignment],
    query_names: &[String],
    substitution_matrices: &[SubstitutionMatrix],
    region_idx: usize,
    args: &Args,
) {
    let score_params = ScoreParams::new(
        alignments.len(),
        args.query_jump_probability,
        args.num_skip_loops_eq_to_jump,
    );

    let matrix_def = MatrixDef::new(alignments, &alignments[0].target_name, query_names);
    let mut confidence_matrix = Matrix::<f64>::new(&matrix_def);
    let mut viterbi_matrix = Matrix::<f64>::new(&matrix_def);
    let mut sources_matrix = Matrix::<usize>::new(&matrix_def);

    windowed_score(
        &mut confidence_matrix,
        alignments,
        substitution_matrices,
        SCORE_WINDOW_SIZE,
        BACKGROUND_WINDOW_SIZE,
    );

    confidence(&mut confidence_matrix);

    windowed_confidence_slow(&mut confidence_matrix);

    // this removes any initial dead space between alignments
    let mut active_cols = confidence_matrix.initial_active_cols();
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
            args.consensus_join_distance,
        );

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

        let (mut iteration_results, cols) =
            segments.get_results_and_active_cols(&matrix_def, region_idx);
        results.append(&mut iteration_results);
        active_cols = cols;

        if args.viz {
            trace_strings.push(segments.trace_soda_string(&matrix_def));
            fragment_strings.push(segments.fragments_soda_string(&matrix_def));
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
        let data = AuroraSodaData::new(
            &confidence_matrix,
            alignments,
            &results,
            trace_strings,
            fragment_strings,
            segment_strings,
            args,
        );

        let html_template = fs::read_to_string("./fixtures/soda/template.html")
            .expect("failed to read template.html");

        let viz_js = fs::read_to_string("./fixtures/soda/viz.js").expect("failed to read viz.js");

        let mut viz_html = html_template.replace(
            "DATA_TARGET",
            &serde_json::to_string(&data).expect("failed to serialize JSON data"),
        );

        viz_html = viz_html.replace("JS_TARGET", &viz_js);

        let file_path = args.viz_output_path.join(format!(
            "{}-{}-{}.html",
            matrix_def.target_name,
            matrix_def.target_start,
            matrix_def.target_start + matrix_def.num_cols
        ));

        let mut file = std::fs::File::create(file_path).expect("failed to create file");

        std::io::Write::write_all(&mut file, viz_html.as_bytes()).expect("failed to write to file");
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
                target_name: "target".to_string(),
                query_name: "query_1".to_string(),
                query_id: 1,
                target_seq: "AAAAAAAAAAGGGGGGGGGG".to_digital_nucleotides(),
                query_seq: "AAAAAAAAAAGGGGGGGGGG".to_digital_nucleotides(),
                target_start: 1,
                target_end: 20,
                query_start: 1,
                query_end: 20,
                strand: Strand::Forward,
                score: 0,
                substitution_matrix: "".to_string(),
                gap_init: -30.0,
                gap_extend: -5.0,
            },
            Alignment {
                target_name: "target".to_string(),
                query_name: "query_2".to_string(),
                query_id: 2,
                target_seq: "AAAAAAAAAAGGGGGGGGGG".to_digital_nucleotides(),
                query_seq: "AAAAAAAAAAGGGGGGGGGG".to_digital_nucleotides(),
                target_start: 21,
                target_end: 40,
                query_start: 1,
                query_end: 20,
                strand: Strand::Forward,
                score: 0,
                substitution_matrix: "".to_string(),
                gap_init: -30.0,
                gap_extend: -5.0,
            },
        ];
    }
}
