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
    windowed_scores::windowed_score,
    BACKGROUND_WINDOW_SIZE, CONSENSUS_JOIN_DISTANCE, SCORE_WINDOW_SIZE, TARGET_JOIN_DISTANCE,
};

pub fn run_pipeline(
    alignments: &[Alignment],
    query_names: &[String],
    substitution_matrices: &[SubstitutionMatrix],
    region_idx: usize,
) {
    let score_params = ScoreParams::new(alignments.len());

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

    // let mut json: Vec<serde_json::Value> = vec![];
    let mut join_id_cnt = 1usize;
    while !active_cols.is_empty() {
        viterbi(
            &confidence_matrix,
            &mut viterbi_matrix,
            &mut sources_matrix,
            &active_cols,
            &score_params,
        );

        let trace = traceback(&viterbi_matrix, &sources_matrix, &active_cols, join_id_cnt);
        join_id_cnt = trace[*active_cols.last().expect("active_cols is empty")].join_id;

        let mut segments = Segments::new(&trace, &confidence_matrix, &active_cols);
        segments.create_links(CONSENSUS_JOIN_DISTANCE, TARGET_JOIN_DISTANCE);
        segments.process_links();
        // json.push(segments.json(&matrix_def));

        let (mut removed_results, cols) =
            segments.get_results_and_active_cols(&matrix_def, region_idx);
        results.append(&mut removed_results);
        active_cols = cols;

        if active_cols.len() == last_num_cols {
            panic!();
        }

        last_num_cols = active_cols.len();
    }

    // let mut file = std::fs::File::create("a.json").expect("failed to create a.json");
    // file.write_all(serde_json::to_string(&json).unwrap().as_bytes())
    //     .expect("failed to write to a.json");

    results.sort_by_key(|r| (r.join_id, r.target_start));
    results.sort_by_key(|r| r.target_start);
    results.retain(|r| r.query_name != "skip");

    Annotation::write(&results, &mut std::io::stdout());
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
