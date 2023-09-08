use crate::{matrix::Matrix, score_params::ScoreParams, CONSENSUS_JOIN_DISTANCE};

use itertools::multizip;

pub fn viterbi(
    confidence_matrix: &Matrix<f64>,
    viterbi_matrix: &mut Matrix<f64>,
    sources_matrix: &mut Matrix<usize>,
    active_cols: &[usize],
    score_params: &ScoreParams,
) {
    let first_col_idx = active_cols[0];

    // initialize the first viterbi column as the
    // log of the first column of confidence values
    for (sparse_row_idx, (confidence_value, viterbi_score, source_value)) in multizip((
        confidence_matrix.col_slice(first_col_idx),
        viterbi_matrix.col_slice_mut(first_col_idx),
        sources_matrix.col_slice_mut(first_col_idx),
    ))
    .enumerate()
    {
        *viterbi_score = confidence_value.ln();
        *source_value = sparse_row_idx;
    }

    // iterate over each of the remaining columns
    for (&col_from_idx, &col_to_idx) in active_cols.iter().zip(active_cols.iter().skip(1)) {
        // TODO: we could keep track of this as they are
        //       being set in the previous loop iteration
        let (sparse_row_idx_of_max_score_in_col_from, max_score_in_prev_col) = viterbi_matrix
            .col_slice(col_from_idx)
            .iter()
            .enumerate()
            // skip the skip state
            .skip(1)
            .fold((0usize, -f64::INFINITY), |acc, (sparse_row_idx, &score)| {
                if score > acc.1 {
                    (sparse_row_idx, score)
                } else {
                    acc
                }
            });

        // to get the skip state score, we only have to
        // compare the cost of looping vs the cost of
        // jumping from the best score in the previous col
        let skip_loop_score = viterbi_matrix.get_skip(col_from_idx) + score_params.skip_loop_score;
        let query_to_skip_score = max_score_in_prev_col + score_params.query_to_skip_score;

        if skip_loop_score > query_to_skip_score {
            viterbi_matrix.set_skip(
                col_to_idx,
                confidence_matrix.get_skip(col_to_idx).ln() + skip_loop_score,
            );
            sources_matrix.set_skip(col_to_idx, 0);
        } else {
            viterbi_matrix.set_skip(
                col_to_idx,
                confidence_matrix.get_skip(col_to_idx).ln() + query_to_skip_score,
            );
            sources_matrix.set_skip(col_to_idx, sparse_row_idx_of_max_score_in_col_from);
        }

        let col_from_length = viterbi_matrix.col_length(col_from_idx);
        let col_to_length = viterbi_matrix.col_length(col_to_idx);

        // precompute the sparse rows by id
        let mut sparse_rows_from_by_id: Vec<Vec<usize>> = vec![vec![]; viterbi_matrix.num_ids()];
        let mut sparse_rows_to_by_id: Vec<Vec<usize>> = vec![vec![]; viterbi_matrix.num_ids()];

        (0..col_from_length).for_each(|sparse_row_idx| {
            let id = viterbi_matrix.id_of_cell_sparse(sparse_row_idx, col_from_idx);
            sparse_rows_from_by_id[id].push(sparse_row_idx);
        });

        (0..col_to_length).for_each(|sparse_row_idx| {
            let id = viterbi_matrix.id_of_cell_sparse(sparse_row_idx, col_to_idx);
            sparse_rows_to_by_id[id].push(sparse_row_idx);
        });

        // now we need to build the column-dependent join transition matrix
        // we are going to over-allocate to make indexing trivial
        let mut join_probs = vec![vec![0.0; col_to_length]; col_from_length];
        (0..col_from_length).for_each(|sparse_row_from_idx| {
            let id = viterbi_matrix.id_of_cell_sparse(sparse_row_from_idx, col_from_idx);
            let consensus_pos_from =
                viterbi_matrix.consensus_position_sparse(sparse_row_from_idx, col_from_idx);

            let sparse_rows_to = &sparse_rows_to_by_id[id];

            debug_assert!({
                let mut others = vec![];
                (0..col_to_length).for_each(|sparse_row_to_idx| {
                    if viterbi_matrix.id_of_cell_sparse(sparse_row_to_idx, col_to_idx) == id {
                        others.push(sparse_row_to_idx);
                    }
                });
                others == *sparse_rows_to
            });

            let mut join_delta_sum = 0.0;
            let strand_from =
                viterbi_matrix.strand_of_cell_sparse(sparse_row_from_idx, col_from_idx);
            sparse_rows_to
                .iter()
                // we need to make sure we don't create a
                // join transition between opposite strands
                .filter(|&&sparse_row_idx| {
                    viterbi_matrix.strand_of_cell_sparse(sparse_row_idx, col_to_idx) == strand_from
                })
                .for_each(|&sparse_row_to_idx| {
                    debug_assert_eq!(
                        viterbi_matrix.id_of_cell_sparse(sparse_row_from_idx, col_from_idx),
                        viterbi_matrix.id_of_cell_sparse(sparse_row_to_idx, col_to_idx),
                    );
                    debug_assert_eq!(
                        viterbi_matrix.strand_of_cell_sparse(sparse_row_from_idx, col_from_idx),
                        viterbi_matrix.strand_of_cell_sparse(sparse_row_to_idx, col_to_idx),
                    );

                    let consensus_pos_to =
                        viterbi_matrix.consensus_position_sparse(sparse_row_to_idx, col_to_idx);

                    let consensus_delta = consensus_pos_from.abs_diff(consensus_pos_to);
                    let join_delta = CONSENSUS_JOIN_DISTANCE.saturating_sub(consensus_delta);
                    join_delta_sum += join_delta as f64;
                    join_probs[sparse_row_from_idx][sparse_row_to_idx] = join_delta as f64;
                });

            // move the probabilities into log space
            if join_delta_sum != 0.0 {
                // if join_delta_sum is greater than zero,
                // that means we have at least one valid transition
                sparse_rows_to.iter().for_each(|&sparse_row_to_idx| {
                    join_probs[sparse_row_from_idx][sparse_row_to_idx] =
                        (join_probs[sparse_row_from_idx][sparse_row_to_idx] / join_delta_sum).ln()
                });
            } else {
                // if join_delta_sum is zero,
                // that means we had no valid transitions
                sparse_rows_to.iter().for_each(|&sparse_row_to_idx| {
                    join_probs[sparse_row_from_idx][sparse_row_to_idx] = -f64::INFINITY;
                });
            }
        });

        let join_scores = join_probs;
        // now we need to pull the scores in
        // for every active cell in the column
        (0..col_to_length)
            // skip the skip state, since
            // we already did it
            .skip(1)
            .for_each(|sparse_row_to_idx| {
                // first figure out what score we'd
                // get if we came from the skip state
                let skip_to_query_score =
                    viterbi_matrix.get_skip(col_from_idx) + score_params.query_to_skip_score;

                // since the sparse row indices aren't the same
                // between columns, we need to find the logical
                // row idx of the current cell we are computing
                let logical_row_to_idx =
                    viterbi_matrix.sparse_to_logical_row_idx(sparse_row_to_idx, col_to_idx);

                let id = viterbi_matrix.id_of_row(logical_row_to_idx);
                let sparse_rows_from = &sparse_rows_from_by_id[id];

                let loop_or_join_tuple = sparse_rows_from.iter().fold(
                    (-f64::INFINITY, 0usize),
                    |acc, &sparse_row_from_idx| {
                        let logical_row_from_idx = viterbi_matrix
                            .sparse_to_logical_row_idx(sparse_row_from_idx, col_from_idx);
                        let num_rows_in_col_from = viterbi_matrix.col_length(col_from_idx);
                        let is_join_transition = logical_row_from_idx != logical_row_to_idx;
                        let score = viterbi_matrix.get_sparse(sparse_row_from_idx, col_from_idx)
                            + join_scores[sparse_row_from_idx][sparse_row_to_idx];

                        if score > acc.0 {
                            (
                                score,
                                sparse_row_from_idx
                                // weird: if this is a join transition (i.e. we came
                                //        from a different logical row), then we add an
                                //        offset equal to the number of rows in the column
                                //        we came from, so that we can decode this
                                //        differently when we are doing a traceback
                                    + (is_join_transition as usize * num_rows_in_col_from),
                            )
                        } else {
                            acc
                        }
                    },
                );

                let query_jump_score = max_score_in_prev_col + score_params.query_jump_score;

                let score_and_source_tuples = [
                    loop_or_join_tuple,
                    (query_jump_score, sparse_row_idx_of_max_score_in_col_from),
                    (skip_to_query_score, 0),
                ];

                let (score, source) = score_and_source_tuples.iter().fold(
                    (-f64::INFINITY, 0usize),
                    |acc, &(score, source)| {
                        if score > acc.0 {
                            (score, source)
                        } else {
                            acc
                        }
                    },
                );

                viterbi_matrix.set_sparse(
                    sparse_row_to_idx,
                    col_to_idx,
                    score
                        + confidence_matrix
                            .get_sparse(sparse_row_to_idx, col_to_idx)
                            .ln(),
                );

                sources_matrix.set_sparse(sparse_row_to_idx, col_to_idx, source)
            });
    }
}

#[derive(Default, Clone, Debug)]
pub struct TraceStep {
    pub row_idx: usize,
    pub sparse_row_idx: usize,
    pub query_id: usize,
    pub join_id: usize,
}

pub type Trace = Vec<TraceStep>;

pub fn traceback(
    viterbi_matrix: &Matrix<f64>,
    sources: &Matrix<usize>,
    active_cols: &[usize],
    join_id_start: usize,
) -> Trace {
    let last_col_idx = *active_cols
        .last()
        .expect("active_cols is empty in call to traceback()");

    let (_, sparse_row_idx_of_max_score) = viterbi_matrix
        .col_slice(last_col_idx)
        .iter()
        .enumerate()
        .fold((-f64::INFINITY, 0usize), |acc, (idx, &score)| {
            if score > acc.0 {
                (score, idx)
            } else {
                acc
            }
        });

    let mut simple_trace = vec![0; viterbi_matrix.num_cols()];
    simple_trace[last_col_idx] = sparse_row_idx_of_max_score;

    //  say we have sources like:
    //
    //    |  0  |  1  |  2  |  3  |  4  |
    //    |-----|-----|-----|-----|-----|
    //  0 |  0  |  0  |  0  |  0  |  0  |
    //  1 |  1  |  1  |  1  |  1  |  1  |
    //  2 |  2  |  1  |  1  |  2  |  2  |
    //
    //  where the value in a cell in column <i> is the sparse
    //  row index of the cell in column <i - 1> that is
    //  the source of the value of the cell in column <i>
    //
    //  if the max value at the end was in
    //  row 2, then the trace starts as:
    //
    //    [     ,     ,     ,     ,  2  ]
    //
    //  prev_col_idx: *, is the column we are labeling
    //       col_idx: $, is already labeled
    //
    //                         *     $
    //    [     ,     ,     ,  2  ,  2  ]
    //
    //                   *     $
    //    [     ,     ,  2  ,  2  ,  2  ]
    //
    //             *     $
    //    [     ,  1  ,  2  ,  2  ,  2  ]
    //
    //       *     $
    //    [  1  ,  1  ,  2  ,  2  ,  2  ]

    active_cols
        .iter()
        .zip(active_cols.iter().skip(1))
        .rev()
        .for_each(|(&prev_col_idx, &col_idx)| {
            // we have col_idx labeled by the previous
            // iteration, so we need to label prev_col_idx
            // with the source of the label in col_idx
            let col_idx_label = simple_trace[col_idx];
            // join transitions are encoded in the sources by
            // adding the sparse column length to the sparse
            // row in the sources matrix, so we need to mod by
            // col_length to get the right index
            let col_length = sources.col_length(col_idx);
            let prev_col_idx_label = sources.get_sparse(col_idx_label % col_length, col_idx);
            simple_trace[prev_col_idx] = prev_col_idx_label;
        });

    let mut trace: Trace = vec![TraceStep::default(); viterbi_matrix.num_cols()];
    let mut join_id_cnt = join_id_start;

    let first_col_idx = active_cols[0];
    let first_col_length = viterbi_matrix.col_length(first_col_idx);
    let first_sparse_row_idx = simple_trace[first_col_idx] % first_col_length;
    let first_col_query_id = viterbi_matrix.id_of_cell_sparse(first_sparse_row_idx, first_col_idx);

    let join_id = if first_col_query_id == 0 {
        0
    } else {
        join_id_cnt
    };

    trace[first_col_idx] = TraceStep {
        row_idx: viterbi_matrix.sparse_to_logical_row_idx(first_sparse_row_idx, first_col_idx),
        sparse_row_idx: first_sparse_row_idx,
        query_id: first_col_query_id,
        join_id,
    };

    active_cols
        .iter()
        .zip(active_cols.iter().skip(1))
        .for_each(|(&prev_col_idx, &col_idx)| {
            let prev_col_length = viterbi_matrix.col_length(prev_col_idx);
            let prev_sparse_row_idx = simple_trace[prev_col_idx];
            let prev_logical_row_idx = viterbi_matrix
                .sparse_to_logical_row_idx(prev_sparse_row_idx % prev_col_length, prev_col_idx);

            let col_length = viterbi_matrix.col_length(col_idx);
            let sparse_row_idx = simple_trace[col_idx] % col_length;
            let logical_row_idx = viterbi_matrix.sparse_to_logical_row_idx(sparse_row_idx, col_idx);
            let query_id = viterbi_matrix.id_of_cell_sparse(sparse_row_idx, col_idx);

            let prev_was_different = prev_logical_row_idx != logical_row_idx;
            let prev_was_join = prev_sparse_row_idx >= prev_col_length;

            if prev_was_different && !prev_was_join && query_id != 0 {
                join_id_cnt += 1;
            }

            let join_id = if query_id == 0 { 0 } else { join_id_cnt };

            trace[col_idx] = TraceStep {
                row_idx: logical_row_idx,
                sparse_row_idx,
                query_id,
                join_id,
            }
        });

    trace
}

#[allow(dead_code)]
pub fn print_viterbi_with_sources(viterbi_matrix: &Matrix<f64>, sources_matrix: &Matrix<usize>) {
    (0..viterbi_matrix.num_rows()).for_each(|row_idx| {
        (0..viterbi_matrix.num_cols()).for_each(|col_idx| {
            //
            if viterbi_matrix.contains_cell(row_idx, col_idx) {
                print!("{:12.3} ", viterbi_matrix.get(row_idx, col_idx));
            } else {
                print!("{:>12.3} ", "x");
            }
        });
        println!();

        (0..sources_matrix.num_cols()).for_each(|col_idx| {
            if sources_matrix.contains_cell(row_idx, col_idx) {
                print!("{:12.3} ", sources_matrix.get(row_idx, col_idx));
            } else {
                print!("{:>12.3} ", "x");
            }
        });
        println!();
        println!();
    });
}
