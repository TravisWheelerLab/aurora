use std::collections::HashMap;

use crate::matrix::Matrix;

///
///
///
///
///
pub fn windowed_confidence_slow(
    matrix: &mut Matrix<f64>,
) -> (HashMap<usize, f64>, HashMap<usize, Vec<f64>>) {
    // TODO: parameterize this
    let half_window_size = 15usize;

    let mut confidence_avg_by_id: HashMap<usize, f64> = HashMap::new();
    let mut confidence_by_id: HashMap<usize, Vec<f64>> = HashMap::new();

    // we keep a buffer of the computed windowed confidences
    // so that we don't overwrite a cell that needs
    // to be used in the computation for the next window
    let mut buffer = vec![0.0; matrix.num_cols()];

    (0..matrix.num_rows())
        // skip the skip state
        .skip(1)
        .for_each(|row_idx| {
            let (row_start, row_end) = matrix.col_range_by_row(row_idx);
            let ali_id = matrix.ali_id(row_idx, row_start);
            (row_start..=row_end).for_each(|center_of_window_col_idx| {
                let window_start_col_idx =
                // saturating subtract to prevent underflow
                (center_of_window_col_idx.saturating_sub(half_window_size))
                // max to prevent going past the beginning of an alignment
                .max(row_start);

                // min(row_end) to prevent going past the end of an alignment
                let window_end_col_idx = (center_of_window_col_idx + half_window_size).min(row_end);
                // the window size won't always be uniform, so we need to compute it every time
                let window_size = window_end_col_idx - window_start_col_idx + 1;

                let window_sum = (window_start_col_idx..=window_end_col_idx)
                    .fold(0.0, |acc, col_idx_in_window| {
                        acc + matrix.get(row_idx, col_idx_in_window)
                    });

                let window_avg = window_sum / window_size as f64;

                buffer[center_of_window_col_idx] = window_avg;
            });

            let mut confidence_sum = 0.0;

            let confidence_vec = confidence_by_id
                .entry(ali_id)
                .or_insert(Vec::with_capacity(row_end - row_start + 1));

            // once we've completed the row, we can copy the buffer into the matrix
            (row_start..=row_end).for_each(|col_idx| {
                matrix.set(row_idx, col_idx, buffer[col_idx]);
                confidence_vec.push(buffer[col_idx]);
                confidence_sum += buffer[col_idx];
            });

            confidence_avg_by_id.insert(ali_id, confidence_sum / (row_end - row_start + 1) as f64);
        });

    confidence_avg_by_id.remove(&0);
    confidence_by_id.remove(&0);

    (confidence_avg_by_id, confidence_by_id)
}
