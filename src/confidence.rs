use crate::matrix::Matrix;

///
///
///
///
pub fn confidence(matrix: &mut Matrix<f64>) -> Vec<f64> {
    let mut skip_confidence_by_col = vec![0.0; matrix.num_cols()];

    for col_idx in 0..matrix.num_cols() {
        let skip_score = matrix.get_mut(0, col_idx);
        *skip_score = skip_score.exp2();

        let mut col_total = *skip_score;

        for score in matrix
            .col_slice_mut(col_idx)
            .iter_mut()
            // .skip(1) to skip the skip state
            .skip(1)
        {
            *score = (*score).exp2();
            col_total += *score;
        }

        let skip_score = matrix.get_mut(0, col_idx);
        *skip_score /= col_total;

        skip_confidence_by_col[col_idx] = *skip_score;

        for score in matrix.col_slice_mut(col_idx).iter_mut() {
            *score /= col_total;
            debug_assert!(!score.is_nan());
        }
    }

    skip_confidence_by_col
}
