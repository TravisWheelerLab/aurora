use crate::matrix::Matrix;

///
///
///
///
pub fn confidence(matrix: &mut Matrix<f64>) {
    for col_idx in 0..matrix.num_cols() {
        let mut col_total = 0.0;

        for score in matrix.col_slice_mut(col_idx).iter_mut() {
            *score = (*score).exp2();
            col_total += *score;
        }

        for score in matrix.col_slice_mut(col_idx).iter_mut() {
            *score /= col_total;
            debug_assert!(!score.is_nan());
        }
    }
}
