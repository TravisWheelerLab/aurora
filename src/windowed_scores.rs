use crate::alignment::{Alignment, VecMap};
use crate::alphabet::{GAP_EXTEND_DIGITAL, GAP_OPEN_DIGITAL, PAD_DIGITAL};
use crate::matrix::Matrix;
use crate::substitution_matrix::SubstitutionMatrix;

///
///
///
///
pub fn background(
    alignments: &[Alignment],
    target_length: usize,
    target_start: usize,
    window_size: usize,
) -> Vec<[f64; 4]> {
    //
    assert!((window_size - 1) % 2 == 0);
    let half_window = (window_size - 1) / 2;

    let mut target_seq = vec![PAD_DIGITAL; target_length + window_size - 1];
    alignments
        .iter()
        .flat_map(|a| {
            a.target_seq
                .iter()
                .filter(|&&c| c != GAP_OPEN_DIGITAL && c != GAP_EXTEND_DIGITAL)
                .enumerate()
                .map(|(idx, c)| (idx + a.target_start - target_start + half_window, c))
        })
        .for_each(|(col_idx, &char)| target_seq[col_idx] = char);

    let mut window_totals = [0.0; 4];
    target_seq[0..window_size].iter().for_each(|&c| match c {
        c if c < 4 => window_totals[c as usize] += 1.0,
        // TODO: actually be smart about ambiguity
        _ => window_totals.iter_mut().for_each(|v| *v += 0.25),
    });

    let mut background = vec![[0.0, 0.0, 0.0, 0.0]; target_length];
    (0..4).for_each(|i| background[0][i] = window_totals[i] / window_size as f64);

    let start = half_window + 1;
    let end = target_seq.len() - half_window;
    (start..end)
        .map(|idx| {
            (
                idx - half_window,
                target_seq[idx - half_window],
                target_seq[idx + half_window],
            )
        })
        .for_each(|(target_idx, removed_char, added_char)| {
            match removed_char {
                c if c < 4 => window_totals[c as usize] -= 1.0,
                // TODO: actually be smart about ambiguity
                _ => window_totals.iter_mut().for_each(|v| *v -= 0.25),
            }
            match added_char {
                c if c < 4 => window_totals[c as usize] += 1.0,
                // TODO: actually be smart about ambiguity
                _ => window_totals.iter_mut().for_each(|v| *v += 0.25),
            }
            (0..4).for_each(|i| background[target_idx][i] = window_totals[i] / window_size as f64);
        });

    background
        .iter()
        .map(|b| b.iter().sum::<f64>())
        .for_each(|s| debug_assert!((1.0 - s).abs() < 1e-6));

    background
}

///
///
///
///
pub fn windowed_score(
    matrix: &mut Matrix<f64>,
    alignments: &[Alignment],
    substitution_matrices: &VecMap<SubstitutionMatrix>,
    window_size: usize,
    background_window_size: usize,
) {
    assert!((window_size - 1) % 2 == 0);
    let half_window = (window_size - 1) / 2;
    let num_cols = matrix.num_cols();
    let matrix_target_start = matrix.target_start();

    let background = background(
        alignments,
        num_cols,
        matrix_target_start,
        background_window_size,
    );

    // set the skip score uniformly
    (0..num_cols).for_each(|col_idx| {
        // TODO: need a constant & parameter for this
        matrix.set(0, col_idx, 10.0);
    });

    (0..alignments.len())
        .map(|ali_idx| {
            (
                ali_idx + 1,
                alignments[ali_idx].target_start - matrix_target_start,
                &alignments[ali_idx],
                substitution_matrices.value(alignments[ali_idx].substitution_matrix_id),
            )
        })
        .for_each(|(row_idx, start_col_idx, alignment, substitution_matrix)| {
            let mut scores: Vec<f64> = vec![];
            let mut target_gap_scores: Vec<f64> = vec![];
            let mut target_gap_start: Option<usize> = None;
            let mut query_gap_start: Option<usize> = None;

            alignment
                .target_seq
                .iter()
                .enumerate()
                .zip(&alignment.query_seq)
                .for_each(|((ali_pos_idx, &target_char), &query_char)| {
                    let compact_idx = scores.len();
                    let target_idx = alignment.target_start + compact_idx - matrix_target_start;

                    let complexity_adjusted_scores =
                        substitution_matrix.scores_with_background(background[target_idx]);

                    match target_char {
                        GAP_OPEN_DIGITAL => target_gap_start = Some(ali_pos_idx),
                        GAP_EXTEND_DIGITAL => {}
                        _ => {
                            scores.push(
                                complexity_adjusted_scores[target_char as usize]
                                    [query_char as usize],
                            );
                            if let Some(idx) = target_gap_start {
                                let gap_len = ali_pos_idx - idx;
                                target_gap_scores.push(substitution_matrix.gap_score(gap_len));
                                target_gap_start = None;
                            } else {
                                target_gap_scores.push(0.0);
                            }
                        }
                    }

                    match query_char {
                        GAP_OPEN_DIGITAL => query_gap_start = Some(compact_idx),
                        GAP_EXTEND_DIGITAL => {}
                        _ => {
                            if let Some(idx) = query_gap_start {
                                let gap_len = compact_idx - idx;
                                let gap_score =
                                    substitution_matrix.gap_score(gap_len) / gap_len as f64;
                                scores[idx..compact_idx]
                                    .iter_mut()
                                    .for_each(|b| *b += gap_score);
                                query_gap_start = None;
                            }
                        }
                    }
                });

            let mut window_score = scores[0..half_window].iter().sum::<f64>()
                + target_gap_scores[0..half_window].iter().sum::<f64>();
            matrix.set(row_idx, start_col_idx, window_score);

            let len = scores.len();
            let mut prev_left_idx = 0usize;
            let mut prev_right_idx = half_window;
            (1..len)
                .map(|idx| {
                    (
                        idx + start_col_idx,
                        idx.saturating_sub(half_window),
                        (len - 1).min(idx + half_window),
                    )
                })
                .for_each(|(col_idx, left_idx, right_idx)| {
                    if left_idx != prev_left_idx {
                        window_score -= scores[left_idx - 1];
                        window_score -= target_gap_scores[left_idx];
                    }
                    if right_idx != prev_right_idx {
                        window_score += scores[right_idx];
                        window_score += target_gap_scores[right_idx];
                    }

                    matrix.set(row_idx, col_idx, window_score);
                    prev_left_idx = left_idx;
                    prev_right_idx = right_idx;
                });
        });
}
