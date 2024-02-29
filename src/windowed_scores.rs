use anyhow::Context;
use thiserror::Error;

use crate::alignment::{Alignment, TandemRepeat, VecMap};
use crate::alphabet::{
    NucleotideByteUtils, GAP_EXTEND_DIGITAL, GAP_OPEN_DIGITAL, N_DIGITAL, PAD_DIGITAL,
};
use crate::matrix::Matrix;
use crate::substitution_matrix::{AlignmentScore, SubstitutionMatrix};

fn build_target_seq_from_alignments(
    alignments: &[Alignment],
    target_start: usize,
    target_length: usize,
) -> Vec<u8> {
    let mut target_seq = vec![N_DIGITAL; target_length];

    alignments.iter().for_each(|ali| {
        ali.target_seq
            .iter()
            // skip all of the target gaps
            .filter(|&&b| b != GAP_OPEN_DIGITAL && b != GAP_EXTEND_DIGITAL)
            // then if we enumerate, we get indices
            // relative to the target sequence
            .enumerate()
            // then we can do a little math to place those indicies
            // relative to the target byte vector we're building
            .map(|(idx, b)| (idx + ali.target_start - target_start, b))
            .for_each(|(idx, &byte)| target_seq[idx] = byte);
    });

    target_seq
}

pub struct Background {
    pub target_start: usize,
    pub target_end: usize,
    pub target_seq: Vec<u8>,
    pub frequencies: Vec<[f64; 4]>,
}

impl Background {
    pub fn from_alignments(alignments: &[Alignment]) -> Self {
        //
        //
        Self {
            target_start: todo!(),
            target_end: todo!(),
            target_seq: todo!(),
            frequencies: todo!(),
        }
    }
}

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

#[derive(Error, Debug)]
pub enum GapError {
    #[error("sequence is length 0")]
    ZeroLengthSequence,
    #[error("gap open found in already opened gap")]
    DoubleOpen,
    #[error("gap extend found before a gap open")]
    ExtendBeforeOpen,
    #[error("alignment starts with a target gap")]
    StartsWithGap,
    #[error("alignment ends with a target gap")]
    EndsWithGap,
}

fn locate_target_gaps(alignment: &Alignment) -> anyhow::Result<Vec<f64>> {
    let error_context = || {
        format!(
            "target seq: {}",
            alignment.target_seq.to_debug_utf8_string()
        )
    };

    let first_target_char = *alignment
        .target_seq
        .first()
        .ok_or(GapError::ZeroLengthSequence)?;

    if first_target_char == GAP_OPEN_DIGITAL || first_target_char == GAP_EXTEND_DIGITAL {
        return Err(GapError::StartsWithGap).with_context(error_context);
    }

    let target_length = alignment.target_end - alignment.target_start + 1;
    let mut gaps = vec![0.0; target_length - 1];
    let mut target_pos = 0usize;
    let mut gap_start: Option<usize> = None;

    alignment
        .target_seq
        .iter()
        .enumerate()
        .try_for_each(|(ali_pos, &target_byte)| {
            match (gap_start, target_byte) {
                (None, GAP_OPEN_DIGITAL) => {
                    gap_start = Some(ali_pos);
                }
                (None, GAP_EXTEND_DIGITAL) => {
                    return Err(GapError::ExtendBeforeOpen).with_context(error_context)
                }
                (None, _) => target_pos += 1,
                (Some(_), GAP_EXTEND_DIGITAL) => {}
                (Some(_), GAP_OPEN_DIGITAL) => {
                    return Err(GapError::DoubleOpen).with_context(error_context)
                }
                (Some(start), _) => {
                    let size = ali_pos - start;
                    gaps[target_pos - 1] = size as f64;
                    target_pos += 1;
                    gap_start = None;
                }
            };
            Ok(())
        })?;

    if gap_start.is_some() {
        Err(GapError::EndsWithGap).with_context(error_context)
    } else {
        Ok(gaps)
    }
}

fn locate_query_gaps(alignment: &Alignment) -> anyhow::Result<Vec<f64>> {
    let error_context = || format!("query seq: {}", alignment.target_seq.to_debug_utf8_string());

    let first_query_char = *alignment
        .query_seq
        .first()
        .ok_or(GapError::ZeroLengthSequence)?;

    if first_query_char == GAP_OPEN_DIGITAL || first_query_char == GAP_EXTEND_DIGITAL {
        return Err(GapError::StartsWithGap).with_context(error_context);
    }

    let target_length = alignment.target_end - alignment.target_start + 1;
    let mut gaps = vec![0.0; target_length];
    let mut gap_start: Option<usize> = None;
    alignment
        .query_seq
        .iter()
        .zip(&alignment.target_seq)
        .filter(|(_, &target_byte)| {
            target_byte != GAP_OPEN_DIGITAL && target_byte != GAP_EXTEND_DIGITAL
        })
        .enumerate()
        .try_for_each(|(target_pos, (&query_byte, _))| {
            match (gap_start, query_byte) {
                (None, GAP_OPEN_DIGITAL) => {
                    gap_start = Some(target_pos);
                }
                (None, GAP_EXTEND_DIGITAL) => {
                    return Err(GapError::ExtendBeforeOpen).with_context(error_context)
                }
                (None, _) => {}
                (Some(_), GAP_OPEN_DIGITAL) => {
                    return Err(GapError::DoubleOpen).with_context(error_context)
                }
                (Some(_), GAP_EXTEND_DIGITAL) => {}
                (Some(start), _) => {
                    let size = target_pos - start;
                    gaps[start..target_pos]
                        .iter_mut()
                        .for_each(|g| *g = size as f64);
                    gap_start = None;
                }
            }
            Ok(())
        })?;

    if gap_start.is_some() {
        Err(GapError::EndsWithGap).with_context(error_context)
    } else {
        Ok(gaps)
    }
}

#[allow(dead_code)]
pub fn windowed_score_2(
    alignment: &Alignment,
    substitution_matrix: &impl AlignmentScore,
    window_size: usize,
    background_window_size: usize,
) -> anyhow::Result<Vec<f64>> {
    let target_length = alignment.target_end - alignment.target_start + 1;

    let mut target_gaps = locate_target_gaps(alignment).unwrap();
    let mut query_gaps = locate_query_gaps(alignment).unwrap();

    debug_assert_eq!(target_length - 1, target_gaps.len());
    debug_assert_eq!(target_length, query_gaps.len());

    let target_gap_fn = |gap_length: &f64| {
        if *gap_length == 0.0 {
            0.0
        } else {
            substitution_matrix.gap_open() + (gap_length - 1.0) * substitution_matrix.gap_extend()
        }
    };

    let query_gap_fn = |gap_length: &f64| {
        if *gap_length == 0.0 {
            0.0
        } else {
            (substitution_matrix.gap_open() + (gap_length - 1.0) * substitution_matrix.gap_extend())
                / gap_length
        }
    };

    target_gaps.iter_mut().for_each(|g| *g = target_gap_fn(g));
    query_gaps.iter_mut().for_each(|g| *g = query_gap_fn(g));

    let scores_without_gaps: Vec<f64> = alignment
        .target_seq
        .iter()
        .zip(&alignment.query_seq)
        .filter(|(&t, _)| t != GAP_OPEN_DIGITAL && t != GAP_EXTEND_DIGITAL)
        .map(|(&t, &q)| match q {
            GAP_OPEN_DIGITAL | GAP_EXTEND_DIGITAL => 0.0,
            _ => substitution_matrix.score(t, q),
        })
        .collect();

    let mut scores = vec![0.0; target_length];

    let half_window = (window_size - 1) / 2;
    let mut left = 0usize;
    let mut right = (target_length - 1).min(half_window);

    scores[0] = scores_without_gaps[left..=right].iter().sum::<f64>()
        + target_gaps[left..right].iter().sum::<f64>()
        + query_gaps[left..=right].iter().sum::<f64>();

    (1..target_length).for_each(|idx| {
        left = idx.saturating_sub(half_window);
        right = (target_length - 1).min(idx + half_window);
        scores[idx] = scores[idx - 1];

        if left != 0 {
            scores[idx] -= scores_without_gaps[left - 1];
            scores[idx] -= target_gaps[left - 1];
            scores[idx] -= query_gaps[left - 1];
        }

        if right != target_length {
            scores[idx] += scores_without_gaps[right];
            scores[idx] += target_gaps[right - 1];
            scores[idx] += query_gaps[right];
        }
    });

    Ok(scores)
}

///
///
///
///
pub fn windowed_score(
    matrix: &mut Matrix<f64>,
    alignments: &[Alignment],
    tandem_repeats: &[TandemRepeat],
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
                // +1 for the skip state
                ali_idx + 1,
                alignments[ali_idx].target_start - matrix_target_start,
                &alignments[ali_idx],
                substitution_matrices.get(alignments[ali_idx].substitution_matrix_id),
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

            let ali_length = scores.len();
            let mut window_score = scores[0..half_window.min(ali_length)].iter().sum::<f64>()
                + target_gap_scores[0..half_window.min(ali_length)]
                    .iter()
                    .sum::<f64>();
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

    tandem_repeats
        .iter()
        .enumerate()
        .for_each(|(repeat_idx, repeat)| {
            let row_idx = repeat_idx + alignments.len() + 1;
            let col_start = repeat.target_start - matrix_target_start;
            let col_end = repeat.target_end - matrix_target_start;

            let repeat_len = repeat.target_end - repeat.target_start + 1;
            let mut scores_padded = vec![0.0; repeat_len + window_size];

            repeat.scores.iter().enumerate().for_each(|(idx, score)| {
                scores_padded[idx + half_window] = *score;
            });

            debug_assert_eq!(col_end - col_start + 1, repeat_len);
            let mut window_score = scores_padded[0..window_size].iter().sum::<f64>();
            (col_start + 1..=col_end)
                .zip(1..repeat_len)
                .map(|(col_idx, repeat_idx)| (col_idx, repeat_idx - 1, repeat_idx + window_size))
                .for_each(|(col_idx, left_idx, right_idx)| {
                    window_score -= scores_padded[left_idx];
                    window_score += scores_padded[right_idx];
                    matrix.set(row_idx, col_idx, window_score);
                });
        });
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::substitution_matrix::SimpleSubstitutionMatrix;

    #[test]
    pub fn test_locate_target_gaps() -> anyhow::Result<()> {
        let ali: Vec<Alignment> = [
            format!("{}\n{}", "AAAAA", "AAAAA"),
            format!("{}\n{}", "A-A-A-A-A", "AAAAAAAAA"),
            format!("{}\n{}", "A-A-+++A-A", "AAAAAAAAAA"),
            format!("{}\n{}", "A-AAAAAA-A", "AAA-+++AAA"),
        ]
        .iter()
        .map(|v| Alignment::from_str(v))
        .collect();

        let gaps: Vec<Vec<f64>> = vec![
            vec![0.0, 0.0, 0.0, 0.0],
            vec![1.0, 1.0, 1.0, 1.0],
            vec![1.0, 4.0, 1.0],
            vec![1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
        ];

        ali.iter().zip(gaps).try_for_each(|(a, g)| {
            assert_eq!(locate_target_gaps(a)?, g);
            Ok(())
        })
    }

    #[test]
    pub fn test_locate_query_gaps() -> anyhow::Result<()> {
        let ali: Vec<Alignment> = [
            format!("{}\n{}", "AAAAA", "AAAAA"),
            format!("{}\n{}", "AAAAAAAAA", "A-A-A-A-A",),
            format!("{}\n{}", "AAAAAAAAAA", "A-A-+++A-A"),
            format!("{}\n{}", "AAA-+++AAA", "A-AAAAAA-A"),
        ]
        .iter()
        .map(|v| Alignment::from_str(v))
        .collect();

        let gaps: Vec<Vec<f64>> = vec![
            vec![0.0, 0.0, 0.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0],
            vec![0.0, 1.0, 0.0, 4.0, 4.0, 4.0, 4.0, 0.0, 1.0, 0.0],
            vec![0.0, 1.0, 0.0, 0.0, 1.0, 0.0],
        ];

        ali.iter().zip(gaps).try_for_each(|(a, g)| {
            assert_eq!(locate_query_gaps(a)?, g);
            Ok(())
        })
    }

    #[test]
    pub fn test_build_target_seq_from_alignments() {
        let mut ali: Vec<Alignment> = [
            format!("{}\n{}", "AAAAA", "AAAAA"),
            format!("{}\n{}", "CCCCC", "CCCCC"),
            format!("{}\n{}", "GGGGG", "GGGGG"),
            format!("{}\n{}", "TTTTT", "TTTTT"),
        ]
        .iter()
        .map(|v| Alignment::from_str(v))
        .collect();

        let starts = [10, 20, 30, 40];
        let ends = [14, 24, 34, 44];
        (0..4).for_each(|i| {
            ali[i].target_start = starts[i];
            ali[i].target_end = ends[i];
        });

        let target_start = *starts.first().unwrap();
        let target_end = ends.last().unwrap();
        let target_length = target_end - target_start + 1;
        let seq = build_target_seq_from_alignments(&ali, target_start, target_length);

        let correct = [
            0, 0, 0, 0, 0, 6, 6, 6, 6, 6, 1, 1, 1, 1, 1, 6, 6, 6, 6, 6, 2, 2, 2, 2, 2, 6, 6, 6, 6,
            6, 3, 3, 3, 3, 3,
        ];
        assert_eq!(seq, correct);
    }

    #[test]
    pub fn test_windowed_score() -> anyhow::Result<()> {
        let s = SimpleSubstitutionMatrix {
            match_score: 2.0,
            sub_score: 1.0,
            gap_open_score: -2.0,
            gap_extend_score: -1.0,
        };

        let a = Alignment::from_str(&format!("{}\n{}", "A-AAAAAA-A", "AAA-+++AAA"));
        let correct = vec![2.0, 0.75, -0.5, -3.75, -3.75, -0.5, 0.75, 2.0];
        let scores = windowed_score_2(&a, &s, 3, crate::BACKGROUND_WINDOW_SIZE)?;
        assert_eq!(scores, correct);

        let correct = vec![-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0];
        let scores = windowed_score_2(
            &a,
            &s,
            crate::SCORE_WINDOW_SIZE,
            crate::BACKGROUND_WINDOW_SIZE,
        )?;
        assert_eq!(scores, correct);

        Ok(())
    }
}
