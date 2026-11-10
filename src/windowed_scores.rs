use anyhow::Context;
use itertools::Itertools;
use thiserror::Error;

use crate::alignment::{Alignment, TandemRepeat};
use crate::alphabet::{
    NucleotideByteUtils, GAP_EXTEND_DIGITAL, GAP_OPEN_DIGITAL, NUCLEOTIDE_WEIGHTS, PAD_DIGITAL,
};
use crate::matrix::{Matrix, MatrixDef};
use crate::substitution_matrix::{AlignmentScore, SubstitutionMatrix};
use crate::util::VecMap;
use crate::viz::RegionAdjudicationSodaWriter;

pub fn build_target_seq_from_alignments(
    alignments: &[Alignment],
    target_start: usize,
    target_length: usize,
) -> Vec<u8> {
    let mut target_seq = vec![PAD_DIGITAL; target_length];

    alignments.iter().for_each(|ali| {
        ali.target_aligned_sequence()
            // skip all of the target gaps
            .filter(|&b| b != GAP_OPEN_DIGITAL && b != GAP_EXTEND_DIGITAL)
            // then if we enumerate, we get indices
            // relative to the target sequence
            .enumerate()
            // then we can do a little math to place those indicies
            // relative to the target byte vector we're building
            .map(|(idx, b)| (idx + ali.target_start - target_start, b))
            .for_each(|(idx, byte)| target_seq[idx] = byte);
    });

    target_seq
}

pub trait BackgroundFrequencies {
    fn frequencies_at(&self, target_pos: usize) -> [f64; 4];
}

impl BackgroundFrequencies for Background<'_> {
    fn frequencies_at(&self, target_pos: usize) -> [f64; 4] {
        debug_assert!(target_pos >= self.target_start);

        self.frequencies[target_pos - self.target_start]
    }
}

#[allow(dead_code)]
pub struct DummyBackground {}

impl BackgroundFrequencies for DummyBackground {
    #[allow(unused_variables)]
    fn frequencies_at(&self, target_pos: usize) -> [f64; 4] {
        [f64::NAN, f64::NAN, f64::NAN, f64::NAN]
    }
}

// TODO: move to a new file
// Allow extra properties for now...
#[allow(dead_code)]
pub struct Background<'a> {
    pub target_start: usize,
    pub target_end: usize,
    pub target_seq: &'a [u8],
    pub frequencies: Vec<[f64; 4]>,
}

impl<'a> Background<'a> {
    pub fn new(
        target_seq: &'a [u8],
        target_start: usize,
        target_length: usize,
        window_size: usize,
    ) -> Self {
        // these are the overall average frequencies
        // in the sequence, and they get used as the
        // frequencies for missing positions (* byte)
        let mut avg_frequencies = [0.0; 4];

        // first sum up the characters
        target_seq.iter().map(|&b| b as usize).for_each(|b| {
            avg_frequencies[0] += NUCLEOTIDE_WEIGHTS[b][0];
            avg_frequencies[1] += NUCLEOTIDE_WEIGHTS[b][1];
            avg_frequencies[2] += NUCLEOTIDE_WEIGHTS[b][2];
            avg_frequencies[3] += NUCLEOTIDE_WEIGHTS[b][3];
        });

        // take the average of the frequencies
        (0..4).for_each(|idx| avg_frequencies[idx] /= target_length as f64);

        let target_end = target_start + target_length - 1;
        let half_window = (window_size - 1) / 2;

        // for every position in the target,
        // this holds the frequencies within
        // the window centered at that position
        let mut frequencies_by_position = vec![[0.0; 4]; target_length];

        // the counts of each base for the current window
        let mut window_counts = [0.0; 4];

        // indices of the left & right window bounds
        let mut left = 0usize;
        let mut right = (target_length - 1).min(half_window);
        // size of the window
        let mut size = (right - left + 1) as f64;

        // fill the first window counts
        target_seq[left..=right]
            .iter()
            .map(|&b| b as usize)
            .for_each(|b| {
                window_counts[0] += NUCLEOTIDE_WEIGHTS[b][0];
                window_counts[1] += NUCLEOTIDE_WEIGHTS[b][1];
                window_counts[2] += NUCLEOTIDE_WEIGHTS[b][2];
                window_counts[3] += NUCLEOTIDE_WEIGHTS[b][3];
            });

        // fill the first position frequencies
        (0..4).for_each(|char_idx| {
            frequencies_by_position[0][char_idx] = window_counts[char_idx] / size
        });

        (1..target_length).for_each(|target_idx| {
            let prev_left = left;
            let prev_right = right;

            left = target_idx.saturating_sub(half_window);
            right = (target_length - 1).min(target_idx + half_window);
            size = (right - left + 1) as f64;

            // if we've moved far enough to the right
            // such that our left bound changes
            if left != prev_left {
                let removed_byte = target_seq[left - 1] as usize;
                window_counts[0] -= NUCLEOTIDE_WEIGHTS[removed_byte][0];
                window_counts[1] -= NUCLEOTIDE_WEIGHTS[removed_byte][1];
                window_counts[2] -= NUCLEOTIDE_WEIGHTS[removed_byte][2];
                window_counts[3] -= NUCLEOTIDE_WEIGHTS[removed_byte][3];
            }

            // if we haven't moved far enough to the right
            // such that our right bound doesn't change
            if right != prev_right {
                let added_byte = target_seq[right] as usize;
                window_counts[0] += NUCLEOTIDE_WEIGHTS[added_byte][0];
                window_counts[1] += NUCLEOTIDE_WEIGHTS[added_byte][1];
                window_counts[2] += NUCLEOTIDE_WEIGHTS[added_byte][2];
                window_counts[3] += NUCLEOTIDE_WEIGHTS[added_byte][3];
            }

            (0..4).for_each(|char_idx| {
                frequencies_by_position[target_idx][char_idx] = window_counts[char_idx] / size
            });
        });

        Self {
            target_start,
            target_end,
            target_seq,
            frequencies: frequencies_by_position,
        }
    }
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
            alignment
                .target_aligned_sequence()
                .collect_vec()
                .to_debug_utf8_string()
        )
    };

    let first_target_char = alignment
        .target_aligned_sequence()
        .next()
        .ok_or(GapError::ZeroLengthSequence)?;

    if first_target_char == GAP_OPEN_DIGITAL || first_target_char == GAP_EXTEND_DIGITAL {
        return Err(GapError::StartsWithGap).with_context(error_context);
    }

    let target_length = alignment.target_end - alignment.target_start + 1;
    let mut gaps = vec![0.0; target_length - 1];
    let mut target_pos = 0usize;
    let mut gap_start: Option<usize> = None;

    alignment
        .target_aligned_sequence()
        .enumerate()
        .try_for_each(|(ali_pos, target_byte)| {
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
    let error_context = || {
        format!(
            "query seq: {}",
            alignment
                .query_aligned_sequence()
                .collect_vec()
                .to_debug_utf8_string()
        )
    };

    let first_query_char = alignment
        .query_aligned_sequence()
        .next()
        .ok_or(GapError::ZeroLengthSequence)?;

    if first_query_char == GAP_OPEN_DIGITAL || first_query_char == GAP_EXTEND_DIGITAL {
        return Err(GapError::StartsWithGap).with_context(error_context);
    }

    let target_length = alignment.target_end - alignment.target_start + 1;
    let mut gaps = vec![0.0; target_length];
    let mut gap_start: Option<usize> = None;
    alignment
        .query_aligned_sequence()
        .zip(alignment.target_aligned_sequence())
        .filter(|&(_, target_byte)| {
            target_byte != GAP_OPEN_DIGITAL && target_byte != GAP_EXTEND_DIGITAL
        })
        .enumerate()
        .try_for_each(|(target_pos, (query_byte, _))| {
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

fn copy_to_matrix<T: std::fmt::Display + Copy + std::default::Default>(
    matrix: &mut Matrix<T>,
    row: usize,
    vals: &[T],
) {
    let (start, end) = matrix.col_range_by_row(row);

    debug_assert_eq!(end - start + 1, vals.len());

    (start..=end)
        .zip(vals.iter())
        .for_each(|(col_idx, &score)| {
            matrix.set(row, col_idx, score);
        });
}

pub fn dump_debug_windowed_score_info(
    matrix_def: &MatrixDef,
    viz: &mut RegionAdjudicationSodaWriter,
    alignments: &[Alignment],
    tandem_repeats: &[TandemRepeat],
    substitution_matrices: &VecMap<SubstitutionMatrix>,
    background: &impl BackgroundFrequencies,
    skip_state_percent_id: f64,
) -> anyhow::Result<()> {
    let mut matrix: Matrix<f64> = Matrix::new(matrix_def);
    let iter_alignments = || {
        alignments
            .iter()
            .enumerate()
            .map(|(idx, a)| (idx, a, substitution_matrices.get(a.substitution_matrix_id)))
            // map the vec enumeration index to a row index
            .map(|(ali_idx, a, m)| (ali_idx + 1, a, m))
    };

    for (row_idx, ali, sub_matrix) in iter_alignments() {
        let mut target_gaps = target_gap_scores(ali, sub_matrix)?;
        target_gaps.push(0.0);
        copy_to_matrix(&mut matrix, row_idx, &target_gaps);
    }
    viz.write_matrix("_target_gap_costs", &matrix)?;

    for (row_idx, ali, sub_matrix) in iter_alignments() {
        let target_gaps = query_gap_scores(ali, sub_matrix)?;
        copy_to_matrix(&mut matrix, row_idx, &target_gaps);
    }
    viz.write_matrix("_query_gap_costs", &matrix)?;

    // set the skip score...
    let skip_scores = get_skip_state_scores(
        &mut matrix,
        alignments,
        substitution_matrices,
        background,
        skip_state_percent_id,
    );
    copy_to_matrix(&mut matrix, 0, &skip_scores);

    for (row_idx, ali, sub_matrix) in iter_alignments() {
        let al_scores = aligned_position_scores(ali, sub_matrix, background);
        copy_to_matrix(&mut matrix, row_idx, &al_scores);
    }

    for (row_idx, repeat) in tandem_repeats
        .iter()
        .enumerate()
        // map the vec enumeration index to a row index
        .map(|(idx, r)| (idx + 1 + alignments.len(), r))
    {
        copy_to_matrix(&mut matrix, row_idx, &repeat.scores);
    }

    viz.write_matrix("_aligned_costs", &matrix)?;

    Ok(())
}

fn get_skip_state_scores(
    matrix: &mut Matrix<f64>,
    alignments: &[Alignment],
    substitution_matrices: &VecMap<SubstitutionMatrix>,
    background: &impl BackgroundFrequencies,
    skip_state_percent_id: f64,
) -> Vec<f64> {
    let mut skip_scores: Vec<f64> = vec![0.0; matrix.num_cols()];

    let matrix_target_start = matrix.target_start();

    for (_row_idx, ali, sub_matrix) in alignments
        .iter()
        .enumerate()
        .map(|(idx, a)| (idx, a, substitution_matrices.get(a.substitution_matrix_id)))
        // map the vec enumeration index to a row index
        .map(|(ali_idx, a, m)| (ali_idx + 1, a, m))
    {
        for (idx, nucl) in ali
            .target_aligned_sequence()
            .filter(|&v| !matches!(v, GAP_EXTEND_DIGITAL | GAP_OPEN_DIGITAL))
            .enumerate()
        {
            let bg = background.frequencies_at(ali.target_start + idx);
            let mat_col = ali.target_start + idx - matrix_target_start;
            let target_letter_distribution = NUCLEOTIDE_WEIGHTS[nucl as usize];
            let mut expected_score = 0.0;

            for (t_i, t_prob) in target_letter_distribution.iter().enumerate() {
                let non_match_sum: f64 = bg
                    .iter()
                    .enumerate()
                    .filter_map(|(i, p)| if i != t_i { Some(p) } else { None })
                    .sum();

                for (q_i, q_prob) in bg.iter().enumerate() {
                    let weight = if t_i == q_i {
                        skip_state_percent_id
                    } else {
                        (1.0 - skip_state_percent_id) * (q_prob / non_match_sum)
                    };

                    let pair_score = sub_matrix.score_with_background(t_i as u8, q_i as u8, &bg);

                    expected_score += t_prob * weight * pair_score;
                }
            }

            skip_scores[mat_col] = skip_scores[mat_col].max(expected_score);
        }
    }

    skip_scores
}

fn init_skip_state(
    matrix: &mut Matrix<f64>,
    alignments: &[Alignment],
    substitution_matrices: &VecMap<SubstitutionMatrix>,
    background: &impl BackgroundFrequencies,
    window_size: usize,
    skip_state_percent_id: f64,
    skip_state_shift: f64,
) {
    let skip_scores = get_skip_state_scores(
        matrix,
        alignments,
        substitution_matrices,
        background,
        skip_state_percent_id,
    );
    let mut skip_windowed_scores: Vec<f64> = vec![0.0; matrix.num_cols()];

    windowed_sum(&skip_scores, &mut skip_windowed_scores, window_size);

    for skip_score in skip_windowed_scores.iter_mut() {
        *skip_score += skip_state_shift;
    }
    copy_to_matrix(matrix, 0, &skip_windowed_scores);
}

pub fn windowed_score(
    matrix: &mut Matrix<f64>,
    alignments: &[Alignment],
    tandem_repeats: &[TandemRepeat],
    substitution_matrices: &VecMap<SubstitutionMatrix>,
    background: &impl BackgroundFrequencies,
    window_size: usize,
    skip_state_percent_id: f64,
    skip_state_shift: f64,
) -> anyhow::Result<()> {
    // set the skip score uniformly
    init_skip_state(
        matrix,
        alignments,
        substitution_matrices,
        background,
        window_size,
        skip_state_percent_id,
        skip_state_shift,
    );

    for (row_idx, ali, sub_matrix) in alignments
        .iter()
        .enumerate()
        .map(|(idx, a)| (idx, a, substitution_matrices.get(a.substitution_matrix_id)))
        // map the vec enumeration index to a row index
        .map(|(ali_idx, a, m)| (ali_idx + 1, a, m))
    {
        let windowed_score = windowed_score_alignment(ali, sub_matrix, window_size, background)?;
        copy_to_matrix(matrix, row_idx, &windowed_score);
    }

    for (row_idx, repeat) in tandem_repeats
        .iter()
        .enumerate()
        // map the vec enumeration index to a row index
        .map(|(idx, r)| (idx + 1 + alignments.len(), r))
    {
        let windowed_score = windowed_score_tandem_repeat(repeat, window_size)?;
        copy_to_matrix(matrix, row_idx, &windowed_score);
    }

    Ok(())
}

fn windowed_sum(src: &[f64], dst: &mut [f64], window_size: usize) {
    debug_assert_eq!(src.len(), dst.len());
    debug_assert!(src.len() >= 1);

    let half_window = (window_size - 1) / 2;

    let mut left = 0usize;
    let mut right = (dst.len() - 1).min(half_window);

    dst[0] = src[left..=right].iter().sum();

    (1..dst.len()).for_each(|idx| {
        let prev_left = left;
        let prev_right = right;

        left = idx.saturating_sub(half_window);
        right = (dst.len() - 1).min(idx + half_window);
        dst[idx] = dst[idx - 1];

        if left != prev_left {
            dst[idx] -= src[left - 1];
        }

        if right != prev_right {
            dst[idx] += src[right];
        }
    });
}

fn windowed_score_tandem_repeat(
    repeat: &TandemRepeat,
    window_size: usize,
) -> anyhow::Result<Vec<f64>> {
    let mut windowed_scores = vec![0.0; repeat.scores.len()];
    windowed_sum(&repeat.scores, &mut windowed_scores, window_size);

    Ok(windowed_scores)
}

fn target_gap_scores(
    alignment: &Alignment,
    substitution_matrix: &impl AlignmentScore,
) -> anyhow::Result<Vec<f64>> {
    let mut target_gaps = locate_target_gaps(alignment)?;

    fn target_gap_score(gap_length: &f64, substitution_matrix: &impl AlignmentScore) -> f64 {
        if *gap_length == 0.0 {
            0.0
        } else {
            substitution_matrix.gap_open() + (gap_length - 1.0) * substitution_matrix.gap_extend()
        }
    }

    target_gaps
        .iter_mut()
        .for_each(|g| *g = target_gap_score(g, substitution_matrix));

    Ok(target_gaps)
}

fn query_gap_scores(
    alignment: &Alignment,
    substitution_matrix: &impl AlignmentScore,
) -> anyhow::Result<Vec<f64>> {
    let mut query_gaps = locate_query_gaps(alignment)?;

    fn query_gap_score(gap_length: &f64, substitution_matrix: &impl AlignmentScore) -> f64 {
        if *gap_length == 0.0 {
            0.0
        } else {
            (substitution_matrix.gap_open() + (gap_length - 1.0) * substitution_matrix.gap_extend())
            // There is a query entry for every nucleotide in the gap, so we have to divide by the length...
                / gap_length
        }
    }

    query_gaps
        .iter_mut()
        .for_each(|g| *g = query_gap_score(g, substitution_matrix));

    Ok(query_gaps)
}

fn aligned_position_scores(
    alignment: &Alignment,
    substitution_matrix: &impl AlignmentScore,
    background: &impl BackgroundFrequencies,
) -> Vec<f64> {
    alignment
        .target_aligned_sequence()
        .zip(alignment.query_aligned_sequence())
        .filter(|&(t, _)| t != GAP_OPEN_DIGITAL && t != GAP_EXTEND_DIGITAL)
        .enumerate()
        .map(|(target_idx, (t, q))| match q {
            GAP_OPEN_DIGITAL | GAP_EXTEND_DIGITAL => 0.0,
            _ => {
                let target_pos = target_idx + alignment.target_start;
                let frequencies = background.frequencies_at(target_pos);
                substitution_matrix.score_with_background(t, q, &frequencies)
            }
        })
        .collect()
}

fn windowed_score_alignment(
    alignment: &Alignment,
    substitution_matrix: &impl AlignmentScore,
    window_size: usize,
    background: &impl BackgroundFrequencies,
) -> anyhow::Result<Vec<f64>> {
    let ali_length = alignment.target_end - alignment.target_start + 1;

    let target_gaps = target_gap_scores(alignment, substitution_matrix)?;
    let query_gaps = query_gap_scores(alignment, substitution_matrix)?;

    debug_assert_eq!(ali_length - 1, target_gaps.len());
    debug_assert_eq!(ali_length, query_gaps.len());

    let scores_without_gaps: Vec<f64> =
        aligned_position_scores(alignment, substitution_matrix, background);

    debug_assert_eq!(ali_length, scores_without_gaps.len());

    let mut windowed_scores = vec![0.0; ali_length];

    let half_window = (window_size - 1) / 2;

    let mut left = 0usize;
    let mut right = (ali_length - 1).min(half_window);

    windowed_scores[0] = scores_without_gaps[left..=right].iter().sum::<f64>()
        + target_gaps[left..right].iter().sum::<f64>()
        + query_gaps[left..=right].iter().sum::<f64>();

    (1..ali_length).for_each(|idx| {
        let prev_left = left;
        let prev_right = right;

        left = idx.saturating_sub(half_window);
        right = (ali_length - 1).min(idx + half_window);
        windowed_scores[idx] = windowed_scores[idx - 1];

        if left != prev_left {
            windowed_scores[idx] -= scores_without_gaps[left - 1];
            windowed_scores[idx] -= target_gaps[left - 1];
            windowed_scores[idx] -= query_gaps[left - 1];
        }

        if right != prev_right {
            windowed_scores[idx] += scores_without_gaps[right];
            windowed_scores[idx] += target_gaps[right - 1];
            windowed_scores[idx] += query_gaps[right];
        }
    });

    Ok(windowed_scores)
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
            println!("{:?}", locate_query_gaps(a));
            assert_eq!(locate_query_gaps(a)?, g);
            Ok(())
        })
    }

    #[test]
    pub fn test_build_target_seq_from_alignments() {
        let starts = [10, 20, 30, 40];
        let ali: Vec<Alignment> = [
            format!("{}\n{}", "AAAAA", "AAAAA"),
            format!("{}\n{}", "CCCCC", "CCCCC"),
            format!("{}\n{}", "GGGGG", "GGGGG"),
            format!("{}\n{}", "TTTTT", "TTTTT"),
        ]
        .iter()
        .zip(starts.iter())
        .map(|(v, &start)| Alignment::from_str_with_target_offset(v, start))
        .collect();

        let target_start = ali.first().unwrap().target_start;
        let target_end = ali.last().unwrap().target_end;
        let target_length = target_end - target_start + 1;
        let seq = build_target_seq_from_alignments(&ali, target_start, target_length);

        let correct = [
            0, 0, 0, 0, 0, 18, 18, 18, 18, 18, 1, 1, 1, 1, 1, 18, 18, 18, 18, 18, 2, 2, 2, 2, 2,
            18, 18, 18, 18, 18, 3, 3, 3, 3, 3,
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

        // use a dummy background that
        // won't adjust the score
        let b = DummyBackground {};

        let a = Alignment::from_str(&format!("{}\n{}", "A-AAAAAA-A", "AAA-+++AAA"));
        let correct = vec![2.0, 0.75, -0.5, -3.75, -3.75, -0.5, 0.75, 2.0];
        let scores = windowed_score_alignment(&a, &s, 3, &b)?;
        assert_eq!(scores, correct);

        let correct = vec![-1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0];
        let scores = windowed_score_alignment(&a, &s, 31, &b)?;
        assert_eq!(scores, correct);

        Ok(())
    }

    #[test]
    pub fn test_background() -> anyhow::Result<()> {
        let ali: Vec<Alignment> = [
            format!("{}\n{}", "AAAAA", "AAAAA"),
            format!("{}\n{}", "CCCCC", "CCCCC"),
            format!("{}\n{}", "GGGGG", "GGGGG"),
            format!("{}\n{}", "TTTTT", "TTTTT"),
        ]
        .iter()
        .zip([10, 15, 20, 25].iter())
        .map(|(v, &start)| Alignment::from_str_with_target_offset(v, start))
        .collect();

        let target_start = ali.first().unwrap().target_start - 5;
        let target_end = ali.last().unwrap().target_end + 5;
        let target_length = target_end - target_start + 1;

        let target_seq = build_target_seq_from_alignments(&ali, target_start, target_length);
        let background = Background::new(&target_seq, target_start, target_length, 5);

        // *****AAAAACCCCCGGGGGTTTTT*****
        let correct = vec![
            // ***
            [0.25, 0.25, 0.25, 0.25],
            // ****
            [0.25, 0.25, 0.25, 0.25],
            // *****
            [0.25, 0.25, 0.25, 0.25],
            // ****A
            [0.4, 0.2, 0.2, 0.2],
            // ***AA
            [0.55, 0.15, 0.15, 0.15],
            // **AAA
            [0.7, 0.1, 0.1, 0.1],
            // *AAAA
            [0.85, 0.05, 0.05, 0.05],
            // AAAAA
            [1.0, 0.0, 0.0, 0.0],
            // AAAAC
            [0.8, 0.2, 0.0, 0.0],
            // AAACC
            [0.6, 0.4, 0.0, 0.0],
            // AACCC
            [0.4, 0.6, 0.0, 0.0],
            // ACCCC
            [0.2, 0.8, 0.0, 0.0],
            // CCCCC
            [0.0, 1.0, 0.0, 0.0],
            // CCCCG
            [0.0, 0.8, 0.2, 0.0],
            // CCCGG
            [0.0, 0.6, 0.4, 0.0],
            // CCGGG
            [0.0, 0.4, 0.6, 0.0],
            // CGGGG
            [0.0, 0.2, 0.8, 0.0],
            // GGGGG
            [0.0, 0.0, 1.0, 0.0],
            // GGGGT
            [0.0, 0.0, 0.8, 0.2],
            // GGGTT
            [0.0, 0.0, 0.6, 0.4],
            // GGTTT
            [0.0, 0.0, 0.4, 0.6],
            // GTTTT
            [0.0, 0.0, 0.2, 0.8],
            // TTTTT
            [0.0, 0.0, 0.0, 1.0],
            // TTTT*
            [0.05, 0.05, 0.05, 0.85],
            // TTT**
            [0.1, 0.1, 0.1, 0.7],
            // TT***
            [0.15, 0.15, 0.15, 0.55],
            // T****
            [0.2, 0.2, 0.2, 0.4],
            // *****
            [0.25, 0.25, 0.25, 0.25],
            // ****
            [0.25, 0.25, 0.25, 0.25],
            // ***
            [0.25, 0.25, 0.25, 0.25],
        ];

        correct
            .into_iter()
            .zip(background.frequencies)
            .for_each(|(a, b)| assert_eq!(a, b));

        Ok(())
    }
}
