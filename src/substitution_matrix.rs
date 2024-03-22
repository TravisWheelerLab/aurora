use std::{
    fs::File,
    io::{BufRead, BufReader, Read},
    path::Path,
};

use crate::alphabet::{
    ALIGNMENT_ALPHABET_STR, GAP_EXTEND_DIGITAL, GAP_OPEN_DIGITAL, NUCLEOTIDE_ALPHABET_UTF8,
    STR_TO_DIGITAL_NUCLEOTIDE,
};

pub trait AlignmentScore {
    fn score(&self, target_char: u8, query_char: u8) -> f64;
    fn score_with_background(&self, target_char: u8, query_char: u8, frequencies: &[f64; 4])
        -> f64;
    fn gap_open(&self) -> f64;
    fn gap_extend(&self) -> f64;
}

pub struct SimpleSubstitutionMatrix {
    pub match_score: f64,
    pub sub_score: f64,
    pub gap_open_score: f64,
    pub gap_extend_score: f64,
}

impl AlignmentScore for SimpleSubstitutionMatrix {
    fn score(&self, target_char: u8, query_char: u8) -> f64 {
        match (target_char, query_char) {
            (GAP_OPEN_DIGITAL, _) | (_, GAP_OPEN_DIGITAL) => self.gap_open_score,
            (GAP_EXTEND_DIGITAL, _) | (_, GAP_EXTEND_DIGITAL) => self.gap_extend_score,
            _ => {
                if target_char == query_char {
                    self.match_score
                } else {
                    self.sub_score
                }
            }
        }
    }

    #[allow(unused_variables)]
    fn score_with_background(
        &self,
        target_char: u8,
        query_char: u8,
        frequencies: &[f64; 4],
    ) -> f64 {
        self.score(target_char, query_char)
    }

    fn gap_open(&self) -> f64 {
        self.gap_open_score
    }

    fn gap_extend(&self) -> f64 {
        self.gap_extend_score
    }
}

pub struct SubstitutionMatrix {
    pub name: String,
    pub gap_open_score: f64,
    pub gap_extend_score: f64,
    pub scores: [[f64; 14]; 14],
    pub core_ratios: [[f64; 4]; 4],
}

impl std::fmt::Debug for SubstitutionMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "score matrix: {}", self.name)?;
        writeln!(f)?;

        writeln!(f, "core ratios (P(ab) / P(b)):")?;
        writeln!(f, "   {:6} {:6} {:6} {:6}", "A", "C", "G", "T")?;

        #[allow(clippy::needless_range_loop)]
        for i in 0..4 {
            write!(f, " {}", ALIGNMENT_ALPHABET_STR[i])?;
            for j in 0..4 {
                write!(f, " {:0.4}", self.core_ratios[i][j])?;
            }
            writeln!(f)?;
        }
        writeln!(f)?;

        // gap scores
        writeln!(f, "gap open:   {:?}", self.gap_open_score)?;
        writeln!(f, "gap extend: {:?}", self.gap_extend_score)?;
        writeln!(f)?;

        // unscaled scores
        writeln!(f, "unscaled scores ln(P(ab) / P(a)P(b)):")?;
        write!(f, "       ")?;
        #[allow(clippy::needless_range_loop)]
        for i in 0..12 {
            write!(f, "{:7} ", ALIGNMENT_ALPHABET_STR[i])?;
        }
        writeln!(f)?;

        #[allow(clippy::needless_range_loop)]
        for i in 0..12 {
            write!(f, " {}", ALIGNMENT_ALPHABET_STR[i])?;
            for j in 0..12 {
                write!(f, " {:7.4}", self.scores[i][j])?;
            }
            writeln!(f)?;
        }

        Ok(())
    }
}

impl PartialEq for SubstitutionMatrix {
    fn eq(&self, other: &Self) -> bool {
        self.name == other.name
    }
}

enum ParserState {
    Header,
    Chars,
    Scores,
}

impl AlignmentScore for SubstitutionMatrix {
    fn score_with_background(
        &self,
        target_char: u8,
        query_char: u8,
        frequencies: &[f64; 4],
    ) -> f64 {
        let t = target_char as usize;
        let q = query_char as usize;
        if t < 4 && q < 4 {
            (self.core_ratios[t][q] / frequencies[t]).ln()
        } else {
            self.scores[t][q]
        }
    }

    fn score(&self, target_char: u8, query_char: u8) -> f64 {
        self.scores[target_char as usize][query_char as usize]
    }

    fn gap_open(&self) -> f64 {
        self.gap_open_score
    }

    fn gap_extend(&self) -> f64 {
        self.gap_extend_score
    }
}

impl SubstitutionMatrix {
    pub fn new(
        name: &str,
        lambda: f64,
        gap_open: f64,
        gap_extend: f64,
        background_freqs_i: [f64; 4],
        target_background_frequencies: [f64; 4],
        original_scores: [[f64; 14]; 14],
    ) -> Self {
        let mut unscaled_scores = original_scores;

        unscaled_scores
            .iter_mut()
            .flat_map(|row| row.iter_mut())
            .for_each(|v| *v *= lambda);

        // multiply out the target background frequencies
        // for only the core nucleotides: [A, C, G, T]
        // **NOTE: we're going to leave the ambiguity characters alone
        //         because their alignmnent scores in Arian's matrices
        //         are essentially magic numbers that cannot be related
        //         to the core log odds ratios
        let mut core_ratios = [[0.0; 4]; 4];
        for target_char in 0..4 {
            for query_char in 0..4 {
                core_ratios[target_char][query_char] = unscaled_scores[target_char][query_char]
                    .exp()
                    * target_background_frequencies[target_char]
            }
        }

        Self {
            name: name.to_string(),
            gap_open_score: gap_open * lambda,
            gap_extend_score: gap_extend * lambda,
            scores: unscaled_scores,
            core_ratios,
        }
    }

    pub fn scores_with_background(&self, background: [f64; 4]) -> [[f64; 14]; 14] {
        let mut scores = self.scores;
        (0..4).for_each(|target_char| {
            (0..4).for_each(|query_char| {
                let score = &mut scores[target_char][query_char];
                *score = self.core_ratios[target_char][query_char] / background[target_char];
                *score = score.ln();
            });
        });

        scores
            .iter()
            .for_each(|row| row.iter().for_each(|s| debug_assert!(s.is_finite())));

        scores
    }

    pub fn gap_score(&self, gap_len: usize) -> f64 {
        let ext_len = gap_len.saturating_sub(1) as f64;
        self.gap_open_score + ext_len * self.gap_extend_score
    }

    #[allow(dead_code)]
    pub fn from_file(path: impl AsRef<Path>) -> Vec<SubstitutionMatrix> {
        let file = File::open(path).expect("failed to open matrix file");
        SubstitutionMatrix::parse(file)
    }

    pub fn parse<R: Read>(matrix_buf: R) -> Vec<SubstitutionMatrix> {
        let mut matrices = vec![];

        let buf_reader = BufReader::new(matrix_buf);
        let mut state = ParserState::Header;

        let mut lines: Vec<String> = buf_reader
            .lines()
            .map(|l| l.expect("failed to read line"))
            .filter(|l| !l.is_empty())
            .collect();

        // add a single blank line to the
        // end to serve as a sentinel
        lines.push("".to_string());

        let line_tokens: Vec<Vec<&str>> = lines
            .iter()
            .map(|l| l.split_whitespace().collect::<Vec<&str>>())
            .collect();

        let mut name = "".to_string();
        let mut lambda = 0.0;
        let mut gap_open = 0.0;
        let mut gap_extend = 0.0;
        let mut background_freqs_i = [0.0, 0.0, 0.0, 0.0];
        let mut background_freqs_j = [0.0, 0.0, 0.0, 0.0];
        let mut chars: Vec<String> = vec![];
        let mut scores_vec: Vec<Vec<f64>> = vec![];
        line_tokens
            .iter()
            .zip(line_tokens.iter().skip(1))
            .for_each(|(tokens, next_tokens)| {
                match state {
                    ParserState::Header => match tokens.first() {
                        Some(&"#matrix") => {
                            name = tokens[1].to_string();
                        }
                        Some(&"#gap-open") => {
                            gap_open = tokens[1].parse::<f64>().expect("failed to parse float");
                        }
                        Some(&"#gap-ext") => {
                            gap_extend = tokens[1].parse::<f64>().expect("failed to parse float");
                        }
                        Some(&"#lambda") => {
                            lambda = tokens[1].parse::<f64>().expect("failed to parse float");
                        }
                        Some(&"#fi") => {
                            let freqs: Vec<f64> = tokens[1..=4]
                                .iter()
                                .map(|&f| f.parse::<f64>().expect("failed to parse float"))
                                .collect();
                            background_freqs_i
                                .iter_mut()
                                .zip(freqs)
                                .for_each(|(a, b)| *a = b);
                        }
                        Some(&"#fj") => {
                            let freqs: Vec<f64> = tokens[1..=4]
                                .iter()
                                .map(|&f| f.parse::<f64>().expect("failed to parse float"))
                                .collect();
                            background_freqs_j
                                .iter_mut()
                                .zip(freqs)
                                .for_each(|(a, b)| *a = b);
                        }
                        _ => panic!(),
                    },
                    ParserState::Chars => {
                        chars = tokens.iter().map(|t| t.to_string()).collect();
                    }
                    ParserState::Scores => scores_vec.push(
                        tokens
                            .iter()
                            .map(|t| t.parse::<f64>().expect("failed to parse float"))
                            .collect(),
                    ),
                }

                // TODO: refactor this to get rid of the closure
                let mut add_matrix = || {
                    assert_eq!(scores_vec.len(), chars.len());
                    let mut scores = [[0.0; 14]; 14];

                    let char_indices: Vec<usize> = chars
                        .iter()
                        .map(|c| {
                            (*STR_TO_DIGITAL_NUCLEOTIDE.get(c).expect("invalid char")) as usize
                        })
                        .collect();

                    scores_vec
                        .iter()
                        .enumerate()
                        .map(|(row_in_file, r)| (char_indices[row_in_file], r))
                        .for_each(|(row_idx, row)| {
                            row.iter()
                                .enumerate()
                                .map(|(col_in_file, r)| (char_indices[col_in_file], r))
                                .for_each(|(col_idx, val)| {
                                    // NOTE: in the RepeatMasker cross_match
                                    //       matrices, the target is the
                                    //       columns and the query is the rows
                                    //
                                    // todo: need to have a command line flag to
                                    //       indicate what the orientation is

                                    scores[col_idx][row_idx] = *val;
                                });
                        });

                    matrices.push(SubstitutionMatrix::new(
                        &name,
                        lambda,
                        gap_open,
                        gap_extend,
                        background_freqs_i,
                        background_freqs_j,
                        scores,
                    ));
                    scores_vec = vec![];
                };

                state = match next_tokens.first() {
                    Some(token) if token.starts_with('#') => {
                        if let ParserState::Scores = state {
                            add_matrix()
                        }
                        ParserState::Header
                    }
                    Some(token) if token.parse::<f64>().is_err() => ParserState::Chars,
                    Some(token) if token.parse::<f64>().is_ok() => ParserState::Scores,
                    None => {
                        add_matrix();
                        return;
                    }
                    _ => panic!(),
                };
            });

        matrices
    }
}

#[cfg(test)]
mod tests {
    use crate::{
        alignment::Alignment,
        alphabet::NucleotideByteUtils,
        windowed_scores::{Background, BackgroundFrequencies},
    };

    use super::*;

    #[test]
    pub fn test_substitution_matrix_parse() -> anyhow::Result<()> {
        let matrix_buf = "
            #matrix 14p35g.matrix
            #gap-open -35
            #gap-ext -6
            #lambda 0.1284
            #fi 0.316249 0.183751 0.183751 0.316249
            #fj 0.326882 0.173118 0.173118 0.326882
              A   R   G   C   Y   T   K   M   S   W   N   X
              8   0 -10 -18 -19 -21 -15  -4 -14  -6  -1 -30
              3   3  12 -17 -18 -19  -9  -8  -8  -9  -1 -30
             -7   2  12 -16 -16 -17  -2 -11  -1 -12  -1 -30
            -17 -16 -16  12   2  -7 -11  -2  -1 -12  -1 -30
            -19 -18 -17  12   0   3  -8  -9  -8  -9  -1 -30
            -21 -19 -18 -10   0   8  -4 -15 -14  -6  -1 -30
            -14  -8  -2 -13  -8  -4  -3 -13  -8  -9  -1 -30
             -4  -8 -13  -2  -8 -14 -13  -3  -8  -9  -1 -30
            -12  -7  -1  -1  -7 -12  -7  -7  -1 -12  -1 -30
             -6 -10 -14 -14 -10  -6 -10 -10 -14  -6  -1 -30
             -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
            -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30"
            .as_bytes();

        // let lambda = 0.1284;

        //  ACGT scores with:
        //    rows = query
        //    cols = target
        //     A   C   G   T
        //  A   8 -18 -10 -21
        //  C -17  12 -16  -7
        //  G  -7 -16  12 -17
        //  T -21 -10 -18   8

        // ACGT scores reoriented with:
        //   rows = target
        //   cols = query
        // let core_scores = [
        //     [8.0, -17.0, -7.0, -21.0],
        //     [-18.0, 12.0, -16.0, -10.0],
        //     [-10.0, -16.0, 12.0, -18.0],
        //     [-21.0, -7.0, -17.0, 8.0],
        // ];

        // the core scores are exp(score)
        // let core_scores_unscaled: [[f64; 4]; 4] = [
        //     [1.0272, -2.1828, -0.8988, -2.6964],
        //     [-2.3112, 1.5408, -2.0544, -1.2840],
        //     [-1.2840, -2.0544, 1.5408, -2.3112],
        //     [-2.6964, -0.8988, -2.1828, 1.0272],
        // ];

        // fi = rows = query
        // fj = cols = target
        // let target_background_frequencies = [0.326882, 0.173118, 0.173118, 0.326882];
        // let query_background_frequencies = [0.316249, 0.183751, 0.183751, 0.316249];

        // the core ratios are P(ab) / P(a)P(b)
        // where:
        //   P(ab): target frequency: aligning residue a to residue b
        //   P(a):  background frequency of residue a (target)
        //   P(a):  background frequency of residue b (query)
        // let core_ratios = [
        //     [0.2888, 0.0068, 0.0244, 0.0070],
        //     [0.0054, 0.1485, 0.0041, 0.0152],
        //     [0.0152, 0.0041, 0.1485, 0.0054],
        //     [0.0070, 0.0244, 0.0068, 0.2888],
        // ];

        // the core ratios without target background are P(ab) / P(b)
        let core_ratios_without_target_background: [[f64; 4]; 4] = [
            [0.9131, 0.0368, 0.1331, 0.0220],
            [0.0172, 0.8082, 0.0222, 0.0479],
            [0.0479, 0.0222, 0.8082, 0.0172],
            [0.0220, 0.1331, 0.0368, 0.9131],
        ];

        let matrix_vec = SubstitutionMatrix::parse(matrix_buf);
        let matrix = matrix_vec.first().unwrap();

        matrix
            .core_ratios
            .iter()
            .zip(core_ratios_without_target_background)
            .for_each(|(a, b)| {
                a.iter()
                    .zip(b)
                    .for_each(|(&x, y)| assert_eq!(format!("{x:6.4}"), format!("{y:6.4}")));
            });

        // -35 * 0.1284 = -4.494
        assert_eq!(matrix.gap_open_score, -4.494);
        // -6  * 0.1284 = -0.7704
        assert_eq!(matrix.gap_extend_score, -0.7704);

        Ok(())
    }

    #[test]
    pub fn test_score() -> anyhow::Result<()> {
        let matrix_buf = "
            #matrix 14p35g.matrix
            #gap-open -35
            #gap-ext -6
            #lambda 0.1284
            #fi 0.316249 0.183751 0.183751 0.316249
            #fj 0.326882 0.173118 0.173118 0.326882
              A   R   G   C   Y   T   K   M   S   W   N   X
              8   0 -10 -18 -19 -21 -15  -4 -14  -6  -1 -30
              3   3  12 -17 -18 -19  -9  -8  -8  -9  -1 -30
             -7   2  12 -16 -16 -17  -2 -11  -1 -12  -1 -30
            -17 -16 -16  12   2  -7 -11  -2  -1 -12  -1 -30
            -19 -18 -17  12   0   3  -8  -9  -8  -9  -1 -30
            -21 -19 -18 -10   0   8  -4 -15 -14  -6  -1 -30
            -14  -8  -2 -13  -8  -4  -3 -13  -8  -9  -1 -30
             -4  -8 -13  -2  -8 -14 -13  -3  -8  -9  -1 -30
            -12  -7  -1  -1  -7 -12  -7  -7  -1 -12  -1 -30
             -6 -10 -14 -14 -10  -6 -10 -10 -14  -6  -1 -30
             -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
            -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30"
            .as_bytes();

        let matrix_vec = SubstitutionMatrix::parse(matrix_buf);
        let matrix = matrix_vec.first().unwrap();

        let correct: [[f64; 12]; 12] = [
            [
                1.0272, -2.1828, -0.8988, -2.6964, -1.7976, -0.5136, -0.1284, 0.3852, -1.5408,
                -0.7704, -3.8520, -2.4396,
            ],
            [
                -2.3112, 1.5408, -2.0544, -1.2840, -1.6692, -0.2568, -0.1284, -2.1828, -0.1284,
                -1.7976, -3.8520, 1.5408,
            ],
            [
                -1.2840, -2.0544, 1.5408, -2.3112, -0.2568, -1.6692, -0.1284, 1.5408, -0.1284,
                -1.7976, -3.8520, -2.1828,
            ],
            [
                -2.6964, -0.8988, -2.1828, 1.0272, -0.5136, -1.7976, -0.1284, -2.4396, -1.5408,
                -0.7704, -3.8520, 0.3852,
            ],
            [
                -1.9260, -1.4124, -0.2568, -0.5136, -0.3852, -1.6692, -0.1284, -1.1556, -0.8988,
                -1.2840, -3.8520, -1.0272,
            ],
            [
                -0.5136, -0.2568, -1.4124, -1.9260, -1.6692, -0.3852, -0.1284, -1.0272, -0.8988,
                -1.2840, -3.8520, -1.1556,
            ],
            [
                -0.1284, -0.1284, -0.1284, -0.1284, -0.1284, -0.1284, -0.1284, -0.1284, -0.1284,
                -0.1284, -3.8520, -0.1284,
            ],
            [
                0.0000, -2.0544, 0.2568, -2.4396, -1.0272, -1.0272, -0.1284, 0.3852, -0.8988,
                -1.2840, -3.8520, -2.3112,
            ],
            [
                -1.7976, -0.1284, -0.1284, -1.7976, -1.0272, -1.0272, -0.1284, -1.0272, -0.1284,
                -1.7976, -3.8520, -1.0272,
            ],
            [
                -0.7704, -1.5408, -1.5408, -0.7704, -1.1556, -1.1556, -0.1284, -1.1556, -1.5408,
                -0.7704, -3.8520, -1.1556,
            ],
            [
                -3.8520, -3.8520, -3.8520, -3.8520, -3.8520, -3.8520, -3.8520, -3.8520, -3.8520,
                -3.8520, -3.8520, -3.8520,
            ],
            [
                -2.4396, 0.2568, -2.0544, 0.0000, -1.0272, -1.0272, -0.1284, -2.3112, -0.8988,
                -1.2840, -3.8520, 0.0000,
            ],
        ];

        (0..12).for_each(|target_char| {
            (0..12).for_each(|query_char| {
                let score = matrix.score(target_char, query_char);
                assert_eq!(
                    format!("{:6.4}", score),
                    format!("{:6.4}", correct[target_char as usize][query_char as usize])
                );
            });
        });

        Ok(())
    }

    #[test]
    pub fn test_background_adjusted_score() -> anyhow::Result<()> {
        let matrix_buf = "
            #matrix 14p35g.matrix
            #gap-open -35
            #gap-ext -6
            #lambda 0.1284
            #fi 0.316249 0.183751 0.183751 0.316249
            #fj 0.326882 0.173118 0.173118 0.326882
              A   R   G   C   Y   T   K   M   S   W   N   X
              8   0 -10 -18 -19 -21 -15  -4 -14  -6  -1 -30
              3   3  12 -17 -18 -19  -9  -8  -8  -9  -1 -30
             -7   2  12 -16 -16 -17  -2 -11  -1 -12  -1 -30
            -17 -16 -16  12   2  -7 -11  -2  -1 -12  -1 -30
            -19 -18 -17  12   0   3  -8  -9  -8  -9  -1 -30
            -21 -19 -18 -10   0   8  -4 -15 -14  -6  -1 -30
            -14  -8  -2 -13  -8  -4  -3 -13  -8  -9  -1 -30
             -4  -8 -13  -2  -8 -14 -13  -3  -8  -9  -1 -30
            -12  -7  -1  -1  -7 -12  -7  -7  -1 -12  -1 -30
             -6 -10 -14 -14 -10  -6 -10 -10 -14  -6  -1 -30
             -1  -1  -1  -1  -1  -1  -1  -1  -1  -1  -1 -30
            -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30 -30"
            .as_bytes();

        let matrix_vec = SubstitutionMatrix::parse(matrix_buf);
        let matrix = matrix_vec.first().unwrap();

        let mut ali: Vec<Alignment> = [
            format!("{}\n{}", "AAAAA", "ACGTA"),
            format!("{}\n{}", "CCCCC", "ACGTA"),
            format!("{}\n{}", "GGGGG", "ACGTA"),
            format!("{}\n{}", "TTTTT", "ACGTA"),
        ]
        .iter()
        .map(|v| Alignment::from_str(v))
        .collect();

        let starts = [10, 15, 20, 25];
        let ends = [14, 19, 24, 29];
        (0..4).for_each(|i| {
            ali[i].target_start = starts[i];
            ali[i].target_end = ends[i];
        });

        let target_start = starts[0] - 5;
        let target_end = ends.last().unwrap() + 5;
        let target_length = target_end - target_start + 1;

        let background = Background::new(&ali, target_start, target_length, 5);

        let correct: [[f64; 5]; 4] = [
            // A-A     A-C      A-G      A-T     A-A
            [0.2657, -3.1384, -2.0170, -3.5914, 0.4199],
            // C-A     C-C      C-G      C-T     C-A
            [-3.5542, 0.0102, -3.8082, -2.8146, -3.5542],
            // G-A     G-C      G-G      G-T     G-A
            [-2.5270, -3.5850, -0.2130, -3.8418, -2.5270],
            // T-A     T-C      T-G      T-T     T-A
            [-3.3037, -1.7938, -3.3010, 0.0716, -3.4579],
        ];

        ali.iter().enumerate().for_each(|(ali_idx, ali)| {
            ali.target_seq
                .iter()
                .zip(&ali.query_seq)
                .enumerate()
                .map(|(idx, bytes)| (idx, idx + ali.target_start, bytes))
                .for_each(|(ali_pos, target_pos, (&target_char, &query_char))| {
                    println!(
                        "{target_pos}: {} - {}",
                        target_char.to_utf8_string(),
                        query_char.to_utf8_string()
                    );

                    let frequencies = background.frequencies_at(target_pos);

                    let score = matrix.score_with_background(target_char, query_char, &frequencies);
                    assert_eq!(
                        format!("{:6.4}", score),
                        format!("{:6.4}", correct[ali_idx][ali_pos])
                    );
                });
        });

        Ok(())
    }
}
