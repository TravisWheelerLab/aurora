use std::{
    fs::File,
    io::{BufRead, BufReader, Read},
    path::Path,
};

use crate::alphabet::STR_TO_DIGITAL_NUCLEOTIDE;

pub trait SubstitutionMatrixSliceExt {
    fn retrieve(&self, matrix_name: &str) -> &SubstitutionMatrix;
}

impl SubstitutionMatrixSliceExt for &[SubstitutionMatrix] {
    fn retrieve(&self, matrix_name: &str) -> &SubstitutionMatrix {
        self.iter()
            .find(|m| m.name == matrix_name)
            .unwrap_or_else(|| panic!("no matrix: {}", matrix_name))
    }
}

pub struct SubstitutionMatrix {
    pub name: String,
    pub lambda: f64,
    pub background_freqs_i: [f64; 4],
    pub background_freqs_j: [f64; 4],
    pub gap_open: f64,
    pub gap_extend: f64,
    pub original_scores: [[f64; 14]; 14],
    pub unscaled_scores: [[f64; 14]; 14],
    pub core_ratios: [[f64; 4]; 4],
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

impl SubstitutionMatrix {
    pub fn new(
        name: &str,
        lambda: f64,
        gap_open: f64,
        gap_extend: f64,
        background_freqs_i: [f64; 4],
        background_freqs_j: [f64; 4],
        original_scores: [[f64; 14]; 14],
    ) -> Self {
        let mut unscaled_scores = original_scores;

        // unscale and exponentiate everything
        unscaled_scores
            .iter_mut()
            .flat_map(|row| row.iter_mut())
            .for_each(|v| *v *= lambda);

        // multiply out the target background frequencies
        // for only the core nucleotides: [A, C, G, T]
        let mut core_ratios = [[0.0; 4]; 4];
        (0..4).for_each(|target_char| {
            (0..4).for_each(|query_char| {
                core_ratios[target_char][query_char] =
                    unscaled_scores[target_char][query_char].exp() * background_freqs_j[target_char]
            });
        });

        Self {
            name: name.to_string(),
            lambda,
            background_freqs_i,
            background_freqs_j,
            gap_open: gap_open * lambda,
            gap_extend: gap_extend * lambda,
            original_scores,
            unscaled_scores,
            core_ratios,
        }
    }

    pub fn scores_with_background(&self, background: [f64; 4]) -> [[f64; 14]; 14] {
        let mut scores = self.unscaled_scores;
        (0..4).for_each(|target_char| {
            (0..4).for_each(|query_char| {
                let score = &mut scores[target_char][query_char];
                *score = self.core_ratios[target_char][query_char] / background[target_char];
                *score = score.ln();
            });
        });
        scores
    }

    pub fn gap_score(&self, gap_len: usize) -> f64 {
        let ext_len = gap_len.saturating_sub(1) as f64;
        self.gap_open + ext_len * self.gap_extend
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
                        .map(|(i, r)| (char_indices[i], r))
                        .for_each(|(row_idx, row)| {
                            row.iter()
                                .enumerate()
                                .map(|(i, r)| (char_indices[i], r))
                                .for_each(|(col_idx, val)| {
                                    scores[row_idx][col_idx] = *val;
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
