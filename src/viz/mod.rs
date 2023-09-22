mod bed;
mod block;

use bed::*;
use block::*;

use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::Path,
};

use itertools::Itertools;
use serde::Serialize;

use crate::{
    alignment::Alignment,
    alphabet::{
        NucleotideByteUtils, DASH_UTF8, GAP_EXTEND_DIGITAL, GAP_OPEN_DIGITAL,
        NUCLEOTIDE_ALPHABET_UTF8, PAD_DIGITAL, SPACE_UTF8,
    },
    matrix::Matrix,
    results::Annotation,
};

///
///
///
///
#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct AuroraSodaData {
    target_start: usize,
    target_end: usize,
    target_seq: String,
    aurora_ann: Vec<BlockGroup>,
    reference_ann: Vec<BlockGroup>,
    alignments_ann: Vec<String>,
}

impl Alignment {
    pub fn soda_string(&self) -> String {
        let mut matches: Vec<u8> = vec![];
        let mut substitutions: Vec<u8> = vec![];
        let mut inserts: Vec<u8> = vec![];

        self.target_seq
            .iter()
            .zip(&self.query_seq)
            .for_each(|(&t, &q)| {
                match t {
                    GAP_OPEN_DIGITAL | GAP_EXTEND_DIGITAL => {
                        //
                    }
                    _ => {
                        if t == q {
                            matches.push(NUCLEOTIDE_ALPHABET_UTF8[q as usize]);
                            substitutions.push(SPACE_UTF8);
                            inserts.push(SPACE_UTF8);
                        } else if q == GAP_OPEN_DIGITAL || q == GAP_EXTEND_DIGITAL {
                            matches.push(SPACE_UTF8);
                            substitutions.push(SPACE_UTF8);
                            inserts.push(DASH_UTF8);
                        } else {
                            matches.push(SPACE_UTF8);
                            substitutions.push(NUCLEOTIDE_ALPHABET_UTF8[q as usize]);
                            inserts.push(SPACE_UTF8);
                        }
                    }
                }
            });

        let match_string = String::from_utf8(matches).unwrap();
        let substitution_string = String::from_utf8(substitutions).unwrap();
        let insert_string = String::from_utf8(inserts).unwrap();
        debug_assert_eq!(self.target_end - self.target_start + 1, match_string.len());

        format!(
            "{},{},{},{},{}",
            match_string,
            substitution_string,
            insert_string,
            self.target_start,
            self.target_end + 1,
        )
    }
}

///
///
///
///
impl AuroraSodaData {
    pub fn new(
        confidence_matrix: &Matrix<f64>,
        alignments: &[Alignment],
        annotations: &[Annotation],
        reference_bed_path: impl AsRef<Path>,
    ) -> Self {
        let target_start = confidence_matrix.target_start();
        let target_length = confidence_matrix.num_cols();
        let target_end = target_start + target_length - 1;

        let mut target_seq_digital_bytes = vec![PAD_DIGITAL; target_length];
        alignments
            .iter()
            .flat_map(|a| {
                a.target_seq
                    .iter()
                    .filter(|&&c| c != GAP_OPEN_DIGITAL && c != GAP_EXTEND_DIGITAL)
                    .enumerate()
                    .map(|(idx, c)| (idx + a.target_start - target_start, c))
            })
            .for_each(|(col_idx, &char)| target_seq_digital_bytes[col_idx] = char);

        let target_seq = target_seq_digital_bytes.into_utf8_string();

        let unique_join_ids: Vec<usize> = annotations.iter().map(|a| a.join_id).unique().collect();

        let aurora_ann: Vec<BlockGroup> = unique_join_ids
            .iter()
            .map(|&id| {
                BlockGroup::from_joined_annotations(
                    &mut annotations
                        .iter()
                        .filter(|&a| a.join_id == id)
                        .collect::<Vec<&Annotation>>(),
                )
            })
            .collect();

        let mut overlapping_bed = vec![];
        let reference_bed_file = File::open(reference_bed_path).expect("failed to open");
        let reader = BufReader::new(reference_bed_file);
        for line in reader.lines() {
            match line {
                Ok(line) => {
                    let tokens: Vec<&str> = line.split_whitespace().collect();

                    let chrom = tokens[0];
                    let thick_start = tokens[6].parse::<usize>().expect("failed to parse usize");
                    let thick_end = tokens[7].parse::<usize>().expect("failed to parse usize");
                    if thick_start < target_end && thick_end > target_start {
                        overlapping_bed.push(BedRecord::from_tokens(&tokens));
                    }
                }
                Err(_) => panic!("failed to read line"),
            }
        }

        let reference_ann = overlapping_bed
            .iter()
            .map(BlockGroup::from_bed_record)
            .collect::<Vec<BlockGroup>>();

        let alignments_ann: Vec<String> = alignments.iter().map(|a| a.soda_string()).collect();

        Self {
            target_start,
            target_end,
            target_seq,
            aurora_ann,
            reference_ann,
            alignments_ann,
        }
    }
}
