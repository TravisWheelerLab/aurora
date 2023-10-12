mod bed;
mod block;

use bed::*;
use block::*;

use std::{
    fs::File,
    io::{BufRead, BufReader},
};

use itertools::Itertools;
use serde::Serialize;

use crate::{
    alignment::Alignment,
    alphabet::{
        NucleotideByteUtils, ALIGNMENT_ALPHABET_UTF8, GAP_EXTEND_DIGITAL, GAP_OPEN_DIGITAL,
        PAD_DIGITAL, SPACE_UTF8,
    },
    matrix::{Matrix, MatrixDef},
    results::Annotation,
    segments::Segments,
    Args,
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
    alignment_strings: Vec<String>,
    trace_strings: Vec<String>,
    fragment_strings: Vec<String>,
    segment_strings: Vec<String>,
}

impl Alignment {
    pub fn soda_string(&self, row: usize) -> String {
        let mut green_bytes: Vec<u8> = vec![];
        let mut orange_bytes: Vec<u8> = vec![];

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
                            green_bytes.push(ALIGNMENT_ALPHABET_UTF8[q as usize]);
                            orange_bytes.push(SPACE_UTF8);
                        } else {
                            green_bytes.push(SPACE_UTF8);
                            orange_bytes.push(ALIGNMENT_ALPHABET_UTF8[q as usize]);
                        }
                    }
                }
            });

        let green_string = String::from_utf8(green_bytes).unwrap();
        let orange_string = String::from_utf8(orange_bytes).unwrap();
        debug_assert_eq!(self.target_end - self.target_start + 1, green_string.len());

        format!(
            "{},{},{},{},{},{},{},{}",
            green_string,
            orange_string,
            self.target_start,
            self.target_end + 1,
            todo!(),
            row,
            self.query_id,
            self.strand,
        )
    }
}

impl Segments {
    pub fn trace_soda_string(&self, matrix_def: &MatrixDef) -> String {
        (0..self.num_segments)
            // only look at the segments removed
            .filter(|&idx| self.marked_for_removal[idx])
            // (segments have multiple fragments if there was a join)
            .flat_map(|idx| &self.trace_fragments[idx])
            .map(|f| {
                format!(
                    "{},{},{},{}, {}",
                    f.start_col_idx,
                    f.end_col_idx,
                    f.query_id,
                    f.row_idx,
                    matrix_def.query_names[f.query_id]
                )
            })
            .join("|")
    }

    pub fn fragments_soda_string(&self, matrix_def: &MatrixDef) -> String {
        (0..self.num_segments)
            .filter(|&idx| self.marked_for_removal[idx])
            .flat_map(|idx| &self.fragments[idx])
            .map(|f| {
                format!(
                    "{},{},{},{:5.4},{},{},{},{}",
                    f.start_col_idx,
                    f.end_col_idx,
                    f.row_idx,
                    f.avg_confidence,
                    matrix_def.query_names[f.query_id],
                    f.consensus_start,
                    f.consensus_end,
                    f.strand
                )
            })
            .join("|")
    }

    pub fn segments_soda_string(&self) -> String {
        (0..self.num_segments)
            .filter(|&idx| self.marked_for_removal[idx])
            .map(|idx| &self.segments[idx])
            .map(|s| format!("{},{}", s.start_col_idx, s.end_col_idx))
            .join("|")
    }
}

///
///
///
///
impl AuroraSodaData {
    // TODO: these parameters need a refactor
    pub fn new(
        confidence_matrix: &Matrix<f64>,
        alignments: &[Alignment],
        annotations: &[Annotation],
        trace_strings: Vec<String>,
        fragment_strings: Vec<String>,
        segment_strings: Vec<String>,
        args: &Args,
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

        let target_name = "";
        if let (Some(path), Some(&offset)) = (
            &args.viz_reference_bed_path,
            args.viz_reference_bed_index.get(target_name),
        ) {
            let file = File::open(path).expect("failed to open reference bed");
            let reader = BufReader::new(file);

            reader
                .lines()
                .skip(offset)
                .map(|l| l.expect("failed to read line"))
                .for_each(|line| {
                    let tokens: Vec<&str> = line.split_whitespace().collect();

                    let target = tokens[0];
                    if target != target_name {
                        return;
                    }

                    let thick_start = tokens[6].parse::<usize>().expect("failed to parse usize");
                    let thick_end = tokens[7].parse::<usize>().expect("failed to parse usize");
                    if thick_start < target_end && thick_end > target_start {
                        overlapping_bed.push(BedRecord::from_tokens(&tokens));
                    }
                });
        }

        let reference_ann = overlapping_bed
            .iter()
            .map(BlockGroup::from_bed_record)
            .collect::<Vec<BlockGroup>>();

        let alignment_strings: Vec<String> = alignments
            .iter()
            .enumerate()
            .map(|(idx, a)| a.soda_string(idx + 1))
            .collect();

        Self {
            target_start,
            target_end,
            target_seq,
            aurora_ann,
            reference_ann,
            alignment_strings,
            trace_strings,
            fragment_strings,
            segment_strings,
        }
    }
}
