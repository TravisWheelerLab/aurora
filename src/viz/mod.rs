mod bed;
mod block;

use bed::*;
use block::*;

use std::{
    collections::HashMap,
    fs::{self, File},
    io::{BufRead, BufReader},
    path::Path,
};

use itertools::Itertools;
use serde::Serialize;

use crate::{
    alignment::{Alignment, Strand, VecMap},
    alphabet::{
        NucleotideByteUtils, ALIGNMENT_ALPHABET_UTF8, GAP_EXTEND_DIGITAL, GAP_OPEN_DIGITAL,
        PAD_DIGITAL, SPACE_UTF8,
    },
    collapse::{Assembly, AssemblyGroup},
    results::Annotation,
    viterbi::TraceSegment,
    Args,
};

pub fn write_soda_html(
    data: &impl Serialize,
    template_path: impl AsRef<Path>,
    js_path: impl AsRef<Path>,
    out_path: impl AsRef<Path>,
) {
    let template = fs::read_to_string(template_path).expect("failed to read template");
    let js = fs::read_to_string(js_path).expect("failed to read js");

    let mut viz_html = template.replace(
        "DATA_TARGET",
        &serde_json::to_string(data).expect("failed to serialize JSON data"),
    );

    viz_html = viz_html.replace("JS_TARGET", &js);

    let mut file = std::fs::File::create(out_path).expect("failed to create file");

    std::io::Write::write_all(&mut file, viz_html.as_bytes()).expect("failed to write to file");
}

impl Alignment {
    pub fn soda_string(&self, row: usize, query_name: &str) -> String {
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
            query_name,
            row,
            self.query_id,
            self.strand,
        )
    }
}

///
///
///
///
#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct AdjudicationSodaData2 {
    target_start: usize,
    target_end: usize,
    target_seq: String,
    aurora_ann: Vec<BlockGroup>,
    reference_ann: Vec<BlockGroup>,
    alignment_strings: Vec<String>,
    assembly_strings: Vec<String>,
    tandem_repeat_strings: Vec<String>,
    conclusive_trace_strings: Vec<String>,
    ambiguous_trace_strings: Vec<String>,
    resolved_assembly_rows: Vec<Vec<usize>>,
    unresolved_assembly_rows: Vec<Vec<usize>>,
    competed_assembly_rows: Vec<Vec<usize>>,
}

impl AdjudicationSodaData2 {
    // TODO: refactor these args
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        group: &AssemblyGroup,
        query_names: &VecMap<String>,
        annotations: &[Annotation],
        trace_conclusive: Vec<Vec<TraceSegment>>,
        trace_ambiguous: Vec<Vec<TraceSegment>>,
        resolved_assembly_rows: Vec<Vec<usize>>,
        unresolved_assembly_rows: Vec<Vec<usize>>,
        competed_assembly_rows: Vec<Vec<usize>>,
        args: &Args,
    ) -> Self {
        let target_start = group.target_start;
        let target_end = group.target_end;
        let target_length = target_end - target_start + 1;

        let mut target_seq_digital_bytes = vec![PAD_DIGITAL; target_length];
        group
            .assemblies
            .iter()
            .flat_map(|assembly| &assembly.alignments)
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
            .collect_vec();

        let alignment_strings = group
            .assemblies
            .iter()
            .enumerate()
            .flat_map(|(idx, assembly)| {
                assembly
                    .alignments
                    .iter()
                    .map(move |a| a.soda_string(idx + 1, query_names.get(a.query_id)))
            })
            .collect_vec();

        let assembly_strings = group
            .assemblies
            .iter()
            .enumerate()
            .map(|(assembly_idx, assembly)| {
                format!(
                    "{},{},{},{},{}",
                    assembly.target_start,
                    assembly.target_end,
                    assembly.query_id,
                    assembly.alignments.len(),
                    assembly_idx + 1,
                )
            })
            .collect_vec();

        let tandem_repeat_strings = group
            .tandem_repeats
            .iter()
            .enumerate()
            .map(|(repeat_idx, repeat)| {
                format!(
                    "{},{},{},{}",
                    repeat.target_start,
                    repeat.target_end,
                    repeat.consensus_pattern,
                    repeat_idx + group.assemblies.len() + 1,
                )
            })
            .collect_vec();

        let conclusive_trace_strings = trace_conclusive
            .iter()
            .map(|t| {
                t.iter()
                    .map(|seg| {
                        format!(
                            "{},{},{},{}",
                            seg.col_start, seg.col_end, seg.query_id, seg.row_idx
                        )
                    })
                    .join("|")
            })
            .collect_vec();

        let ambiguous_trace_strings = trace_ambiguous
            .iter()
            .map(|t| {
                t.iter()
                    .map(|seg| {
                        format!(
                            "{},{},{},{}",
                            seg.col_start, seg.col_end, seg.query_id, seg.row_idx
                        )
                    })
                    .join("|")
            })
            .collect_vec();

        Self {
            target_start,
            target_end,
            target_seq,
            aurora_ann,
            reference_ann,
            alignment_strings,
            assembly_strings,
            tandem_repeat_strings,
            conclusive_trace_strings,
            ambiguous_trace_strings,
            resolved_assembly_rows,
            unresolved_assembly_rows,
            competed_assembly_rows,
        }
    }
}

///
///
///
///
#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct AssemblySodaData {
    target_start: usize,
    target_end: usize,
    consensus_start: usize,
    consensus_end: usize,
    consensus_ali_strings: Vec<String>,
    target_ali_strings: Vec<String>,
    consensus_assembly_strings: Vec<Vec<String>>,
    target_assembly_strings: Vec<Vec<String>>,
    links: Vec<String>,
    prev: usize,
    next: usize,
    suffix: String,
}

impl AssemblySodaData {
    pub fn new(
        assemblies: &[Assembly],
        query_ids: &[usize],
        links: Vec<String>,
        confidence: &HashMap<usize, f64>,
    ) -> Self {
        let query_id = assemblies[0].query_id;
        let strand = assemblies[0].strand;

        let target_assembly_strings = assemblies
            .iter()
            .map(|a| {
                a.alignments
                    .iter()
                    .map(|ali| {
                        format!(
                            "{},{},{},{}",
                            ali.id,
                            ali.target_start,
                            ali.target_end,
                            confidence.get(&ali.id).unwrap(),
                        )
                    })
                    .collect_vec()
            })
            .collect_vec();

        let consensus_assembly_strings = assemblies
            .iter()
            .map(|a| {
                a.alignments
                    .iter()
                    .map(|ali| match strand {
                        Strand::Forward => format!(
                            "{},{},{},{}",
                            ali.id,
                            ali.query_start,
                            ali.query_end,
                            confidence.get(&ali.id).unwrap(),
                        ),

                        Strand::Reverse => format!(
                            "{},{},{},{}",
                            ali.id,
                            ali.query_end,
                            ali.query_start,
                            confidence.get(&ali.id).unwrap(),
                        ),

                        Strand::Unset => panic!(),
                    })
                    .collect_vec()
            })
            .collect_vec();

        let target_ali_strings = target_assembly_strings.iter().flatten().cloned().collect();

        let consensus_ali_strings = consensus_assembly_strings
            .iter()
            .flatten()
            .cloned()
            .collect();

        let query_id_idx = query_ids
            .iter()
            .position(|id| *id == query_id)
            .expect("failed to find query_id");

        let next_idx = query_id_idx + 1;

        let next = if next_idx < query_ids.len() {
            query_ids[next_idx]
        } else {
            query_ids[0]
        };

        let prev_idx = query_id_idx as isize - 1;
        let prev = if prev_idx > 0 {
            query_ids[prev_idx as usize]
        } else {
            query_ids[query_ids.len() - 1]
        };

        let suffix = match strand {
            Strand::Forward => "fwd".to_string(),
            Strand::Reverse => "rev".to_string(),
            _ => panic!(),
        };

        let target_start = assemblies.iter().map(|a| a.target_start).min().unwrap();
        let target_end = assemblies.iter().map(|a| a.target_end).max().unwrap();

        let (consensus_start, consensus_end) = match strand {
            Strand::Forward => (
                assemblies.iter().map(|a| a.query_start).min().unwrap(),
                assemblies.iter().map(|a| a.query_end).max().unwrap(),
            ),
            Strand::Reverse => (
                assemblies.iter().map(|a| a.query_end).min().unwrap(),
                assemblies.iter().map(|a| a.query_start).max().unwrap(),
            ),

            Strand::Unset => panic!(),
        };

        Self {
            target_start,
            target_end,
            consensus_start,
            consensus_end,
            consensus_ali_strings,
            target_ali_strings,
            consensus_assembly_strings,
            target_assembly_strings,
            links,
            prev,
            next,
            suffix,
        }
    }
}
