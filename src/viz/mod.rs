mod bed;
mod block;

use bed::*;
use block::*;

use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    num::ParseIntError,
    path::Path,
};

use itertools::Itertools;
use serde::Serialize;

use crate::{
    alignment::{Alignment, AlignmentData, Strand},
    alphabet::{
        NucleotideByteUtils, ALIGNMENT_ALPHABET_UTF8, GAP_EXTEND_DIGITAL, GAP_OPEN_DIGITAL,
        SPACE_UTF8,
    },
    annotation::Annotation,
    collapse::{Assembly, AssemblyGroup},
    matrix::Matrix,
    split::SplitResults,
    viterbi::TraceSegment,
    AuroraArgs,
};

const SODA_JS: &str = include_str!("../../fixtures/soda/soda.js");

#[derive(Clone, Debug)]
pub struct VizConstraint {
    pub target_name: String,
    pub target_start: usize,
    pub target_end: usize,
}

impl std::str::FromStr for VizConstraint {
    type Err = ParseIntError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let tokens: Vec<&str> = s.split(':').collect();

        Ok(VizConstraint {
            target_name: tokens[0].to_string(),
            target_start: tokens[1].parse()?,
            target_end: tokens[2].parse()?,
        })
    }
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

pub struct AdjudicationSodaData<'a> {
    group: &'a AssemblyGroup<'a>,
    confidence_matrix: &'a Matrix<'a, f64>,
    alignment_data: &'a AlignmentData,
    target_seq: &'a [u8],
    annotations: Vec<Annotation>,
    split_results: Vec<SplitResults>,
    maybe_constraint: Option<&'a VizConstraint>,
    args: &'a AuroraArgs,
}

impl<'a> AdjudicationSodaData<'a> {
    const TEMPLATE: &'static str = include_str!("../../fixtures/soda/annotations.html");
    const JS: &'static str = include_str!("../../fixtures/soda/annotations.js");

    pub fn new(
        group: &'a AssemblyGroup,
        confidence_matrix: &'a Matrix<'a, f64>,
        alignment_data: &'a AlignmentData,
        target_seq: &'a [u8],
        args: &'a AuroraArgs,
    ) -> Self {
        Self {
            group,
            confidence_matrix,
            alignment_data,
            target_seq,
            annotations: vec![],
            split_results: vec![],
            maybe_constraint: None,
            args,
        }
    }

    pub fn constrain(&mut self, constraint: &'a VizConstraint) {
        self.maybe_constraint = Some(constraint);
    }

    fn constraint(&self) -> VizConstraint {
        match self.maybe_constraint {
            Some(constraint) => constraint.clone(),
            None => VizConstraint {
                target_name: String::default(),
                target_start: 0,
                target_end: usize::MAX,
            },
        }
    }

    pub fn write(&self, path: impl AsRef<Path>) {
        let data = serde_json::json!({
            "targetStart": self.constrained_target_start(),
            "targetEnd": self.constrained_target_end(),
            "targetSeq": self.target_seq(),
            "numQueries": self.num_queries(),
            "auroraAnn": self.aurora_ann(),
            "referenceAnn": self.reference_ann(),
            "alignmentStrings": self.alignment_strings(),
            "assemblyStrings": self.assembly_strings(),
            "tandemRepeatStrings": self.tandem_repeat_strings(),
            "conclusiveTraceStrings": self.conclusive_trace_strings(),
            "ambiguousTraceStrings": self.ambiguous_trace_strings(),
            "resolvedAssemblyRows": self.resolved_assembly_rows(),
            "unresolvedAssemblyRows": self.unresolved_assembly_rows(),
            "competedAssemblyRows": self.competed_assembly_rows(),
            "inactiveSegmentStrings": self.inactive_segment_strings(),
            "confidenceSegmentStrings": self.confidence_segment_strings(),
        });

        let viz_html = Self::TEMPLATE
            .replace("SODA_TARGET", SODA_JS)
            .replace(
                "DATA_TARGET",
                &serde_json::to_string(&data).expect("failed to serialize JSON data"),
            )
            .replace("JS_TARGET", Self::JS);

        let mut file = std::fs::File::create(path).expect("failed to create file");

        std::io::Write::write_all(&mut file, viz_html.as_bytes()).expect("failed to write to file");
    }

    pub fn add(&mut self, results: SplitResults) {
        self.split_results.push(results);
    }

    pub fn set_annotations(&mut self, annotations: Vec<Annotation>) {
        self.annotations = annotations;
    }

    fn target_start(&self) -> usize {
        self.group.target_start
    }

    fn target_end(&self) -> usize {
        self.group.target_end
    }

    fn constrained_target_start(&self) -> usize {
        self.group.target_start.max(self.constraint().target_start)
    }

    fn constrained_target_end(&self) -> usize {
        self.group.target_end.min(self.constraint().target_end)
    }

    fn num_queries(&self) -> usize {
        self.group.assemblies.len()
    }

    fn target_seq(&self) -> String {
        let start_idx = self.constrained_target_start() - self.target_start();
        let end_idx = self.constrained_target_end() - self.target_start();
        self.target_seq[start_idx..=end_idx].to_utf8_string()
    }

    fn aurora_ann(&self) -> Vec<BlockGroup> {
        let unique_join_ids: Vec<usize> = self
            .annotations
            .iter()
            // constraint filter
            .filter(|a| {
                a.target_start <= self.constrained_target_end()
                    && a.target_end >= self.constrained_target_start()
            })
            .map(|a| a.join_id)
            .unique()
            .collect();

        unique_join_ids
            .iter()
            .map(|&id| {
                BlockGroup::from_joined_annotations(
                    &mut self
                        .annotations
                        .iter()
                        .filter(|&a| a.join_id == id)
                        .collect::<Vec<&Annotation>>(),
                    &self.alignment_data.query_lengths,
                )
            })
            .collect()
    }

    fn reference_ann(&self) -> Vec<BlockGroup> {
        let mut overlapping_bed = vec![];

        let target_name = self
            .alignment_data
            .target_name_map
            .get(self.group.target_id);

        if let (Some(path), Some(&offset)) = (
            &self.args.visualization_args.viz_reference_bed_path,
            self.args
                .visualization_args
                .viz_reference_bed_index
                .get(target_name),
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
                    if thick_start < self.target_end() && thick_end > self.target_start() {
                        overlapping_bed.push(BedRecord::from_tokens(&tokens));
                    }
                });
        }

        overlapping_bed
            .iter()
            // constraint filter
            .filter(|b| {
                b.chrom_start <= self.constrained_target_end()
                    && b.chrom_end >= self.constrained_target_start()
            })
            .map(BlockGroup::from_bed_record)
            .collect()
    }

    fn alignment_strings(&self) -> Vec<String> {
        self.group
            .assemblies
            .iter()
            .enumerate()
            // constraint filter
            .filter(|(_, a)| {
                a.target_start <= self.constrained_target_end()
                    && a.target_end >= self.constrained_target_start()
            })
            .flat_map(|(idx, assembly)| {
                assembly.alignments.iter().map(move |a| {
                    a.soda_string(idx + 1, self.alignment_data.query_name_map.get(a.query_id))
                })
            })
            .collect()
    }

    fn assembly_strings(&self) -> Vec<String> {
        self.group
            .assemblies
            .iter()
            .enumerate()
            // constraint filter
            .filter(|(_, a)| {
                a.target_start <= self.constrained_target_end()
                    && a.target_end >= self.constrained_target_start()
            })
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
            .collect()
    }

    fn tandem_repeat_strings(&self) -> Vec<String> {
        self.group
            .tandem_repeats
            .iter()
            .enumerate()
            // constraint filter
            .filter(|(_, r)| {
                r.target_start <= self.constrained_target_end()
                    && r.target_end >= self.constrained_target_start()
            })
            .map(|(repeat_idx, repeat)| {
                format!(
                    "{},{},{},{},{}",
                    repeat.target_start,
                    repeat.target_end,
                    repeat.consensus_pattern,
                    repeat.period,
                    repeat_idx + self.group.assemblies.len() + 1,
                )
            })
            .collect()
    }

    fn trace_string(&self, seg: &TraceSegment) -> String {
        let mut conf = 0.0;
        (seg.col_start..=seg.col_end)
            .for_each(|col_idx| conf += self.confidence_matrix.get(seg.row_idx, col_idx));

        conf /= (seg.col_end - seg.col_start + 1) as f64;

        format!(
            "{},{},{},{},{:3.2}",
            seg.col_start, seg.col_end, seg.query_id, seg.row_idx, conf
        )
    }

    fn conclusive_trace_strings(&self) -> Vec<String> {
        self.split_results
            .iter()
            .map(|r| {
                r.trace_conclusive
                    .iter()
                    // constraint filter
                    .filter(|s| {
                        s.col_start + self.target_start() <= self.constrained_target_end()
                            && s.col_end + self.target_start() >= self.constrained_target_start()
                    })
                    .map(|seg| self.trace_string(seg))
                    .join("|")
            })
            .collect()
    }

    fn ambiguous_trace_strings(&self) -> Vec<String> {
        self.split_results
            .iter()
            .map(|r| {
                r.trace_ambiguous
                    .iter()
                    // constraint filter
                    .filter(|s| {
                        s.col_start + self.target_start() <= self.constrained_target_end()
                            && s.col_end + self.target_start() >= self.constrained_target_start()
                    })
                    .map(|seg| self.trace_string(seg))
                    .join("|")
            })
            .collect()
    }

    fn resolved_assembly_rows(&self) -> Vec<Vec<usize>> {
        self.split_results
            .iter()
            .map(|r| r.resolved_assembly_rows.clone())
            .collect()
    }

    fn unresolved_assembly_rows(&self) -> Vec<Vec<usize>> {
        self.split_results
            .iter()
            .map(|r| r.unresolved_assembly_rows.clone())
            .collect()
    }

    fn competed_assembly_rows(&self) -> Vec<Vec<usize>> {
        self.split_results
            .iter()
            .map(|r| r.competed_assembly_rows.clone())
            .collect()
    }

    fn inactive_segment_strings(&self) -> Vec<Vec<String>> {
        self.split_results
            .iter()
            .map(|r| {
                r.inactive_col_ranges
                    .iter()
                    .map(|s| {
                        format!(
                            "{},{}",
                            s.col_start + self.target_start(),
                            s.col_end + self.target_start()
                        )
                    })
                    .collect_vec()
            })
            .collect()
    }

    fn confidence_segment_strings(&self) -> Vec<Vec<String>> {
        self.split_results
            .iter()
            .map(|r| {
                r.trace_ambiguous
                    .iter()
                    .chain(r.trace_conclusive.iter())
                    // map each trace segment to
                    // target start & end coordinates
                    .map(|seg| {
                        (
                            seg.col_start + self.target_start(),
                            seg.col_end + self.target_start(),
                        )
                    })
                    .flat_map(|(seg_target_start, seg_target_end)| {
                        self.group
                            .assemblies
                            .iter()
                            .enumerate()
                            // constraint filter
                            .filter(|(_, a)| {
                                a.target_start <= self.constrained_target_end()
                                    && a.target_end >= self.constrained_target_start()
                            })
                            // get the row idx of the assembly
                            .map(|(i, a)| (i + 1, a))
                            .flat_map(move |(row_idx, assembly)| {
                                assembly
                                    .alignments
                                    .iter()
                                    .filter(move |ali| {
                                        ali.target_start < seg_target_end
                                            && ali.target_end > seg_target_start
                                    })
                                    .map(move |ali| {
                                        (
                                            row_idx,
                                            ali,
                                            seg_target_start.max(ali.target_start)
                                                - self.target_start(),
                                            seg_target_end.min(ali.target_end)
                                                - self.target_start(),
                                        )
                                    })
                                    .map(move |(row_idx, ali, col_start, col_end)| {
                                        let mut conf = 0.0;
                                        (col_start..=col_end).for_each(|col_idx| {
                                            conf += self.confidence_matrix.get(row_idx, col_idx)
                                        });
                                        conf /= (col_end - col_start + 1) as f64;
                                        format!(
                                            "{},{},{},{:3.2},{},{},{},{},{}",
                                            seg_target_start.max(ali.target_start),
                                            seg_target_end.min(ali.target_end),
                                            row_idx,
                                            conf,
                                            self.confidence_matrix
                                                .consensus_position(row_idx, col_start),
                                            self.confidence_matrix
                                                .consensus_position(row_idx, col_end),
                                            self.alignment_data
                                                .query_lengths
                                                .get(&ali.query_id)
                                                .unwrap(),
                                            self.confidence_matrix.strand_of_row(row_idx),
                                            self.alignment_data.query_name_map.get(ali.query_id),
                                        )
                                    })
                            })
                    })
                    .collect_vec()
            })
            .collect()
    }
}

///
///
///
///
#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub struct AssemblySodaData {
    query_id: usize,
    query_name: String,
    target_start: usize,
    target_end: usize,
    consensus_start: usize,
    consensus_end: usize,
    consensus_ali_strings: Vec<String>,
    target_ali_strings: Vec<String>,
    consensus_assembly_strings: Vec<Vec<String>>,
    target_assembly_strings: Vec<Vec<String>>,
    links: Vec<String>,
    suffix: String,
}

impl AssemblySodaData {
    const TEMPLATE: &'static str = include_str!("../../fixtures/soda/assembly.html");
    const JS: &'static str = include_str!("../../fixtures/soda/assembly.js");

    pub fn new(
        assemblies: &[Assembly],
        links: Vec<String>,
        confidence: &HashMap<usize, f64>,
        alignment_data: &AlignmentData,
    ) -> Self {
        let query_id = assemblies[0].query_id;
        let strand = assemblies[0].strand;
        let query_name = alignment_data.query_name_map.get(query_id).clone();

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
            query_id,
            query_name,
            target_start,
            target_end,
            consensus_start,
            consensus_end,
            consensus_ali_strings,
            target_ali_strings,
            consensus_assembly_strings,
            target_assembly_strings,
            links,
            suffix,
        }
    }

    pub fn write(&self, path: impl AsRef<Path>) {
        let viz_html = Self::TEMPLATE
            .replace("SODA_TARGET", SODA_JS)
            .replace(
                "DATA_TARGET",
                &serde_json::to_string(&self).expect("failed to serialize JSON data"),
            )
            .replace("JS_TARGET", Self::JS);

        let mut file = std::fs::File::create(path).expect("failed to create file");

        std::io::Write::write_all(&mut file, viz_html.as_bytes()).expect("failed to write to file");
    }
}
