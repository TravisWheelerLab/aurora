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
    AuroraArgs, alignment::{Alignment, AlignmentData, Strand}, alphabet::{
        ALIGNMENT_ALPHABET_UTF8, GAP_EXTEND_DIGITAL, GAP_OPEN_DIGITAL, NucleotideByteUtils, SPACE_UTF8
    }, annotation::Annotation, assembly::AssemblyGraph, chunks::ProximityGroup, matrix::Matrix, segments::{BlockType, SegmentedMatrix}, viterbi::{RefinedTraceSegment, TraceSegment}
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
        let mut red_bytes: Vec<Vec<u8>> = vec![];
        let mut red_starts: Vec<usize> = vec![];

        let mut t_pos = self.target_start;
        let mut in_target_gap = false;
        let mut current_red_bytes: Vec<u8> = vec![];
        let mut target_gap_start = 0;

        self.target_seq
            .iter()
            .zip(&self.query_seq)
            // .enumerate()
            // .map(move |(ali_idx, c)| (ali_idx + target_pos, c))
            // .for_each(|(t_idx, (&t, &q))| match t {
            .for_each(|(&t, &q)| match t {
                GAP_OPEN_DIGITAL | GAP_EXTEND_DIGITAL => {
                    if in_target_gap {
                        current_red_bytes.push(ALIGNMENT_ALPHABET_UTF8[q as usize]);
                    } else {
                        in_target_gap = true;
                        target_gap_start = t_pos;
                        current_red_bytes = vec![ALIGNMENT_ALPHABET_UTF8[q as usize]]
                    }
                }
                _ => {
                    t_pos += 1;
                    if in_target_gap {
                        red_bytes.push(current_red_bytes.clone());
                        red_starts.push(target_gap_start);
                    }
                    in_target_gap = false;
                    if t == q {
                        green_bytes.push(ALIGNMENT_ALPHABET_UTF8[q as usize]);
                        orange_bytes.push(SPACE_UTF8);
                    } else {
                        green_bytes.push(SPACE_UTF8);
                        orange_bytes.push(ALIGNMENT_ALPHABET_UTF8[q as usize]);
                    }
                }
            });

        let green_string = String::from_utf8(green_bytes).unwrap();
        let orange_string = String::from_utf8(orange_bytes).unwrap();
        let mut red_string = red_bytes
            .into_iter()
            .zip(red_starts)
            .fold("".to_string(), |acc, (b, s)| {
                format!("{acc}{}:{s}|", String::from_utf8(b).unwrap())
            });

        red_string.pop();

        debug_assert_eq!(self.target_end - self.target_start + 1, green_string.len());

        format!(
            "{},{},{},{},{},{},{},{},{}",
            green_string,
            orange_string,
            red_string,
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
    group: &'a ProximityGroup<'a>,
    confidence_matrix: &'a Matrix<'a, f64>,
    alignment_data: &'a AlignmentData,
    target_seq: &'a [u8],
    annotations: Vec<Annotation>,
    trace: &'a Vec<RefinedTraceSegment>,
    maybe_constraint: Option<&'a VizConstraint>,
    segments: &'a SegmentedMatrix,
    history_counts: &'a [usize],
    links: &'a AssemblyGraph,
    args: &'a AuroraArgs,
}

impl<'a> AdjudicationSodaData<'a> {
    const TEMPLATE: &'static str = include_str!("../../fixtures/soda/annotations.html");
    const JS: &'static str = include_str!("../../fixtures/soda/annotations.js");

    pub fn new(
        group: &'a ProximityGroup,
        confidence_matrix: &'a Matrix<'a, f64>,
        alignment_data: &'a AlignmentData,
        target_seq: &'a [u8],
        trace: &'a Vec<RefinedTraceSegment>,
        segments: &'a SegmentedMatrix,
        history_counts: &'a [usize],
        links: &'a AssemblyGraph,
        args: &'a AuroraArgs,
    ) -> Self {
        Self {
            group,
            confidence_matrix,
            alignment_data,
            target_seq,
            annotations: vec![],
            trace,
            maybe_constraint: None,
            segments,
            history_counts,
            links,
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
            "assemblyStrings": self.assembly_strings(),
            "auroraAnn": self.aurora_ann(),
            "referenceAnn": self.reference_ann(),
            "alignmentStrings": self.alignment_strings(),
            "tandemRepeatStrings": self.tandem_repeat_strings(),
            "conclusiveTraceStrings": self.conclusive_trace_strings(),
            "ambiguousTraceStrings": self.ambiguous_trace_strings(),
            "resolvedAssemblyRows": self.resolved_assembly_rows(),
            "unresolvedAssemblyRows": self.unresolved_assembly_rows(),
            "competedAssemblyRows": self.competed_assembly_rows(),
            "inactiveSegmentStrings": self.inactive_segment_strings(),
            "confidenceSegmentStrings": self.confidence_segment_strings(),
            "historySegments": self.history_segments(),
            "historyBlocks": self.history_blocks(),
            "blockLinks": self.block_links()
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

    fn block_links(&self) -> Vec<Vec<String>> {
        self.links.link_graph
            .iter()
            .map(|links| {
                links.iter().map(|edge| {
                    format!("{},{}", edge.edge_to, edge.weight)
                }).collect()
            }).collect()
    }

    fn history_segments(&self) -> Vec<String> {
        self.segments
            .iter()
            .zip(self.history_counts.iter())
            .map(|(s, h_s)| format!("{},{},{}", s.start_col, s.end_col, h_s))
            .collect()
    }

    fn history_blocks(&self) -> Vec<String> {
        self.segments
            .iter()
            .enumerate()
            .flat_map(|(s_idx, s)| {
                s.blocks.iter().enumerate().map(move |(b_idx, b)| {
                    let q_id = match b.query_id {
                        Some(v) => v.to_string(),
                        _ => (-1).to_string(),
                    };
                    let name = match b.block_type {
                        BlockType::Skip => "Skip".to_string(),
                        BlockType::Alignment => self
                            .alignment_data
                            .query_name_map
                            .get(b.query_id.unwrap())
                            .to_string(),
                        BlockType::TandemRepeat => format!(
                            "repeat#{}",
                            self.group.tandem_repeats[b.row_idx].consensus_pattern
                        ),
                    };

                    format!(
                        "{},{},{},{},{},{},{},{},{}",
                        s_idx,
                        b_idx,
                        b.row_idx,
                        q_id,
                        b.target_start,
                        b.target_end,
                        b.can_join_up_to,
                        b.confidence,
                        name
                    )
                })
            })
            .collect()
    }

    fn assembly_strings(&self) -> Vec<String> {
        self.group
            .alignments
            .iter()
            .enumerate()
            .map(|(idx, ali)| {
                format!(
                    "{},{},{},{},{}",
                    ali.target_start,
                    ali.target_end,
                    ali.query_id,
                    1,
                    idx + 1
                )
            })
            .collect()
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
        self.group.alignments.len()
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
            .alignments
            .iter()
            .enumerate()
            // constraint filter
            .filter(|(_, a)| {
                a.target_start <= self.constrained_target_end()
                    && a.target_end >= self.constrained_target_start()
            })
            .map(|(idx, alignment)| {
                alignment.soda_string(
                    idx + 1,
                    self.alignment_data.query_name_map.get(alignment.query_id),
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
                    repeat_idx + self.group.alignments.len() + 1,
                )
            })
            .collect()
    }

    fn trace_string(&self, seg: &RefinedTraceSegment) -> String {
        let mut conf = 0.0;
        (seg.col_start..=seg.col_end)
            .for_each(|col_idx| conf += self.confidence_matrix.get(seg.row_idx, col_idx));

        conf /= (seg.col_end - seg.col_start + 1) as f64;

        format!(
            "{},{},{},{},{:3.2}",
            seg.col_start,
            seg.col_end,
            seg.query_id.unwrap_or(0),
            seg.row_idx,
            conf
        )
    }

    fn conclusive_trace_strings(&self) -> Vec<String> {
        vec![self
            .trace
            .iter()
            // constraint filter
            .filter(|s| {
                s.col_start + self.target_start() <= self.constrained_target_end()
                    && s.col_end + self.target_start() >= self.constrained_target_start()
            })
            .map(|seg| self.trace_string(seg))
            .join("|")]
    }

    fn ambiguous_trace_strings(&self) -> Vec<String> {
        vec!["".to_string()]
    }

    fn resolved_assembly_rows(&self) -> Vec<Vec<usize>> {
        vec![self.confidence_matrix.initial_active_cols()]
    }

    fn unresolved_assembly_rows(&self) -> Vec<Vec<usize>> {
        vec![vec![]]
    }

    fn competed_assembly_rows(&self) -> Vec<Vec<usize>> {
        vec![vec![]]
    }

    fn inactive_segment_strings(&self) -> Vec<Vec<String>> {
        let active_cols = self.confidence_matrix.initial_active_cols();
        let mut inactive_col_ranges: Vec<(usize, usize)> = vec![];
        active_cols
            .iter()
            .zip(active_cols.iter().skip(1))
            .for_each(|(&a, &b)| {
                if b - 1 != a {
                    inactive_col_ranges.push((a + 1, b - 1));
                }
            });

        vec![inactive_col_ranges
            .iter()
            .map(|(start, end)| {
                format!(
                    "{},{}",
                    start + self.target_start(),
                    end + self.target_start()
                )
            })
            .collect_vec()]
    }

    fn confidence_segment_strings(&self) -> Vec<Vec<String>> {
        let conf_strings = self
            .trace
            .iter()
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
                    .alignments
                    .iter()
                    .enumerate()
                    // constraint filter
                    .filter(move |(_, a)| {
                        a.target_start <= self.constrained_target_end()
                            && a.target_end >= self.constrained_target_start()
                            && a.target_start < seg_target_end
                            && a.target_end > seg_target_start
                    })
                    // get the row idx of the assembly
                    .map(|(i, a)| (i + 1, a))
                    .map(move |(row_idx, ali)| {
                        let col_start =
                            seg_target_start.max(ali.target_start) - self.target_start();
                        let col_end = seg_target_end.min(ali.target_end) - self.target_start();

                        let mut conf = 0.0;
                        (col_start..=col_end).for_each(|col_idx| {
                            conf += self.confidence_matrix.get(row_idx, col_idx)
                        });
                        conf /= (col_end - col_start + 1) as f64;
                        format!(
                            "{},{},{},{:3.2},{},{},{},{},{},{}",
                            seg_target_start.max(ali.target_start),
                            seg_target_end.min(ali.target_end),
                            row_idx,
                            conf,
                            self.confidence_matrix
                                .consensus_position(row_idx, col_start),
                            self.confidence_matrix.consensus_position(row_idx, col_end),
                            self.alignment_data
                                .query_lengths
                                .get(&ali.query_id)
                                .unwrap(),
                            self.confidence_matrix.strand_of_row(row_idx),
                            self.alignment_data.query_name_map.get(ali.query_id),
                            ali.id,
                        )
                    })
            })
            .collect_vec();

        vec![conf_strings]
    }
}

