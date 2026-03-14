mod bed;
mod block;
pub mod debug;
pub mod stats;

use bed::*;
use block::*;

use std::{
    collections::HashMap,
    fs::{self, File},
    io::{self, BufRead, BufReader, Write},
    num::ParseIntError,
    path::{Path, PathBuf},
};

use crate::{
    alignment::{Alignment, AlignmentData},
    alphabet::{
        NucleotideByteUtils, ALIGNMENT_ALPHABET_UTF8, GAP_EXTEND_DIGITAL, GAP_OPEN_DIGITAL,
        SPACE_UTF8,
    },
    annotation::AmbiguousAnnotation,
    assembly::SegmentAssemblyGraph,
    chunks::ProximityGroup,
    history_tracing::{AnnotatedRange, RefinedTraceSegment},
    matrix::Matrix,
    segments::{BlockType, SegmentedMatrix},
    VisualizationArgs,
};
use base64::prelude::*;
use itertools::Itertools;

const SODA_JS: &str = include_str!("../../fixtures/soda/soda.js");
pub const ICON_SVG: &str = include_str!("../../fixtures/soda/icon-opt.svg");
const INDEX_TEMPLATE: &str = include_str!("../../fixtures/soda/index.html");

pub fn write_index_file(
    writer: &mut impl Write,
    alignment_data: &AlignmentData,
    proximity_groups: &[ProximityGroup],
    viz_constraints: &[VizConstraint],
) -> std::io::Result<()> {
    let mut index_links = String::new();

    viz_constraints
        .iter()
        .enumerate()
        .for_each(|(idx, c)| {
            index_links.push_str(&format!(
                "<div class=\"region\" data-target=\"{name}\" data-start=\"{start}\" data-end=\"{end}\"><a href=\"{name}-{start}-{end}.html\">slice {idx} | {name} {start}:{end}</a></div><br>\n",
                name = c.target_name,
                start = c.target_start,
                end = c.target_end,
                idx = idx
            ));
        });

    proximity_groups.iter().enumerate().for_each(|(idx, g)| {
        index_links.push_str(&format!(
            "<div class=\"region\" data-target=\"{name}\" data-start=\"{start}\" data-end=\"{end}\"><a href=\"{idx}/index.html\"><h3>region {idx} | {name} {start}:{end}</h3></a>\n",
            name = alignment_data.target_name_map.get(g.target_id),
            start = g.target_start,
            end = g.target_end,
            idx = idx,
        ));
    });

    writeln!(
        writer,
        "{}",
        INDEX_TEMPLATE.replace("INDEX_LINKS_TARGET", &index_links)
    )
}

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

pub struct AdjudicationSodaWriter {
    viz_path: PathBuf,
    region_idx: usize,
    has_dumped_confidences: bool,
    finished: bool,
    constraints: Vec<VizConstraint>,
}

impl AdjudicationSodaWriter {
    const TEMPLATE: &'static str = include_str!("../../fixtures/soda/annotations.html");
    const JS: &'static str = include_str!("../../fixtures/soda/annotations.js");

    pub fn new(
        proximity_group: &ProximityGroup,
        alignment_data: &AlignmentData,
        viz_path: &impl AsRef<Path>,
        region_idx: usize,
        constraints: &[VizConstraint],
    ) -> Self {
        Self {
            viz_path: viz_path.as_ref().to_path_buf(),
            region_idx: region_idx,
            has_dumped_confidences: false,
            finished: false,
            constraints: constraints
                .iter()
                .filter(|c| {
                    &c.target_name
                        == alignment_data
                            .target_name_map
                            .get(proximity_group.target_id)
                })
                .filter(|c| {
                    c.target_start < proximity_group.target_end
                        && c.target_end > proximity_group.target_start
                })
                .cloned()
                .collect_vec(),
        }
    }

    pub fn write_confidences(&mut self, confidence_matrix: &Matrix<f64>) -> io::Result<()> {
        self.write_confidences_internal(Some(confidence_matrix))
    }

    fn write_confidences_internal(
        &mut self,
        confidence_matrix: Option<&Matrix<f64>>,
    ) -> io::Result<()> {
        if self.has_dumped_confidences {
            return Err(io::Error::other("Attempted to write confidences twice!"));
        }
        // Extract first part of html template...
        let html_start = Self::TEMPLATE
            .split_once("FILE_SPLIT_POINT")
            .ok_or(io::Error::other(
                "HTML template is broken, missing split point.",
            ))?
            .0;

        // Make the directory for the region if it does not exist...
        let viz_dir = self.viz_path.join(format!("{}", self.region_idx));
        fs::create_dir_all(&viz_dir)?;

        self.write_confidences_single(
            &viz_dir.join("index.html"),
            html_start,
            "../",
            confidence_matrix,
        )?;

        for constraint in self.constraints.iter() {
            self.write_confidences_single(
                &self.viz_path.join(format!(
                    "{}-{}-{}.html",
                    constraint.target_name, constraint.target_start, constraint.target_end
                )),
                html_start,
                "",
                confidence_matrix,
            )?;
        }

        self.has_dumped_confidences = true;
        Ok(())
    }

    fn write_confidences_single(
        &self,
        path: &PathBuf,
        html_start: &str,
        relative_path: &str,
        confidence_matrix: Option<&Matrix<f64>>,
    ) -> io::Result<()> {
        let html_start = html_start
            .replace("REGION_INDEX", &self.region_idx.to_string())
            .replace("RELATIVE_PATH_TARGET", relative_path)
            .replace("SODA_TARGET", SODA_JS)
            .replace("JS_TARGET", Self::JS)
            .replace(
                "CONFIDENCE_TARGET",
                &self.confidence_json(confidence_matrix)?,
            );

        let mut file = File::create(path)?;
        file.write_all(html_start.as_bytes())?;

        Ok(())
    }

    fn confidence_json(&self, confidence_matrix: Option<&Matrix<f64>>) -> io::Result<String> {
        let json = serde_json::json!(self.alignment_confidences(confidence_matrix));
        let string = serde_json::to_string(&json)?;
        Ok(string)
    }

    fn alignment_confidences(
        &self,
        confidence_matrix: Option<&Matrix<f64>>,
    ) -> Option<Vec<String>> {
        if let Some(val) = confidence_matrix {
            let region_start = val.def.target_start;

            return Some(
                val.def
                    .col_range_by_logical_row
                    .iter()
                    .enumerate()
                    .map(|(row, &(start, end))| {
                        let max_arr = (start..=end)
                            .map(|col| {
                                *val.data[col]
                                    .iter()
                                    .max_by(|a, b| a.total_cmp(b))
                                    .unwrap_or(&0.0)
                            })
                            .collect_vec();
                        let seq_arr = (start..=end)
                            .zip(max_arr.iter())
                            .flat_map(|(col, max_score)| {
                                ((val.get(row, col).ln() - max_score.ln()) as f32).to_le_bytes()
                            })
                            .collect_vec();
                        let max_arr_enc = max_arr
                            .iter()
                            .flat_map(|v| (v.ln() as f32).to_le_bytes())
                            .collect_vec();

                        format!(
                            "{},{},{},{}",
                            region_start + start,
                            region_start + end,
                            BASE64_STANDARD.encode(seq_arr),
                            BASE64_STANDARD.encode(max_arr_enc)
                        )
                    })
                    .collect_vec(),
            );
        }

        None
    }

    pub fn write(&mut self, args: AdjudicationSodaDataArgs) -> io::Result<()> {
        if !self.has_dumped_confidences {
            self.write_confidences_internal(None);
        }

        if self.finished {
            return Err(io::Error::other("Already fully written the visual!"));
        }

        let soda_data = AdjudicationSodaData::new(args);
        let html_end = Self::TEMPLATE
            .split_once("FILE_SPLIT_POINT")
            .ok_or(io::Error::other(
                "HTML template is broken, missing split point.",
            ))?
            .1;

        self.write_single(
            &self
                .viz_path
                .join(format!("{}/index.html", self.region_idx)),
            html_end,
            &soda_data,
            None,
        )?;

        for constraint in self.constraints.clone().iter() {
            self.write_single(
                &self.viz_path.join(format!(
                    "{}-{}-{}.html",
                    constraint.target_name, constraint.target_start, constraint.target_end
                )),
                html_end,
                &soda_data,
                Some(constraint),
            )?;
        }

        self.finished = true;
        Ok(())
    }

    fn write_single(
        &mut self,
        path: &PathBuf,
        template: &str,
        args: &AdjudicationSodaData,
        constraint: Option<&VizConstraint>,
    ) -> io::Result<()> {
        args.constrain(constraint);

        let result = template.replace("DATA_TARGET", &serde_json::to_string(&args.to_json())?);
        let mut file = File::options().create(false).append(true).open(path)?;

        file.write_all(result.as_bytes())?;

        Ok(())
    }
}

struct AdjudicationSodaData<'a> {
    group: &'a ProximityGroup<'a>,
    alignment_confidences: &'a [f64],
    active_columns: &'a [(usize, usize)],
    alignment_data: &'a AlignmentData,
    target_seq: &'a [u8],
    annotations: &'a [AmbiguousAnnotation],
    trace: &'a Vec<RefinedTraceSegment>,
    maybe_constraint: Option<VizConstraint>,
    segments: &'a SegmentedMatrix,
    history_counts: &'a [usize],
    links: &'a SegmentAssemblyGraph,
    viz_args: &'a VisualizationArgs,
}

pub struct AdjudicationSodaDataArgs<'a> {
    pub group: &'a ProximityGroup<'a>,
    pub alignment_confidences: &'a [f64],
    pub active_columns: &'a [(usize, usize)],
    pub alignment_data: &'a AlignmentData,
    pub annotations: &'a [AmbiguousAnnotation],
    pub target_seq: &'a [u8],
    pub trace: &'a Vec<RefinedTraceSegment>,
    pub segments: &'a SegmentedMatrix,
    pub history_counts: &'a [usize],
    pub links: &'a SegmentAssemblyGraph,
    pub viz_args: &'a VisualizationArgs,
}

impl<'a> AdjudicationSodaData<'a> {
    pub fn new(args: AdjudicationSodaDataArgs<'a>) -> Self {
        Self {
            group: args.group,
            alignment_confidences: args.alignment_confidences,
            active_columns: args.active_columns,
            alignment_data: args.alignment_data,
            target_seq: args.target_seq,
            annotations: args.annotations,
            trace: args.trace,
            maybe_constraint: None,
            segments: args.segments,
            history_counts: args.history_counts,
            links: args.links,
            viz_args: args.viz_args,
        }
    }

    pub fn constrain(&mut self, constraint: Option<&VizConstraint>) {
        self.maybe_constraint = constraint.cloned();
    }

    fn constraint(&self) -> VizConstraint {
        match &self.maybe_constraint {
            Some(constraint) => constraint.clone(),
            None => VizConstraint {
                target_name: String::default(),
                target_start: 0,
                target_end: usize::MAX,
            },
        }
    }

    pub fn to_json(&self) -> serde_json::Value {
        serde_json::json!({
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
        })
    }

    fn history_segments(&self) -> Vec<String> {
        self.segments
            .iter()
            .zip(self.history_counts.iter())
            .map(|(s, h_s)| format!("{},{},{}", s.start_col, s.end_col, h_s))
            .collect()
    }

    fn history_blocks(&self) -> Vec<String> {
        let mut links_per_block: HashMap<(usize, usize), Vec<(usize, usize, f64)>> = HashMap::new();

        for (&(a, b), &w) in self.links.link_graph.iter() {
            links_per_block
                .entry(a)
                .or_default()
                .push((b.0, b.1, w.weight));
            links_per_block
                .entry(b)
                .or_default()
                .push((a.0, a.1, w.weight));
        }

        self.segments
            .iter()
            .enumerate()
            .flat_map(|(s_idx, s)| {
                s.blocks
                    .iter()
                    .enumerate()
                    .map(|(b_idx, b)| {
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
                                self.group.tandem_repeats
                                    [b.row_idx - self.group.alignments.len() - 1]
                                    .consensus_pattern
                            ),
                        };

                        let links = links_per_block
                            .get(&(s_idx, b.row_idx))
                            .iter()
                            .flat_map(|v| v.iter().map(|&v| format!("{}:{}:{}", v.0, v.1, v.2)))
                            .join(";");

                        format!(
                            "{},{},{},{},{},{},{},{},{},{},{}",
                            s_idx,
                            b_idx,
                            b.row_idx,
                            q_id,
                            b.col_start,
                            b.col_end,
                            b.can_join_up_to,
                            b.avg_confidence,
                            b.alignment_score,
                            name,
                            links,
                        )
                    })
                    .collect_vec()
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
                let a_bound = a.get_target_bounds();
                a_bound.0 <= self.constrained_target_end()
                    && a_bound.1 >= self.constrained_target_start()
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
                        .collect::<Vec<&AmbiguousAnnotation>>(),
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
            &self.viz_args.viz_reference_bed_path,
            self.viz_args.viz_reference_bed_index.get(target_name),
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

    fn trace_string(&self, seg: &AnnotatedRange) -> String {
        format!(
            "{},{},{},{},{:3.2}",
            seg.col_start,
            seg.col_end,
            seg.query_id.unwrap_or(0),
            seg.row_idx,
            seg.avg_confidence
        )
    }

    fn conclusive_trace_strings(&self) -> Vec<String> {
        vec![self
            .trace
            .iter()
            // constraint filter
            .filter(|s| {
                let bound = s.max_bounds();

                bound.0 + self.target_start() <= self.constrained_target_end()
                    && bound.1 + self.target_start() >= self.constrained_target_start()
            })
            .flat_map(|seg| seg.annotated.iter().map(|v| self.trace_string(v)))
            .join("|")]
    }

    fn ambiguous_trace_strings(&self) -> Vec<String> {
        vec!["".to_string()]
    }

    fn resolved_assembly_rows(&self) -> Vec<Vec<usize>> {
        let columns = self
            .active_columns
            .iter()
            .flat_map(|&(start, end)| (start..=end))
            .collect_vec();
        vec![columns]
    }

    fn unresolved_assembly_rows(&self) -> Vec<Vec<usize>> {
        vec![vec![]]
    }

    fn competed_assembly_rows(&self) -> Vec<Vec<usize>> {
        vec![vec![]]
    }

    fn inactive_segment_strings(&self) -> Vec<Vec<String>> {
        let mut inactive_col_ranges: Vec<(usize, usize)> = vec![];
        self.active_columns
            .iter()
            .zip(self.active_columns.iter().skip(1))
            .for_each(|(&a, &b)| {
                if b.0 - 1 != a.1 {
                    inactive_col_ranges.push((a.1 + 1, b.0 - 1));
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
                let bounds = seg.max_bounds();
                (
                    bounds.0 + self.target_start(),
                    bounds.1 + self.target_start(),
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
                        let conf = self.alignment_confidences[row_idx];

                        format!(
                            "{},{},{},{:3.2},{},{},{},{},{},{}",
                            seg_target_start.max(ali.target_start),
                            seg_target_end.min(ali.target_end),
                            row_idx,
                            conf,
                            ali.query_start,
                            ali.query_end,
                            self.alignment_data
                                .query_lengths
                                .get(&ali.query_id)
                                .unwrap(),
                            ali.strand,
                            self.alignment_data.query_name_map.get(ali.query_id),
                            ali.id,
                        )
                    })
            })
            .collect_vec();

        vec![conf_strings]
    }
}
