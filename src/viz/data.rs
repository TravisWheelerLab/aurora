use crate::{
    alignment::AlignmentData,
    alphabet::NucleotideByteUtils,
    annotation::AmbiguousAnnotation,
    assembly::SegmentAssemblyGraph,
    chunks::ProximityGroup,
    history_tracing::{AnnotatedRange, RefinedTraceSegment},
    segments::{BlockType, SegmentedMatrix},
    viz::bed::BedRecord,
    viz::block::BlockGroup,
    viz::VizConstraint,
};
use itertools::Itertools;
use std::io;
use std::io::{BufRead, BufReader};
use std::{collections::HashMap, path::PathBuf};
use std::{
    fs::File,
    io::{Seek, SeekFrom},
};

pub struct FullAdjudicationSodaDataArgs<'a> {
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
    pub viz_bed_path: Option<&'a PathBuf>,
    pub viz_bed_offset: Option<(u64, usize)>,
}

pub struct AdjudicationSodaData<'a> {
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
    viz_bed_path: Option<&'a PathBuf>,
    viz_bed_offset: Option<(u64, usize)>,
}

impl<'a> AdjudicationSodaData<'a> {
    pub fn new(args: FullAdjudicationSodaDataArgs<'a>) -> Self {
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
            viz_bed_path: args.viz_bed_path,
            viz_bed_offset: args.viz_bed_offset,
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

    pub fn to_json(&self) -> io::Result<serde_json::Value> {
        Ok(serde_json::json!({
            "targetStart": self.constrained_target_start(),
            "targetEnd": self.constrained_target_end(),
            "targetSeq": self.target_seq(),
            "numQueries": self.num_queries(),
            "assemblyStrings": self.assembly_strings(),
            "auroraAnn": self.aurora_ann(),
            "referenceAnn": self.reference_ann()?,
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
        }))
    }

    fn history_segments(&self) -> Vec<String> {
        self.segments
            .iter()
            .zip(self.history_counts.iter())
            .map(|(s, h_s)| format!("{},{},{}", s.start_col, s.end_col, h_s))
            .collect()
    }

    fn history_blocks(&self) -> Vec<String> {
        type BlockLocation = (usize, usize); // segment index, dense block index.
        type BlockEdge = (usize, usize, f64); // other segment index, other dense block index, edge weight

        let mut links_per_block: HashMap<BlockLocation, Vec<BlockEdge>> = HashMap::new();

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

    fn reference_ann(&self) -> io::Result<Vec<BlockGroup>> {
        let mut overlapping_bed = vec![];

        let target_name = self
            .alignment_data
            .target_name_map
            .get(self.group.target_id);

        if let (Some(path), Some(offset)) = (self.viz_bed_path, self.viz_bed_offset) {
            let file = File::open(path).expect("failed to open reference bed");
            let mut reader = BufReader::new(file);

            reader.seek(SeekFrom::Start(offset.0))?;

            for raw_line in reader.lines().take(offset.1) {
                let raw_line = raw_line?;
                let line = raw_line.trim();

                if line.is_empty() {
                    continue;
                }

                let tokens: Vec<&str> = line.split_whitespace().collect();

                let target = tokens[0];
                if target != target_name {
                    break;
                }

                let thick_start = tokens[6].parse::<usize>().expect("failed to parse usize");
                let thick_end = tokens[7].parse::<usize>().expect("failed to parse usize");
                if thick_start < self.target_end() && thick_end > self.target_start() {
                    overlapping_bed
                        .push(BedRecord::from_tokens(&tokens).map_err(|e| io::Error::other(e))?);
                }
            }
        }

        Ok(overlapping_bed
            .iter()
            // constraint filter
            .filter(|b| {
                b.chrom_start <= self.constrained_target_end()
                    && b.chrom_end >= self.constrained_target_start()
            })
            .map(BlockGroup::from_bed_record)
            .collect())
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
            .flat_map(|&(start, end)| start..=end)
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
