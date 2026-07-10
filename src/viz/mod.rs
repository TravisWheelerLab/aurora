mod bed;
mod block;
mod data;
pub mod debug;
mod stats;

use bed::*;
use block::*;
use data::*;
use stats::*;

use std::{
    collections::{BTreeMap, HashMap},
    fs::{self, File},
    io::{self, BufRead, BufWriter, Write},
    num::ParseIntError,
    path::{Path, PathBuf},
};

use zstd::encode_all;

use anyhow::{anyhow, Context};

use crate::{
    alignment::{Alignment, AlignmentData},
    alphabet::{ALIGNMENT_ALPHABET_UTF8, GAP_EXTEND_DIGITAL, GAP_OPEN_DIGITAL, SPACE_UTF8},
    annotation::AmbiguousAnnotation,
    assembly::SegmentAssemblyGraph,
    chunks::ProximityGroup,
    history_tracing::RefinedTraceSegment,
    matrix::Matrix,
    segments::SegmentedMatrix,
    util::VecMap,
};
use base64::prelude::*;
use itertools::Itertools;

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

pub struct SodaVizWriter {
    viz_path: PathBuf,
    viz_ref_bed_path: Option<PathBuf>,
    bed_index: Option<Vec<(u64, usize)>>,
    constraints: Vec<VizConstraint>,
    hide_simple_regions: bool,
}

/// Check if there are any overlapping alignments in a proximity group...
fn is_simple_proximity_group(group: &ProximityGroup) -> bool {
    // Start -> End
    let mut alignment_starts: BTreeMap<usize, usize> = BTreeMap::new();

    fn check_and_add(tree: &mut BTreeMap<usize, usize>, al_start: usize, al_end: usize) -> bool {
        let before = tree.range(..=al_start).next_back();
        let after = tree.range(al_start + 1..).next();

        if let Some((_, &b_end)) = before {
            if al_start <= b_end {
                return false;
            }
        }

        if let Some((&af_start, _)) = after {
            if al_end >= af_start {
                return false;
            }
        }

        tree.insert(al_start, al_end);
        true
    }

    for al in group.alignments.iter() {
        if !check_and_add(&mut alignment_starts, al.target_start, al.target_end) {
            return false;
        }
    }

    for tr in group.tandem_repeats.iter() {
        if !check_and_add(&mut alignment_starts, tr.target_start, tr.target_end) {
            return false;
        }
    }

    true
}

impl SodaVizWriter {
    pub const SODA_JS: &str = include_str!("../../fixtures/soda/soda.js");
    pub const ICON_SVG: &'static str = include_str!("../../fixtures/soda/icon-opt.svg");
    const INDEX_TEMPLATE: &'static str = include_str!("../../fixtures/soda/index.html");
    const FZSTD_JS: &'static str = include_str!("../../fixtures/soda/fzstd.js");
    const HTML_TEMPLATE: &'static str = include_str!("../../fixtures/soda/annotations.html");
    const JS: &'static str = include_str!("../../fixtures/soda/annotations.js");

    pub fn new(
        proximity_groups: &[ProximityGroup],
        target_name_map: &VecMap<String>,
        viz_path: &impl AsRef<Path>,
        viz_ref_bed_path: Option<&impl AsRef<Path>>,
        constraints: &[VizConstraint],
        hide_simple_regions: bool,
    ) -> anyhow::Result<Self> {
        let viz_path_buf = viz_path.as_ref().to_path_buf();
        let viz_ref_bed_buf = viz_ref_bed_path
            .map(|v| v.as_ref().to_path_buf().canonicalize())
            .transpose()?;

        if let Result::Ok(metadata) = fs::metadata(&viz_path_buf) {
            if metadata.is_dir() {
                // TODO: real error
                return Err(anyhow!(
                    "directory: '{}' already exists",
                    viz_path_buf.to_str().unwrap_or("?")
                ));
            }
        }

        let bed_index = match &viz_ref_bed_buf {
            Some(path) => {
                let mut target_groups: HashMap<&String, (usize, usize)> = HashMap::new();

                for (i, grp) in proximity_groups.iter().enumerate() {
                    let entry = target_groups
                        .entry(target_name_map.get(grp.target_id))
                        .or_insert((i, i + 1));
                    entry.0 = entry.0.min(i);
                    entry.1 = entry.1.max(i + 1);
                }

                let file = File::open(path).context(format!(
                    "failed to open viz reference bed file: '{}'",
                    path.to_str().unwrap_or("?")
                ))?;
                let reader = io::BufReader::new(file);

                let mut chrom_list = vec![String::from("sentinel")];
                let mut prev_start = 0usize;

                let mut line_start = 0u64;

                let mut group_end: usize = 0;
                let mut group_offset: usize = 0;

                // File offset, first line (inclusive), last line (exclusive).
                let mut index: Vec<Option<(u64, usize, usize)>> =
                    vec![None; proximity_groups.len()];
                reader
                    .lines()
                    .enumerate()
                    .try_for_each(|(line_num, res_line)| {
                        let line = res_line?;
                        let trimmed_line = line.trim();

                        if trimmed_line.is_empty() {
                            line_start += line.len() as u64 + 1;
                            return Ok(());
                        }

                        let line_num_info = || format!("failed to read line {}", line_num + 1);

                        let tokens: Vec<&str> = trimmed_line.split_whitespace().collect();

                        if tokens.len() < 3 {
                            return Result::Err(anyhow!("line doesn't have at least 3 columns!"))
                                .with_context(line_num_info);
                        }

                        let chrom = tokens[0].to_string();
                        let start = tokens[1].parse::<usize>().with_context(line_num_info)?;
                        let end = tokens[2].parse::<usize>().with_context(line_num_info)?;

                        let last_chrom = chrom_list.last().context("chrom list is empty")?;

                        if chrom == *last_chrom {
                            if prev_start > start {
                                return Result::Err(anyhow!("bed file is unsorted"));
                            }
                        } else if !chrom_list.contains(&chrom) {
                            chrom_list.push(chrom.clone());
                            let range = target_groups.get(&chrom).unwrap_or(&(0, 0));
                            group_offset = range.0;
                            group_end = range.1;
                        } else {
                            return Result::Err(anyhow!("bed file is unsorted"));
                        }

                        while group_offset < group_end
                            && (start > proximity_groups[group_offset].target_end)
                        {
                            group_offset += 1;
                        }

                        if group_offset < group_end
                            && end >= proximity_groups[group_offset].target_start
                        {
                            let entry =
                                index[group_offset].get_or_insert((line_start, line_num, line_num));
                            entry.2 = line_num + 1;
                        }

                        prev_start = start;
                        line_start += line.len() as u64 + 1; // Include the \n

                        Ok(())
                    })
                    .context(format!(
                        "failed to parse bed file: '{}'",
                        path.to_str().unwrap_or("?")
                    ))?;
                Some(index)
            }
            None => None,
        };

        fs::create_dir_all(&viz_path)?;

        Ok(Self {
            viz_path: viz_path_buf,
            viz_ref_bed_path: viz_ref_bed_buf,
            bed_index: bed_index.map(|v| {
                v.iter()
                    .map(|v| match v {
                        Some((file_start, first_line, last_line)) => {
                            (*file_start, last_line - first_line)
                        }
                        None => (0, 0),
                    })
                    .collect_vec()
            }),
            constraints: constraints.into(),
            hide_simple_regions,
        })
    }

    fn write_index_file(
        writer: &mut impl Write,
        regions: &[ProximityGroup],
        target_name_map: &VecMap<String>,
        viz_constraints: &[VizConstraint],
        exclude_simple_viz: bool,
    ) -> std::io::Result<()> {
        let mut index_links = String::new();

        viz_constraints
            .iter()
            .enumerate()
            .for_each(|(idx, c)| {
                index_links.push_str(&format!(
                    "<div class=\"region\" data-target=\"{name}\" data-start=\"{start}\" data-end=\"{end}\"><a href=\"{name}-{start}-{end}.html\"><h3>slice {idx} | {name} {start}:{end}</h3></a></div>\n",
                    name = c.target_name,
                    start = c.target_start,
                    end = c.target_end,
                    idx = idx
                ));
            });

        regions.iter().enumerate().for_each(|(idx, group)| {
            if exclude_simple_viz && is_simple_proximity_group(group) {
                return;
            }

            index_links.push_str(&format!(
                "<div class=\"region\" data-target=\"{name}\" data-start=\"{start}\" data-end=\"{end}\"><a href=\"{idx}/index.html\"><h3>region {idx} | {name} {start}:{end}</h3></a></div>\n",
                name = target_name_map.get(*&group.target_id),
                start = group.target_start,
                end = group.target_end,
                idx = idx,
            ));
        });

        writeln!(
            writer,
            "{}",
            Self::INDEX_TEMPLATE.replace("INDEX_LINKS_TARGET", &index_links)
        )
    }

    pub fn new_region(
        &self,
        proximity_group: &ProximityGroup,
        alignment_data: &AlignmentData,
        region_idx: usize,
    ) -> Option<RegionAdjudicationSodaWriter> {
        if self.hide_simple_regions && is_simple_proximity_group(proximity_group) {
            return None;
        }

        Some(RegionAdjudicationSodaWriter::new(
            proximity_group,
            alignment_data,
            &self.viz_path,
            region_idx,
            &self.constraints,
            self.viz_ref_bed_path.as_ref(),
            self.bed_index.as_ref().map(|b| b[region_idx]),
        ))
    }

    pub fn finalize(
        &self,
        proximity_groups: &[ProximityGroup],
        annotations: &[(usize, Vec<AmbiguousAnnotation>)],
        target_name_map: &VecMap<String>,
        query_lengths: &HashMap<usize, usize>,
    ) -> io::Result<()> {
        let mut target_stats_writer =
            BufWriter::new(File::create(self.viz_path.join("target_stats.html"))?);
        write_target_statistics(
            &mut target_stats_writer,
            &annotations,
            self.viz_ref_bed_path.as_ref(),
            self.bed_index.as_ref(),
        )?;
        let mut family_stats_writer =
            BufWriter::new(File::create(self.viz_path.join("family_stats.html"))?);
        write_family_statistics(&mut family_stats_writer, &annotations, &query_lengths)?;
        let mut inv_stats_writer =
            BufWriter::new(File::create(self.viz_path.join("inversion_stats.html"))?);
        write_inversion_statistics(&mut inv_stats_writer, &annotations)?;
        let mut icon_file = File::create(self.viz_path.join("icon.svg"))?;
        icon_file.write_all(Self::ICON_SVG.as_bytes())?;

        let mut index_file = BufWriter::new(File::create(self.viz_path.join("index.html"))?);
        Self::write_index_file(
            &mut index_file,
            proximity_groups,
            target_name_map,
            &self.constraints,
            self.hide_simple_regions,
        )?;

        let mut js_file = File::create(self.viz_path.join("annotations.js"))?;
        js_file.write_all(
            Self::JS
                .replace("HTML_TARGET", Self::HTML_TEMPLATE)
                .replace("FZSTD_TARGET", Self::FZSTD_JS)
                .replace("SODA_TARGET", Self::SODA_JS)
                .as_bytes(),
        )?;

        Ok(())
    }
}

pub struct RegionAdjudicationSodaWriter {
    viz_path: PathBuf,
    region_idx: usize,
    has_dumped_confidences: bool,
    finished: bool,
    constraints: Vec<VizConstraint>,
    viz_bed_path: Option<PathBuf>,
    viz_bed_offset: Option<(u64, usize)>,
}

fn to_safe_compressed_string(data: &str) -> io::Result<String> {
    let bytes = encode_all(data.as_bytes(), 1)?;
    Ok(BASE64_STANDARD.encode(bytes))
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
}

impl RegionAdjudicationSodaWriter {
    const DATA_TEMPLATE: &'static str = include_str!("../../fixtures/soda/annotation_data.html");

    fn new(
        proximity_group: &ProximityGroup,
        alignment_data: &AlignmentData,
        viz_path: &impl AsRef<Path>,
        region_idx: usize,
        constraints: &[VizConstraint],
        viz_bed_path: Option<&impl AsRef<Path>>,
        viz_bed_offset: Option<(u64, usize)>,
    ) -> Self {
        Self {
            viz_path: viz_path.as_ref().to_path_buf(),
            region_idx,
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
            viz_bed_path: viz_bed_path.map(|v| v.as_ref().to_path_buf()),
            viz_bed_offset: viz_bed_offset,
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
        let html_start = Self::DATA_TEMPLATE
            .split_once("FILE_SPLIT_POINT")
            .ok_or(io::Error::other(
                "HTML template is broken, missing split point.",
            ))?
            .0;

        // Make the directory for the region if it does not exist...
        let viz_dir = self.viz_path.join(format!("{}", self.region_idx));
        fs::create_dir_all(&viz_dir)?;

        self.write_confidences_single(&viz_dir.join("index.html"), html_start, confidence_matrix)?;

        for constraint in self.constraints.iter() {
            self.write_confidences_single(
                &self.viz_path.join(format!(
                    "{}-{}-{}.html",
                    constraint.target_name, constraint.target_start, constraint.target_end
                )),
                html_start,
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
        confidence_matrix: Option<&Matrix<f64>>,
    ) -> io::Result<()> {
        let html_start = html_start
            .replace("REGION_INDEX", &self.region_idx.to_string())
            .replace(
                "CONFIDENCES_TARGET",
                &to_safe_compressed_string(&self.confidence_json(confidence_matrix)?)?,
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
            self.write_confidences_internal(None)?;
        }

        if self.finished {
            return Err(io::Error::other("Already fully written the visual!"));
        }

        let mut soda_data = AdjudicationSodaData::new(FullAdjudicationSodaDataArgs {
            group: args.group,
            alignment_confidences: args.alignment_confidences,
            active_columns: args.active_columns,
            alignment_data: args.alignment_data,
            annotations: args.annotations,
            target_seq: args.target_seq,
            trace: args.trace,
            segments: args.segments,
            history_counts: args.history_counts,
            links: args.links,
            viz_bed_path: self.viz_bed_path.as_ref(),
            viz_bed_offset: self.viz_bed_offset,
        });
        let html_end = Self::DATA_TEMPLATE
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
            "../annotations.js",
            &mut soda_data,
            None,
        )?;

        for constraint in self.constraints.iter() {
            self.write_single(
                &self.viz_path.join(format!(
                    "{}-{}-{}.html",
                    constraint.target_name, constraint.target_start, constraint.target_end
                )),
                html_end,
                "annotations.js",
                &mut soda_data,
                Some(constraint),
            )?;
        }

        self.finished = true;
        Ok(())
    }

    fn write_single(
        &self,
        path: &PathBuf,
        template: &str,
        js_path: &str,
        args: &mut AdjudicationSodaData,
        constraint: Option<&VizConstraint>,
    ) -> io::Result<()> {
        args.constrain(constraint);

        let result = template
            .replace(
                "DATA_TARGET",
                &to_safe_compressed_string(&serde_json::to_string(&args.to_json()?)?)?,
            )
            .replace("JS_PATH_TARGET", js_path);
        let mut file = File::options().create(false).append(true).open(path)?;

        file.write_all(result.as_bytes())?;

        Ok(())
    }
}
