use crate::{
    alignment::Strand,
    annotation::AmbiguousAnnotation,
    assembly::{
        block_target_distance, relative_consensus_distance, ConsensusDistanceNormalization,
    },
    segments::{Block, Unordered},
    viz::{bed::BedRecord, SodaVizWriter},
};
use core::str;
use itertools::Itertools;
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufRead, BufReader, Read, Seek, SeekFrom, Write},
    path::Path,
};

const TABLE_HTML: &str = include_str!("../../fixtures/soda/table.html");

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Default)]
struct FamilyInfo {
    pub occurrences: usize,
    pub coverage: usize,
    pub kimura80_values: Unordered<Vec<f64>>,
    pub join_consensus_dist: Unordered<Vec<f64>>,
    pub nojoin_consensus_dist: Unordered<Vec<f64>>,
    pub join_target_dist: Unordered<Vec<f64>>,
    pub nojoin_target_dist: Unordered<Vec<f64>>,
    pub join_kimura_dist: Unordered<Vec<f64>>,
    pub nojoin_kimura_dist: Unordered<Vec<f64>>,
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
struct InversionInfo {
    pub inversions: usize,
    pub normal_joins: usize,
}

#[derive(Debug, Default, PartialEq, Eq, PartialOrd, Ord)]
struct TargetInfo {
    pub te_count: usize,
    pub te_coverage: usize,
    pub joins: usize,
    pub sr_count: usize,
    pub sr_coverage: usize,
    pub total_coverage: usize,
}

fn write_tsv<const N: usize, A: std::fmt::Display, B: std::fmt::Display>(
    writer: &mut impl Write,
    header: &[A; N],
    data: &[[B; N]],
) -> std::io::Result<()> {
    writeln!(writer, "{}", header.iter().join("\t"))?;

    for row in data.iter() {
        writeln!(writer, "{}", row.iter().join("\t"))?;
    }

    Ok(())
}

fn write_statistics_table_page<const N: usize, A: std::fmt::Display, B: std::fmt::Display>(
    writer: &mut impl Write,
    title: &str,
    header: &[A; N],
    data: &[[B; N]],
) -> std::io::Result<()> {
    let mut tmp_writer = Vec::<u8>::new();

    write_tsv(&mut tmp_writer, header, data)?;

    let table_page = TABLE_HTML
        .replace("PAGE_TITLE", title)
        .replace("SODA_TARGET", SodaVizWriter::SODA_JS)
        .replace(
            "TSV_TARGET",
            str::from_utf8(&tmp_writer).expect("UTF8 decoding failed!"),
        );

    write!(writer, "{}", table_page)?;

    Ok(())
}

pub fn write_target_statistics(
    stats_writer: &mut impl Write,
    region_annotations: &[(usize, Vec<AmbiguousAnnotation>)],
    bed_file: Option<impl AsRef<Path>>,
    bed_index: Option<&Vec<(u64, usize)>>,
) -> std::io::Result<()> {
    let mut aurora_target_stats: HashMap<&String, TargetInfo> = HashMap::new();

    fn summary<'a>(data: impl Iterator<Item = &'a TargetInfo>) -> TargetInfo {
        let mut summary_data = TargetInfo::default();

        for val in data {
            summary_data.joins += val.joins;
            summary_data.sr_count += val.sr_count;
            summary_data.sr_coverage += val.sr_coverage;
            summary_data.te_count += val.te_count;
            summary_data.te_coverage += val.te_coverage;
            summary_data.total_coverage += val.total_coverage;
        }

        summary_data
    }

    fn to_table_row(name: &str, data: &TargetInfo) -> [String; 7] {
        [
            name.to_string(),
            data.te_count.to_string(),
            data.te_coverage.to_string(),
            data.sr_count.to_string(),
            data.sr_coverage.to_string(),
            data.joins.to_string(),
            data.total_coverage.to_string(),
        ]
    }

    for (_region, annotations) in region_annotations.iter() {
        let mut join_ids_seen: HashSet<usize> = HashSet::new();

        for amb_annot in annotations.iter() {
            let is_a_join = !join_ids_seen.insert(amb_annot.join_id);
            let t_start = amb_annot
                .annotations
                .iter()
                .map(|v| v.target_start)
                .min()
                .unwrap_or(0);
            let t_end = amb_annot
                .annotations
                .iter()
                .map(|v| v.target_end)
                .max()
                .unwrap_or(0);
            let coverage = t_end.saturating_sub(t_start);
            let is_sr = amb_annot.annotations.iter().any(|v| v.query_id == 0);

            let entry = aurora_target_stats
                .entry(&amb_annot.target_name)
                .or_insert_with(TargetInfo::default);

            if is_sr {
                entry.sr_count += !is_a_join as usize;
                entry.sr_coverage += coverage;
            } else {
                entry.te_count += !is_a_join as usize;
                entry.te_coverage += coverage;
            }
            entry.joins += is_a_join as usize;
            entry.total_coverage += coverage;
        }
    }

    let mut table_stats = Vec::new();

    table_stats.extend(
        aurora_target_stats
            .iter()
            .map(|(n, v)| to_table_row(&format!("{} Aurora", n), v)),
    );

    table_stats.push(to_table_row(
        "Totals Aurora",
        &summary(aurora_target_stats.values()),
    ));

    if let (Some(bed_path), Some(bed_index)) = (bed_file, bed_index) {
        let mut reader = BufReader::new(File::open(bed_path)?);
        let mut rm_results: HashMap<String, TargetInfo> = HashMap::new();

        for &(offset, total_lines) in bed_index {
            reader.seek(SeekFrom::Start(offset))?;

            for raw_line in reader.by_ref().lines().take(total_lines) {
                let line = raw_line?;

                if line.trim().is_empty() {
                    continue;
                }

                let bed_record = BedRecord::from_str(line.trim())
                    .map_err(|v| std::io::Error::other(v.to_string()))?;

                let is_sr = bed_record.name.to_lowercase().contains("simple_repeat");
                let coverage: usize = bed_record
                    .block_starts
                    .iter()
                    .zip(bed_record.block_sizes)
                    .filter(|v| v.0 >= &0)
                    .map(|v| v.1.abs() as usize)
                    .sum();
                let joins = bed_record
                    .block_starts
                    .iter()
                    .filter(|&v| v >= &0)
                    .count()
                    .saturating_sub(1);

                let entry = rm_results
                    .entry(bed_record.chrom)
                    .or_insert_with(TargetInfo::default);

                if is_sr {
                    entry.sr_count += 1;
                    entry.sr_coverage += coverage;
                } else {
                    entry.te_count += 1;
                    entry.te_coverage += coverage;
                }
                entry.joins += joins;
                entry.total_coverage += coverage;
            }
        }

        table_stats.extend(
            rm_results
                .iter()
                .map(|(k, v)| to_table_row(&format!("{} Repeat Masker", k), v)),
        );

        table_stats.push(to_table_row(
            "Totals Repeat Masker",
            &summary(rm_results.values()),
        ));
    }

    write_statistics_table_page(
        stats_writer,
        "Target Statistics",
        &[
            "Target_Name_string",
            "TE_Count_int",
            "TE_Coverage_(Bases)_int",
            "Simple_Repeat_Count_int",
            "Simple_Repeat_Coverage_(Bases)_int",
            "Joins_int",
            "Total_Coverage_int",
        ],
        &table_stats,
    )?;

    Ok(())
}

pub fn write_family_statistics(
    stats_writer: &mut impl Write,
    region_annotations: &[(usize, Vec<AmbiguousAnnotation>)],
    query_lengths: &HashMap<usize, usize>,
) -> std::io::Result<()> {
    let mut family_stats = HashMap::<&String, FamilyInfo>::new();

    for (_region, annotations) in region_annotations.iter() {
        let mut prior_elem = HashMap::<usize, (&AmbiguousAnnotation, usize)>::new();

        for amb_annot in annotations.iter() {
            for (query_idx, query) in amb_annot.annotations.iter().enumerate() {
                let name = &query.query_name;
                let target_range = query.target_end - query.target_start;
                let (join_cons, nojoin_cons, join_target, nojoin_target, join_div, nojoin_div) =
                    if let Some(&(prior_annot, prior_query_idx)) = prior_elem.get(&query.query_id) {
                        let block_c = Block::from_annotation(amb_annot, query_idx, 0);
                        let block_p = Block::from_annotation(prior_annot, prior_query_idx, 0);
                        let q_len = *query_lengths.get(&query.query_id).unwrap_or(&1);

                        let (c_dist, _join_type) = relative_consensus_distance(
                            &block_p,
                            &block_c,
                            ConsensusDistanceNormalization::WithLength(q_len),
                        );
                        let t_dist = block_target_distance(&block_p, &block_c) as f64;
                        let d_dist = (block_c.kimura80 - block_p.kimura80).abs();

                        if prior_annot.join_id == amb_annot.join_id {
                            (Some(c_dist), None, Some(t_dist), None, Some(d_dist), None)
                        } else {
                            (None, Some(c_dist), None, Some(t_dist), None, Some(d_dist))
                        }
                    } else {
                        (None, None, None, None, None, None)
                    };

                let f = family_stats.entry(name).or_insert_with(FamilyInfo::default);

                f.occurrences += 1;
                f.coverage += target_range;
                f.kimura80_values.0.push(query.kimura80);
                f.join_consensus_dist.0.extend(join_cons.iter());
                f.nojoin_consensus_dist.0.extend(nojoin_cons.iter());
                f.join_target_dist.0.extend(join_target.iter());
                f.nojoin_target_dist.0.extend(nojoin_target.iter());
                f.join_kimura_dist.0.extend(join_div.iter());
                f.nojoin_kimura_dist.0.extend(nojoin_div.iter());

                prior_elem
                    .entry(query.query_id)
                    .insert_entry((amb_annot, query_idx));
            }
        }
    }

    write_statistics_table_page(
        stats_writer,
        "Family Statistics",
        &[
            "Family_string",
            "Occurrences_int",
            "Coverage_int",
            "Kimura80_KDE_violin",
            "Joined_Consensus_Distance_violin",
            "Unjoined_Consensus_Distance_violin",
            "Joined_Target_Distance_violin",
            "Unjoined_Target_Distance_violin",
            "Joined_Kimura_Difference_violin",
            "Unjoined_Kimura_Difference_violin",
        ],
        &family_stats
            .iter()
            .sorted_by(|v1, v2| v2.1.cmp(v1.1))
            .map(|(k, v)| {
                [
                    k.to_string(),
                    v.occurrences.to_string(),
                    v.coverage.to_string(),
                    v.kimura80_values.0.iter().join(":"),
                    v.join_consensus_dist.0.iter().join(":"),
                    v.nojoin_consensus_dist.0.iter().join(":"),
                    v.join_target_dist.0.iter().join(":"),
                    v.nojoin_target_dist.0.iter().join(":"),
                    v.join_kimura_dist.0.iter().join(":"),
                    v.nojoin_kimura_dist.0.iter().join(":"),
                ]
            })
            .collect_vec(),
    )?;

    Ok(())
}

pub fn write_inversion_statistics(
    stats_writer: &mut impl std::io::Write,
    region_annotations: &[(usize, Vec<AmbiguousAnnotation>)],
) -> std::io::Result<()> {
    let mut inversion_stats = HashMap::<usize, InversionInfo>::new();

    for (region, annotations) in region_annotations.iter() {
        let mut inversions = InversionInfo {
            inversions: 0,
            normal_joins: 0,
        };

        let mut prior_strand_val = HashMap::<usize, Strand>::new();

        for amb_annot in annotations.iter() {
            if let Some(strand) = amb_annot.annotations.first().map(|v| &v.strand) {
                if let Some(prior_strand) = prior_strand_val.get(&amb_annot.join_id) {
                    if strand != prior_strand {
                        inversions.inversions += 1;
                    } else {
                        inversions.normal_joins += 1;
                    }
                }

                prior_strand_val.insert(amb_annot.join_id, *strand);
            }
        }

        inversion_stats.insert(*region, inversions);
    }

    write_statistics_table_page(
        stats_writer,
        "Inversion Statistics",
        &["Region_region", "Inversions_int", "Normal_Joins_int"],
        &inversion_stats
            .iter()
            .sorted_by(|v1, v2| v2.1.cmp(v1.1))
            .map(|(k, v)| {
                [
                    k.to_string(),
                    v.inversions.to_string(),
                    v.normal_joins.to_string(),
                ]
            })
            .collect_vec(),
    )?;

    Ok(())
}
