use crate::{
    alignment::Strand,
    annotation::AmbiguousAnnotation,
    assembly::{
        block_target_distance, relative_consensus_distance, ConsensusDistanceNormalization,
    },
    segments::{Block, Unordered},
    viz::SodaVizWriter,
};
use core::str;
use itertools::Itertools;
use std::{collections::HashMap, io::Write};

const TABLE_HTML: &str = include_str!("../../fixtures/soda/table.html");

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
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

                family_stats
                    .entry(name)
                    .and_modify(|f| {
                        f.occurrences += 1;
                        f.coverage += target_range;
                        f.kimura80_values.0.push(query.kimura80);
                        f.join_consensus_dist.0.extend(join_cons.iter());
                        f.nojoin_consensus_dist.0.extend(nojoin_cons.iter());
                        f.join_target_dist.0.extend(join_target.iter());
                        f.nojoin_target_dist.0.extend(nojoin_target.iter());
                        f.join_kimura_dist.0.extend(join_div.iter());
                        f.nojoin_kimura_dist.0.extend(nojoin_div.iter());
                    })
                    .or_insert(FamilyInfo {
                        occurrences: 1,
                        coverage: target_range,
                        kimura80_values: Unordered(vec![query.kimura80]),
                        join_consensus_dist: Unordered(join_cons.iter().copied().collect()),
                        nojoin_consensus_dist: Unordered(nojoin_cons.iter().copied().collect()),
                        join_target_dist: Unordered(join_target.iter().copied().collect()),
                        nojoin_target_dist: Unordered(nojoin_target.iter().copied().collect()),
                        join_kimura_dist: Unordered(join_div.iter().copied().collect()),
                        nojoin_kimura_dist: Unordered(nojoin_div.iter().copied().collect()),
                    });

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
