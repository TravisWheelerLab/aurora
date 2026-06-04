use crate::{
    alignment::Strand, annotation::AmbiguousAnnotation, segments::Unordered, viz::SODA_JS,
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
        .replace("SODA_TARGET", SODA_JS)
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
) -> std::io::Result<()> {
    let mut family_stats = HashMap::<&String, FamilyInfo>::new();

    for (_region, annotations) in region_annotations.iter() {
        for amb_annot in annotations.iter() {
            for query in amb_annot.annotations.iter() {
                let name = &query.query_name;
                let target_range = query.target_end - query.target_start;

                family_stats
                    .entry(name)
                    .and_modify(|f| {
                        f.occurrences += 1;
                        f.coverage += target_range;
                        f.kimura80_values.0.push(query.kimura80);
                    })
                    .or_insert(FamilyInfo {
                        occurrences: 1,
                        coverage: target_range,
                        kimura80_values: Unordered(vec![query.kimura80]),
                    });
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
