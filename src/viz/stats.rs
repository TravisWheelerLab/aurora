use crate::{
    alignment::{AlignmentData, Strand},
    annotation::AmbiguousAnnotation,
};
use itertools::Itertools;
use std::collections::HashMap;

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
struct FamilyInfo {
    pub occurrences: usize,
    pub coverage: usize,
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord)]
struct InversionInfo {
    pub inversions: usize,
    pub normal_joins: usize,
}

pub fn write_family_statistics(
    stats_writer: &mut impl std::io::Write,
    region_annotations: &[(usize, Vec<AmbiguousAnnotation>)],
    alignment_data: &AlignmentData,
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
                    })
                    .or_insert(FamilyInfo {
                        occurrences: 1,
                        coverage: target_range,
                    });
            }
        }
    }

    writeln!(stats_writer, "<html><body><style>table {{ border-collapse: collapse; border: 2px solid rgb(140 140 140); font-family: sans-serif; font-size: 0.8rem; letter-spacing: 1px; }}\nth, td {{ border: 1px solid rgb(160 160 160); padding: 8px 10px; }}</style><table><thead><tr><th>Family</th><th>Occurrences</th><th>Coverage</th></thead><tbody>")?;

    for (k, v) in family_stats.iter().sorted_by(|v1, v2| v2.1.cmp(v1.1)) {
        writeln!(
            stats_writer,
            "<tr><td>{}</td><td>{}</td><td>{}</td></tr>",
            k, v.occurrences, v.coverage
        )?;
    }

    writeln!(stats_writer, "</tbody></body></html>")?;

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

    writeln!(stats_writer, "<html><body><style>table {{ border-collapse: collapse; border: 2px solid rgb(140 140 140); font-family: sans-serif; font-size: 0.8rem; letter-spacing: 1px; }}\nth, td {{ border: 1px solid rgb(160 160 160); padding: 8px 10px; }}</style><table><thead><tr><th>Region</th><th>Inversions</th><th>Normal Joins</th></thead><tbody>")?;

    for (k, v) in inversion_stats.iter().sorted_by(|v1, v2| v2.1.cmp(v1.1)) {
        writeln!(
            stats_writer,
            "<tr><td><a href=\"{}/index.html\">{}</a></td><td>{}</td><td>{}</td></tr>",
            k, k, v.inversions, v.normal_joins
        )?;
    }

    writeln!(stats_writer, "</tbody></body></html>")?;

    Ok(())
}
