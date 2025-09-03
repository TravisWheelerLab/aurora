use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufWriter, Write},
};

use itertools::Itertools;

use crate::{
    alignment::{Alignment, AlignmentData, Strand, TandemRepeat},
    chunks::ProximityGroup,
    viz::AssemblySodaData,
    AnnotationArgs, AuroraArgs,
};

/// The direction of an `Edge` in terms of where
/// `&Alignment` B (value) is in relation to `&Alignment` A (key)
/// in the coordinate space of the chromosome
#[derive(PartialEq, Clone, Copy, Debug)]
pub enum Direction {
    Left,
    Right,
}

#[derive(PartialEq, Clone, Copy, Debug)]
pub struct Edge<'a> {
    pub ali_to: &'a Alignment,
    pub weight: f64,
    pub direction: Direction,
}

pub fn assembly_graph<'a>(
    alignments: &[&'a Alignment],
    args: &AnnotationArgs,
) -> HashMap<&'a Alignment, Vec<Edge<'a>>> {
    // this relies on the alignments being sorted by target start
    // note: this assertion iter will only run in debug mode
    alignments
        .iter()
        .zip(alignments.iter().skip(1))
        .for_each(|(a, b)| {
            debug_assert!(a.target_start <= b.target_start);
        });

    let mut graph: HashMap<&Alignment, Vec<Edge>> =
        alignments.iter().map(|&a| (a, vec![])).collect();

    alignments.iter().enumerate().for_each(|(a_idx, &a)| {
        alignments[a_idx + 1..].iter().for_each(|&b| {
            // TODO: this is highly suspect, as this should never happen
            //       ?????
            if a == b {
                return;
            }

            let target_distance = b.target_start as isize - a.target_end as isize;

            let within_target_distance_threshold =
                target_distance < args.target_join_distance as isize;

            let consensus_distance = match a.strand {
                Strand::Forward => b.query_start as isize - a.query_end as isize,
                Strand::Reverse => a.query_end as isize - b.query_start as isize,
                Strand::Unset => panic!(),
            };

            // TODO: PARAMETERIZE THIS
            let consensus_is_colinear = consensus_distance > -20;

            // let weight = consensus_distance.abs() as f64;
            let weight = target_distance.abs() as f64;

            if within_target_distance_threshold && consensus_is_colinear {
                let a_edges = graph.entry(a).or_default();
                a_edges.push(Edge {
                    ali_to: b,
                    weight,
                    direction: Direction::Right,
                });

                let b_edges = graph.entry(b).or_default();
                b_edges.push(Edge {
                    ali_to: a,
                    weight,
                    direction: Direction::Left,
                });
            }
        });
    });

    graph
}

pub struct AlignmentRange {
    pub ali_id: usize,
    /// The start column of the usage of the alignment relative to the assembly
    pub assembly_col_start: usize,
    /// The end column of the usage of the alignment relative to the assembly
    pub assembly_col_end: usize,
}

/// Represents graph of compatable alignments on the genome.
/// For each alignment, stores all alignments from the same query in front of it.
pub struct AssemblyGraph<'a> {
    pub fwd_map: HashMap<&'a Alignment, HashSet<&'a Alignment>>,
    pub rev_map: HashMap<&'a Alignment, HashSet<&'a Alignment>>,
}

impl<'a> AssemblyGraph<'a> {
    pub fn new(
        group: &ProximityGroup<'a>,
        confidence_avg_by_id: &HashMap<usize, f64>,
        args: &AuroraArgs,
        alignment_data: &AlignmentData,
    ) -> Self {
        let mut query_ids: Vec<usize> = group
            .alignments
            .iter()
            .map(|a| a.query_id)
            .unique()
            .collect();

        query_ids.sort();

        let mut fwd_map: HashMap<&'a Alignment, HashSet<&'a Alignment>> = HashMap::new();
        let mut rev_map: HashMap<&'a Alignment, HashSet<&'a Alignment>> = HashMap::new();

        let mut query_files: Vec<String> = vec![];
        let mut file_query_ids: Vec<usize> = vec![];
        let mut file_strand: Vec<bool> = vec![];

        query_ids
            .iter()
            // grab the alignments for this ID
            .map(|id| (id, group.alignments.iter().filter(|a| a.query_id == *id)))
            .for_each(|(query_id, alignments)| {
                // split the forward and reverse stranded alignments
                let (fwd_ali, rev_ali): (Vec<&Alignment>, Vec<&Alignment>) = alignments
                    .into_iter()
                    .partition(|a| a.strand == Strand::Forward);

                let fwd_graph = assembly_graph(&fwd_ali, &args.annotation_args);
                let rev_graph = assembly_graph(&rev_ali, &args.annotation_args);

                if args.visualization_args.assembly_viz {
                    let fwd_links = fwd_graph
                        .iter()
                        .flat_map(|(ali_from, edges)| {
                            edges.iter().map(|e| {
                                let ali_to = e.ali_to;
                                let start = (ali_from.target_start + ali_from.target_end) / 2;
                                let end = (ali_to.target_start + ali_to.target_end) / 2;
                                format!("{}-{},{},{}", ali_from.id, ali_to.id, start, end)
                            })
                        })
                        .collect_vec();

                    let rev_links = rev_graph
                        .iter()
                        .flat_map(|(ali_from, edges)| {
                            edges.iter().map(|e| {
                                let ali_to = e.ali_to;
                                let start = (ali_from.target_start + ali_from.target_end) / 2;
                                let end = (ali_to.target_start + ali_to.target_end) / 2;
                                format!("{}-{},{},{}", ali_from.id, ali_to.id, start, end)
                            })
                        })
                        .collect_vec();

                    if !fwd_ali.is_empty() {
                        query_files.push(format!("{}-fwd.html", query_id));
                        file_query_ids.push(*query_id);
                        file_strand.push(true);
                        AssemblySodaData::new(
                            &fwd_ali,
                            fwd_links,
                            confidence_avg_by_id,
                            alignment_data,
                        )
                        .write(
                            args.visualization_args
                                .viz_output_path
                                .join(query_files.last().unwrap()),
                        );
                    }

                    if !rev_ali.is_empty() {
                        query_files.push(format!("{}-rev.html", query_id));
                        file_query_ids.push(*query_id);
                        file_strand.push(false);
                        AssemblySodaData::new(
                            &rev_ali,
                            rev_links,
                            confidence_avg_by_id,
                            alignment_data,
                        )
                        .write(
                            args.visualization_args
                                .viz_output_path
                                .join(query_files.last().unwrap()),
                        );
                    }
                }

                // Extend graph of forward (on query sequence) alignments, we filter to only alignments to the right of each alignment...
                fwd_map.extend(fwd_graph.into_iter().map(|(al, edges)| {
                    return (
                        al,
                        edges
                            .iter()
                            .filter(|e| e.direction == Direction::Right)
                            .map(|e| e.ali_to)
                            .collect(),
                    );
                }));
                rev_map.extend(rev_graph.into_iter().map(|(al, edges)| {
                    return (
                        al,
                        edges
                            .iter()
                            .filter(|e| e.direction == Direction::Right)
                            .map(|e| e.ali_to)
                            .collect(),
                    );
                }));
            });

        if args.visualization_args.assembly_viz {
            let error_msg = "failed to write to assembly index file";
            let viz_args = &args.visualization_args;

            let asm_index_file =
                File::create(viz_args.viz_output_path.join("assembly_index.html")).unwrap();
            let mut asm_index_writer = BufWriter::new(asm_index_file);

            writeln!(&mut asm_index_writer, "<!doctype html>\n<html>\n<body>\n<ul>\n<a href=\"../index.html\">Back</a><br>\n<h1>Assemblies</h1>").expect(error_msg);

            file_query_ids
                .iter()
                .zip(query_files)
                .zip(file_strand)
                .for_each(|((&q_id, file_name), is_fwd)| {
                    writeln!(
                        &mut asm_index_writer,
                        "<li><a href=\"{}\">{} {} (id {})</a></li>",
                        file_name,
                        alignment_data.query_name_map.get(q_id),
                        if is_fwd { "Forward" } else { "Reverse" },
                        q_id
                    )
                    .expect(error_msg);
                });

            writeln!(&mut asm_index_writer, "</ul>\n</body>\n</html>").expect(error_msg);
        }

        return Self { fwd_map, rev_map };
    }

    /// Check if two alignments are compatable, or could possibly be connected with an insertion in the middle.
    pub fn compatable(&self, alignment: &Alignment, alignment_other: &Alignment) -> bool {
        return self.fwd_map[alignment].contains(alignment_other)
            || self.rev_map[alignment].contains(alignment_other);
    }
}
