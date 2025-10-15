use std::{
    collections::{HashMap, HashSet},
    fs::File,
    hash::Hash,
    io::{BufWriter, Write},
};

use itertools::Itertools;

use crate::{
    alignment::{Alignment, AlignmentData, Strand},
    chunks::ProximityGroup,
    viz::AssemblySodaData,
    AnnotationArgs, AuroraArgs,
};

/// The direction of an `Edge` in terms of where
/// `&Alignment` B (value) is in relation to `&Alignment` A (key)
/// in the coordinate space of the chromosome
#[derive(Hash, Eq, PartialEq, Clone, Copy, Debug)]
pub enum Direction {
    Left,
    Right,
}

#[derive(Clone, Copy, Debug)]
pub struct Edge {
    pub edge_to: usize,
    pub weight: f64,
    pub direction: Direction,
}

// Hashing and equivalence only based on to edge field.
impl PartialEq for Edge {
    fn eq(&self, other: &Self) -> bool {
        self.edge_to == other.edge_to
    }
}

impl Eq for Edge {}

impl Hash for Edge {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.edge_to.hash(state);
    }
}

fn link_assemblies(
    graph: &mut [HashSet<Edge>],
    alignments: &[(usize, &Alignment)],
    args: &AnnotationArgs,
) {
    // this relies on the alignments being sorted by target start
    // note: this assertion iter will only run in debug mode
    alignments
        .iter()
        .zip(alignments.iter().skip(1))
        .for_each(|(a, b)| {
            debug_assert!(a.1.target_start <= b.1.target_start);
        });
    // We also rely on the fact that all alignment indexes are actually in the graph!
    alignments
        .iter()
        .for_each(|a| debug_assert!(a.0 < graph.len()));

    alignments
        .iter()
        .enumerate()
        .for_each(|(idx, &(a_idx, a))| {
            alignments[idx + 1..].iter().for_each(|&(b_idx, b)| {
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
                    graph[a_idx].insert(Edge {
                        edge_to: b_idx,
                        weight,
                        direction: Direction::Right,
                    });
                    graph[b_idx].insert(Edge {
                        edge_to: a_idx,
                        weight,
                        direction: Direction::Left,
                    });
                }
            });
        });
}

/// Represents graph of compatable alignments on the genome.
/// For each alignment, stores all alignments from the same query in front of it.
pub struct AssemblyGraph {
    pub fwd_graph: Vec<HashSet<Edge>>,
    pub rev_graph: Vec<HashSet<Edge>>,
}

impl AssemblyGraph {
    pub fn new(
        group: &ProximityGroup,
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

        let mut fwd_graph: Vec<HashSet<Edge>> = vec![HashSet::new(); group.alignments.len()];
        let mut rev_graph: Vec<HashSet<Edge>> = vec![HashSet::new(); group.alignments.len()];

        let mut query_files: Vec<String> = vec![];
        let mut file_query_ids: Vec<usize> = vec![];
        let mut file_strand: Vec<bool> = vec![];

        query_ids
            .iter()
            // grab the alignments for this ID
            .map(|id| {
                (
                    id,
                    group
                        .alignments
                        .iter()
                        .enumerate()
                        .filter(|&(_, a)| a.query_id == *id),
                )
            })
            .for_each(|(query_id, alignments)| {
                // split the forward and reverse stranded alignments
                let (fwd_ali, rev_ali): (Vec<(usize, &Alignment)>, Vec<(usize, &Alignment)>) =
                    alignments
                        .into_iter()
                        .partition(|&(_, a)| a.strand == Strand::Forward);

                link_assemblies(&mut fwd_graph, &fwd_ali, &args.annotation_args);
                link_assemblies(&mut rev_graph, &rev_ali, &args.annotation_args);

                if args.visualization_args.assembly_viz {
                    let fwd_links = fwd_ali
                        .iter()
                        .flat_map(|&(ali_idx, ali_from)| {
                            fwd_graph[ali_idx].iter().map(|e| {
                                let ali_to = &group.alignments[e.edge_to];
                                let start = (ali_from.target_start + ali_from.target_end) / 2;
                                let end = (ali_to.target_start + ali_to.target_end) / 2;
                                format!("{}-{},{},{}", ali_from.id, ali_to.id, start, end)
                            })
                        })
                        .collect_vec();

                    let rev_links = rev_ali
                        .iter()
                        .flat_map(|&(ali_idx, ali_from)| {
                            fwd_graph[ali_idx].iter().map(|e| {
                                let ali_to = &group.alignments[e.edge_to];
                                let start = (ali_from.target_start + ali_from.target_end) / 2;
                                let end = (ali_to.target_start + ali_to.target_end) / 2;
                                format!("{}-{},{},{}", ali_from.id, ali_to.id, start, end)
                            })
                        })
                        .collect_vec();

                    if !fwd_ali.is_empty() {
                        let fwd_ali_only = fwd_ali.iter().map(|&(_, al)| al).collect_vec();
                        query_files.push(format!("{}-fwd.html", query_id));
                        file_query_ids.push(*query_id);
                        file_strand.push(true);
                        AssemblySodaData::new(
                            &fwd_ali_only,
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
                        let rev_ali_only = rev_ali.iter().map(|&(_, al)| al).collect_vec();
                        query_files.push(format!("{}-rev.html", query_id));
                        file_query_ids.push(*query_id);
                        file_strand.push(false);
                        AssemblySodaData::new(
                            &rev_ali_only,
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

        Self {
            fwd_graph,
            rev_graph,
        }
    }
}
