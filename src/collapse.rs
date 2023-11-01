use std::collections::HashMap;

use itertools::Itertools;

use crate::{
    alignment::{Alignment, Strand},
    chunks::ProximityGroup,
    viz::{write_soda_html, AuroraAssemblySodaData},
    Args,
};

///
///
///
///
#[derive(PartialEq, Clone, Copy, Debug)]
pub enum Direction {
    Left,
    Right,
}

#[derive(PartialEq, Clone, Copy, Debug)]
pub struct Edge<'a> {
    pub alignment_to: &'a Alignment,
    pub weight: f64,
    pub direction: Direction,
}

pub fn assembly_graph<'a>(
    alignments: &[&'a Alignment],
    args: &Args,
) -> HashMap<&'a Alignment, Vec<Edge<'a>>> {
    // this relies on the alignments being sorted by target start
    alignments
        .iter()
        .zip(alignments.iter().skip(1))
        .for_each(|(a, b)| {
            debug_assert!(a.target_start <= b.target_start);
        });

    let mut hash_graph: HashMap<&Alignment, Vec<Edge>> = HashMap::new();

    alignments.iter().enumerate().for_each(|(a_idx, &a)| {
        alignments[a_idx + 1..].iter().for_each(|&b| {
            if a == b {
                return;
            }

            let within_target_distance =
                b.target_start.saturating_sub(a.target_end) < args.target_join_distance;

            let consensus_distance = match a.strand {
                Strand::Forward => b.query_start as isize - a.query_end as isize,
                Strand::Reverse => a.query_end as isize - b.query_start as isize,
                Strand::Unset => panic!(),
            };

            // TODO: PARAMETERIZE THIS
            let consensus_colinear = consensus_distance > -20;

            let weight = consensus_distance.abs() as f64;

            if within_target_distance && consensus_colinear {
                let a_edges = hash_graph.entry(a).or_default();
                a_edges.push(Edge {
                    alignment_to: b,
                    weight,
                    direction: Direction::Right,
                });

                let b_edges = hash_graph.entry(b).or_default();
                b_edges.push(Edge {
                    alignment_to: a,
                    weight,
                    direction: Direction::Right,
                })
            }
        });
    });

    hash_graph
}

pub fn assembly<'a>(
    graph: &mut HashMap<&'a Alignment, Vec<Edge<'a>>>,
    confidence_avg_by_id: &HashMap<usize, f64>,
) -> Vec<Assembly<'a>> {
    // take all indices of the remaining alignment tuples
    let mut remaining: Vec<&Alignment> = graph.keys().copied().collect_vec();

    // sort them by their respective confidence values
    remaining.sort_by(|a, b| {
        confidence_avg_by_id
            .get(&a.id)
            .unwrap()
            .partial_cmp(confidence_avg_by_id.get(&b.id).unwrap())
            .expect("failed to compare confidences")
    });

    let mut assemblies: Vec<Assembly> = vec![];

    // while we have alignments left to place, start an
    // assembly with the highest remaining confidence
    while let Some(assembly_root) = remaining.pop() {
        // remove the assembly root from all outgoing edges
        graph.values_mut().for_each(|edge_list| {
            edge_list.retain(|e| e.alignment_to != assembly_root);
        });

        let mut assembly_alignments = vec![assembly_root];

        loop {
            // the alignment we started with is the
            // one that's going to maintain edges
            let root_ali = assembly_alignments[0];

            // find the edge with the minimum weight
            let join_ali = match graph[root_ali].iter().min_by(|edge_a, edge_b| {
                edge_a
                    .weight
                    .partial_cmp(&edge_b.weight)
                    .expect("failed to compare edge weights")
            }) {
                Some(edge) => edge.alignment_to,
                // if we don't have an edge,
                // we're done with this assembly
                None => break,
            };

            // add the joined alignment to the assembly
            assembly_alignments.push(join_ali);

            // remove the joined alignment from further consideration
            match remaining.iter().position(|&ali| ali == join_ali) {
                Some(pos) => remaining.remove(pos),
                None => panic!("failed to remove index from remaining"),
            };

            // remove the joined alignment from all outgoing edges
            graph.values_mut().for_each(|edge_list| {
                edge_list.retain(|e| e.alignment_to != join_ali);
            });

            let root_is_on_left = root_ali.target_start < join_ali.target_start;

            if root_is_on_left {
                // if the root is on the left, we need to
                // keep all of it's left-pointing edges
                graph
                    .get_mut(root_ali)
                    .unwrap()
                    .retain(|e| e.direction == Direction::Left);

                // and keep all of the right-pointing joined edges
                let mut new_edges: Vec<Edge> = graph[join_ali]
                    .iter()
                    .filter(|e| e.direction == Direction::Right)
                    .copied()
                    .collect_vec();

                graph.get_mut(root_ali).unwrap().append(&mut new_edges);
            } else {
                // if the root is on the right, we need to
                // keep all of it's right-pointing edges
                graph
                    .get_mut(root_ali)
                    .unwrap()
                    .retain(|e| e.direction == Direction::Right);

                // and keep all of the left-pointing joined edges
                let mut new_edges: Vec<Edge> = graph[join_ali]
                    .iter()
                    .filter(|e| e.direction == Direction::Left)
                    .copied()
                    .collect_vec();

                graph.get_mut(root_ali).unwrap().append(&mut new_edges);
            };
            // we remove ALL edges from the joined alignment
            graph.insert(join_ali, vec![]);
        }
        assemblies.push(Assembly::new(assembly_alignments));
    }

    assemblies
}

///
///
///
///
pub struct Assembly<'a> {
    pub alignments: Vec<&'a Alignment>,
    pub alignment_ranges: Vec<[usize; 2]>,
}

impl<'a> Assembly<'a> {
    fn new(alignments: Vec<&'a Alignment>) -> Self {
        Self {
            alignments,
            alignment_ranges: vec![],
        }
    }
}

///
///
///
///
pub struct AssemblyGroup<'a> {
    pub target_id: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub assemblies: Vec<Assembly<'a>>,
}

impl<'a> AssemblyGroup<'a> {
    pub fn new(
        group: &ProximityGroup,
        confidence_avg_by_id: &HashMap<usize, f64>,
        args: &Args,
    ) -> Self {
        let mut assemblies: Vec<Assembly> = vec![];

        let mut query_ids: Vec<usize> = group
            .alignments
            .iter()
            .map(|a| a.query_id)
            .unique()
            .collect();

        query_ids.sort();

        query_ids
            .iter()
            // grab the alignments for this ID
            .map(|id| (id, group.alignments.iter().filter(|a| a.query_id == *id)))
            .for_each(|(query_id, alignments)| {
                // split the forward and reverse stranded alignments
                let (fwd_ali, rev_ali): (Vec<&Alignment>, Vec<&Alignment>) = alignments
                    .into_iter()
                    .partition(|a| a.strand == Strand::Forward);

                let mut fwd_graph = assembly_graph(&fwd_ali, args);
                let mut rev_graph = assembly_graph(&rev_ali, args);

                // if we are going to generate soda output
                // for the assemblies, we need to store the
                // links before we start messing with the graph
                // let mut fwd_links = vec![];
                // let mut rev_links = vec![];
                // if args.viz {
                //     fwd_links = fwd_ali
                //         .iter()
                //         .enumerate()
                //         .flat_map(|(a_idx, a_t)| {
                //             fwd_graph[a_idx]
                //                 .iter()
                //                 .map(|edge| &fwd_ali[edge.idx_to])
                //                 .map(|b_t| {
                //                     let start = (a_t.target_start + a_t.target_end) / 2;
                //                     let end = (b_t.target_start + b_t.target_end) / 2;
                //                     format!("{}-{},{},{}", a_t.id, b_t.id, start, end)
                //                 })
                //         })
                //         .collect_vec();

                //     rev_links = rev_ali
                //         .iter()
                //         .enumerate()
                //         .flat_map(|(a_idx, a_t)| {
                //             rev_graph[a_idx]
                //                 .iter()
                //                 .map(|edge| &rev_ali[edge.idx_to])
                //                 .map(|b_t| {
                //                     let start = (a_t.target_start + a_t.target_end) / 2;
                //                     let end = (b_t.target_start + b_t.target_end) / 2;
                //                     format!("{}-{},{},{}", a_t.id, b_t.id, start, end)
                //                 })
                //         })
                //         .collect_vec();
                // }

                let fwd_assemblies = assembly(&mut fwd_graph, confidence_avg_by_id);
                let rev_assemblies = assembly(&mut rev_graph, confidence_avg_by_id);

                // check that the sum of assembly lengths is equal
                // to the number of alignments we started with
                debug_assert!({
                    let cnt = fwd_assemblies
                        .iter()
                        .map(|a| a.alignments.len())
                        .sum::<usize>();
                    fwd_ali.len() == cnt
                });

                debug_assert!({
                    let cnt = rev_assemblies
                        .iter()
                        .map(|a| a.alignments.len())
                        .sum::<usize>();
                    rev_ali.len() == cnt
                });

                // if args.viz {
                //     if !fwd_ali.is_empty() {
                //         let fwd_data = AuroraAssemblySodaData::new(
                //             &fwd_ali,
                //             &fwd_assemblies,
                //             &query_ids,
                //             fwd_links,
                //         );

                //         let fwd_path = args.viz_output_path.join(format!("{}-fwd.html", query_id));

                //         write_soda_html(
                //             &fwd_data,
                //             "./fixtures/soda/assembly-template.html",
                //             "./fixtures/soda/assembly.js",
                //             fwd_path,
                //         );
                //     }

                //     if !rev_ali.is_empty() {
                //         let rev_data = AuroraAssemblySodaData::new(
                //             &rev_ali,
                //             &rev_assemblies,
                //             &query_ids,
                //             rev_links,
                //         );

                //         let rev_path = args.viz_output_path.join(format!("{}-rev.html", query_id));

                //         write_soda_html(
                //             &rev_data,
                //             "./fixtures/soda/assembly-template.html",
                //             "./fixtures/soda/assembly.js",
                //             rev_path,
                //         );
                //     }
                // }

                // // map the tuple-index-based assemblies
                // // to row-index-based assemblies
                // assemblies.append(
                //     &mut fwd_assemblies
                //         .iter()
                //         .map(|a| {
                //             a.iter()
                //                 .map(|&tuple_idx| &fwd_ali[tuple_idx])
                //                 .map(|t| t.row_idx)
                //                 .collect_vec()
                //         })
                //         .collect_vec(),
                // );

                // assemblies.append(
                //     &mut rev_assemblies
                //         .iter()
                //         .map(|a| {
                //             a.iter()
                //                 .map(|&tuple_idx| &rev_ali[tuple_idx])
                //                 .map(|t| t.row_idx)
                //                 .collect_vec()
                //         })
                //         .collect_vec(),
                // );
            });

        Self {
            target_id: group.target_id,
            target_start: group.target_start,
            target_end: group.target_end,
            assemblies,
        }
    }
}
