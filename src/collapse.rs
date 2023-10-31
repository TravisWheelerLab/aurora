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
pub struct AlignmentTuple<'a> {
    pub row_idx: usize,
    pub alignment: &'a Alignment,
    pub confidence: f64,
}

///
///
///
///
#[derive(PartialEq, Clone, Copy, Debug)]
pub enum Direction {
    Left,
    Right,
}

///
///
///
///
#[derive(PartialEq, Clone, Copy, Debug)]
pub struct Edge {
    pub idx_to: usize,
    pub weight: f64,
    pub direction: Direction,
}

pub fn assembly_graph(tuples: &[AlignmentTuple], args: &Args) -> Vec<Vec<Edge>> {
    // this relies on the alignments being sorted by target start
    tuples.iter().zip(tuples.iter().skip(1)).for_each(|(a, b)| {
        debug_assert!(a.alignment.target_start <= b.alignment.target_start);
    });

    let mut graph: Vec<Vec<Edge>> = vec![vec![]; tuples.len()];
    (0..tuples.len()).for_each(|a_idx| {
        (a_idx + 1..tuples.len()).for_each(|b_idx| {
            let a = &tuples[a_idx];
            let b = &tuples[b_idx];

            if a.alignment == b.alignment {
                return;
            }

            let within_target_distance = b
                .alignment
                .target_start
                .saturating_sub(a.alignment.target_end)
                < args.target_join_distance;

            let consensus_distance = match a.alignment.strand {
                Strand::Forward => {
                    b.alignment.query_start as isize - a.alignment.query_end as isize
                }
                Strand::Reverse => {
                    a.alignment.query_end as isize - b.alignment.query_start as isize
                }
                Strand::Unset => panic!(),
            };

            // TODO: PARAMETERIZE THIS
            let consensus_colinear = consensus_distance > -20;

            if within_target_distance && consensus_colinear {
                graph[a_idx].push(Edge {
                    idx_to: b_idx,
                    weight: consensus_distance.abs() as f64,
                    direction: Direction::Right,
                });
                graph[b_idx].push(Edge {
                    idx_to: a_idx,
                    weight: consensus_distance.abs() as f64,
                    direction: Direction::Left,
                });
            }
        });
    });

    graph
}

pub fn assembly(graph: &mut [Vec<Edge>], tuples: &[AlignmentTuple]) -> Vec<Vec<usize>> {
    // take all indices of the remaining alignment tuples
    let mut remaining_indices = (0..tuples.len()).collect_vec();

    // sort them by their respective confidence values
    remaining_indices.sort_by(|a, b| {
        tuples[*a]
            .confidence
            .partial_cmp(&tuples[*b].confidence)
            .expect("failed to compare confidence values")
    });

    // we'll keep track of the assembled
    // alignments as vectors of those indices
    let mut assemblies: Vec<Vec<usize>> = vec![];

    // while we have alignments left to place, start an
    // assembly with the highest remaining confidence
    while let Some(assembly_root) = remaining_indices.pop() {
        // remove the assembly root from all outgoing edges
        graph.iter_mut().for_each(|edge_list| {
            edge_list.retain(|e| e.idx_to != assembly_root);
        });

        let mut assembly = vec![assembly_root];

        loop {
            // the alignment we started with is the
            // one that's going to maintain edges
            let root_idx = assembly[0];

            // find the edge with the minimum weight
            let join_idx = match graph[root_idx].iter().min_by(|a, b| {
                a.weight
                    .partial_cmp(&b.weight)
                    .expect("failed to compare edge weights")
            }) {
                Some(edge) => edge.idx_to,
                // if we don't have an edge,
                // we're done with this assembly
                None => break,
            };

            // add the joined alignment to the assembly
            assembly.push(join_idx);

            // remove the joined alignment from further consideration
            match remaining_indices.iter().position(|&idx| idx == join_idx) {
                Some(pos) => remaining_indices.remove(pos),
                None => panic!("failed to remove index from remaining"),
            };

            // remove the joined alignment from all outgoing edges
            graph.iter_mut().for_each(|edge_list| {
                edge_list.retain(|e| e.idx_to != join_idx);
            });

            let root_is_on_left = root_idx < join_idx;

            if root_is_on_left {
                // if the root is on the left, we need to
                // keep all of it's left-pointing edges
                graph[root_idx].retain(|e| e.direction == Direction::Left);

                // and keep all of the right-pointing joined edges
                let mut new_edges: Vec<Edge> = graph[join_idx]
                    .iter()
                    .filter(|e| e.direction == Direction::Right)
                    .copied()
                    .collect_vec();

                graph[root_idx].append(&mut new_edges);
            } else {
                // if the root is on the right, we need to
                // keep all of it's right-pointing edges
                graph[root_idx].retain(|e| e.direction == Direction::Right);

                // and keep all of the left-pointing joined edges
                let mut new_edges: Vec<Edge> = graph[join_idx]
                    .iter()
                    .filter(|e| e.direction == Direction::Left)
                    .copied()
                    .collect_vec();

                graph[root_idx].append(&mut new_edges);
            };

            // we remove ALL edges from the joined alignment
            graph[join_idx] = vec![];
        }
        assemblies.push(assembly);
    }

    assemblies
}

///
///
///
///
pub fn collapse(
    group: &ProximityGroup,
    confidence_avg_by_row: &[f64],
    args: &Args,
) -> Vec<Vec<usize>> {
    let mut assemblies = vec![];

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
        .map(|id| {
            (
                id,
                group
                    .alignments
                    .iter()
                    // enumerate to get the row
                    .enumerate()
                    .map(|(idx, ali)| AlignmentTuple {
                        row_idx: idx + 1,
                        alignment: ali,
                        confidence: confidence_avg_by_row[idx + 1],
                    })
                    .filter(|t| t.alignment.query_id == *id),
            )
        })
        .for_each(|(query_id, alignment_tuples)| {
            // split the forward and reverse stranded alignments
            let (fwd_tuples, rev_tuples): (Vec<_>, Vec<_>) = alignment_tuples
                .into_iter()
                .partition(|a| a.alignment.strand == Strand::Forward);

            let mut fwd_graph = assembly_graph(&fwd_tuples, args);
            let mut rev_graph = assembly_graph(&rev_tuples, args);

            // if we are going to generate soda output
            // for the assemblies, we need to store the
            // links before we start messing with the graph
            let mut fwd_links = vec![];
            let mut rev_links = vec![];
            if args.viz {
                fwd_links = fwd_tuples
                    .iter()
                    .enumerate()
                    .flat_map(|(a_idx, a_t)| {
                        fwd_graph[a_idx]
                            .iter()
                            .map(|edge| &fwd_tuples[edge.idx_to])
                            .map(|b_t| {
                                let start =
                                    (a_t.alignment.target_start + a_t.alignment.target_end) / 2;
                                let end =
                                    (b_t.alignment.target_start + b_t.alignment.target_end) / 2;
                                format!("{}-{},{},{}", a_t.row_idx, b_t.row_idx, start, end)
                            })
                    })
                    .collect_vec();

                rev_links = rev_tuples
                    .iter()
                    .enumerate()
                    .flat_map(|(a_idx, a_t)| {
                        rev_graph[a_idx]
                            .iter()
                            .map(|edge| &rev_tuples[edge.idx_to])
                            .map(|b_t| {
                                let start =
                                    (a_t.alignment.target_start + a_t.alignment.target_end) / 2;
                                let end =
                                    (b_t.alignment.target_start + b_t.alignment.target_end) / 2;
                                format!("{}-{},{},{}", a_t.row_idx, b_t.row_idx, start, end)
                            })
                    })
                    .collect_vec();
            }

            let fwd_assemblies = assembly(&mut fwd_graph, &fwd_tuples);
            let rev_assemblies = assembly(&mut rev_graph, &rev_tuples);

            // check that the sum of assembly lengths is equal
            // to the number of alignments we started with
            debug_assert!({
                let cnt = fwd_assemblies.iter().map(|a| a.len()).sum::<usize>();
                fwd_tuples.len() == cnt
            });

            debug_assert!({
                let cnt = rev_assemblies.iter().map(|a| a.len()).sum::<usize>();
                rev_tuples.len() == cnt
            });

            if args.viz {
                if !fwd_tuples.is_empty() {
                    let fwd_data = AuroraAssemblySodaData::new(
                        &fwd_tuples,
                        &fwd_assemblies,
                        &query_ids,
                        fwd_links,
                    );

                    let fwd_path = args.viz_output_path.join(format!("{}-fwd.html", query_id));

                    write_soda_html(
                        &fwd_data,
                        "./fixtures/soda/assembly-template.html",
                        "./fixtures/soda/assembly.js",
                        fwd_path,
                    );
                }

                if !rev_tuples.is_empty() {
                    let rev_data = AuroraAssemblySodaData::new(
                        &rev_tuples,
                        &rev_assemblies,
                        &query_ids,
                        rev_links,
                    );

                    let rev_path = args.viz_output_path.join(format!("{}-rev.html", query_id));

                    write_soda_html(
                        &rev_data,
                        "./fixtures/soda/assembly-template.html",
                        "./fixtures/soda/assembly.js",
                        rev_path,
                    );
                }
            }

            // map the tuple-index-based assemblies
            // to row-index-based assemblies
            assemblies.append(
                &mut fwd_assemblies
                    .iter()
                    .map(|a| {
                        a.iter()
                            .map(|&tuple_idx| &fwd_tuples[tuple_idx])
                            .map(|t| t.row_idx)
                            .collect_vec()
                    })
                    .collect_vec(),
            );
            assemblies.append(
                &mut rev_assemblies
                    .iter()
                    .map(|a| {
                        a.iter()
                            .map(|&tuple_idx| &rev_tuples[tuple_idx])
                            .map(|t| t.row_idx)
                            .collect_vec()
                    })
                    .collect_vec(),
            );
        });

    assemblies
}
