use itertools::Itertools;

use crate::{
    alignment::{Alignment, AlignmentData, Strand},
    chunks::ProximityGroup,
    viz::write_soda_html,
    Args,
};

///
///
///
///
struct AlignmentTuple<'a> {
    row_idx: usize,
    alignment: &'a Alignment,
    confidence: f64,
}

///
///
///
///
#[derive(PartialEq, Clone, Copy, Debug)]
enum Direction {
    Left,
    Right,
}

///
///
///
///
#[derive(PartialEq, Clone, Copy, Debug)]
struct Edge {
    idx_to: usize,
    weight: f64,
    direction: Direction,
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

            let consensus_distance =
                b.alignment.query_start as isize - a.alignment.query_end as isize;

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

///
///
///
///
pub fn collapse(
    group: &ProximityGroup,
    confidence_avg_by_row: &[f64],
    region_idx: usize,
    args: &Args,
) {
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
                .filter(|t| t.alignment.query_id == *id)
        })
        .for_each(|alignment_tuples| {
            // split the forward and reverse stranded alignments
            let (fwd, rev): (Vec<_>, Vec<_>) = alignment_tuples
                .into_iter()
                .partition(|a| a.alignment.strand == Strand::Forward);

            // graph forward, reverse
            // graph jsons before messing with them
            // assembly algorithm should be the same for forward/reverse

            // ---------
            // viz stuff
            // ---------
            let links = fwd
                .iter()
                .enumerate()
                .flat_map(|(a_idx, a_t)| {
                    graph[a_idx]
                        .iter()
                        .map(|edge| &fwd[edge.idx_to])
                        .map(|b_t| {
                            let start = (a_t.alignment.target_start + a_t.alignment.target_end) / 2;
                            let end = (b_t.alignment.target_start + b_t.alignment.target_end) / 2;
                            serde_json::json!({
                                "id": format!("{}-{}", a_t.row_idx, b_t.row_idx),
                                "start": start,
                                "end": end,
                            })
                        })
                })
                .collect_vec();
            // ---------
            // ---------

            // take all indices of the remaining alignment tuples
            let mut remaining_indices = (0..fwd.len()).collect_vec();

            // sort them by their respective confidence values
            remaining_indices.sort_by(|a, b| {
                fwd[*a]
                    .confidence
                    .partial_cmp(&fwd[*b].confidence)
                    .expect("failed to compare confidence values")
            });

            // we'll keep track of the assembled
            // alignments as vectors of those indices
            let mut assemblies: Vec<Vec<usize>> = vec![];

            // while we have alignments left to place
            while !remaining_indices.is_empty() {
                // start an assembly with the highest remaining confidence
                let assembly_root = remaining_indices.pop().unwrap();

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

            debug_assert!({
                let cnt = assemblies.iter().map(|a| a.len()).sum::<usize>();
                fwd.len() == cnt
            });

            // ---------
            // viz stuff
            // ---------

            let fwd_consensus_ali = fwd
                .iter()
                .map(|t| {
                    serde_json::json!({
                        "id": t.row_idx.to_string(),
                        "start": t.alignment.query_start,
                        "end": t.alignment.query_end,
                        "strand": t.alignment.strand,
                        "conf": t.confidence,
                    })
                })
                .collect_vec();

            let fwd_target_ali = fwd
                .iter()
                .map(|t| {
                    serde_json::json!({
                        "id": t.row_idx.to_string(),
                        "start": t.alignment.target_start,
                        "end": t.alignment.target_end,
                        "strand": t.alignment.strand,
                        "conf": t.confidence,
                    })
                })
                .collect_vec();

            let assemblies_target = assemblies
                .iter()
                .map(|a| {
                    a.iter()
                        .map(|idx| {
                            serde_json::json!({
                                "id": fwd[*idx].row_idx.to_string(),
                                "start": fwd[*idx].alignment.target_start,
                                "end": fwd[*idx].alignment.target_end,
                                "conf": fwd[*idx].confidence,
                            })
                        })
                        .collect_vec()
                })
                .collect_vec();

            let assemblies_consensus = assemblies
                .iter()
                .map(|a| {
                    a.iter()
                        .map(|idx| {
                            serde_json::json!({
                                "id": fwd[*idx].row_idx.to_string(),
                                "start": fwd[*idx].alignment.query_start,
                                "end": fwd[*idx].alignment.query_end,
                                "conf": fwd[*idx].confidence,
                            })
                        })
                        .collect_vec()
                })
                .collect_vec();

            let query_id_idx = query_ids
                .iter()
                .position(|id| id == query_id)
                .expect("failed to find query_id");

            let next_idx = query_id_idx + 1;

            let next = if next_idx < query_ids.len() {
                query_ids[next_idx]
            } else {
                query_ids[0]
            };

            let prev_idx = query_id_idx as isize - 1;
            let prev = if prev_idx > 0 {
                query_ids[prev_idx as usize]
            } else {
                query_ids[query_ids.len() - 1]
            };

            let fwd_data = serde_json::json!({
                "consensus": fwd_consensus_ali,
                "target": fwd_target_ali,
                "assembliesConsensus": assemblies_consensus,
                "assembliesTarget": assemblies_target,
                "links": links,
                "prev": prev,
                "next": next,
                "suffix": "fwd",
            });

            let html_path = viz_path.join(format!("{}-fwd.html", query_id));

            write_soda_html(
                &fwd_data,
                "./fixtures/soda/assembly-template.html",
                "./fixtures/soda/assembly.js",
                html_path,
            );
        });
}
