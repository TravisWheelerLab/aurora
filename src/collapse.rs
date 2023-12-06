use std::collections::HashMap;

use itertools::Itertools;

use crate::{
    alignment::{Alignment, Strand, TandemRepeat},
    chunks::ProximityGroup,
    viz::{write_soda_html, AssemblySodaData},
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
    pub ali_to: &'a Alignment,
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

    let mut graph: HashMap<&Alignment, Vec<Edge>> =
        alignments.iter().map(|&a| (a, vec![])).collect();

    alignments.iter().enumerate().for_each(|(a_idx, &a)| {
        alignments[a_idx + 1..].iter().for_each(|&b| {
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

pub fn assembly<'a>(
    graph: &mut HashMap<&'a Alignment, Vec<Edge<'a>>>,
    confidence_by_id: &HashMap<usize, Vec<f64>>,
) -> Vec<Assembly<'a>> {
    // take all indices of the remaining alignment tuples
    let mut remaining: Vec<&Alignment> = graph.keys().copied().collect_vec();

    // sort the edge lists by edge weight
    graph
        .values_mut()
        .for_each(|edge_list| edge_list.sort_by(|a, b| a.weight.partial_cmp(&b.weight).unwrap()));

    // sort the remaining alignments by their minimum edge weights
    remaining.sort_by(|a, b| {
        let x = match graph.get(a).unwrap().first() {
            Some(e) => e.weight,
            None => f64::INFINITY,
        };

        let y = match graph.get(b).unwrap().first() {
            Some(e) => e.weight,
            None => f64::INFINITY,
        };

        x.partial_cmp(&y).unwrap()
    });

    remaining.reverse();

    let mut assemblies: Vec<Assembly> = vec![];

    // while we have alignments left to place, start an
    // assembly with the highest remaining confidence
    while let Some(assembly_root) = remaining.pop() {
        // remove the assembly root from all outgoing edges
        graph.values_mut().for_each(|edge_list| {
            edge_list.retain(|e| e.ali_to != assembly_root);
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
                Some(edge) => edge.ali_to,
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
                edge_list.retain(|e| e.ali_to != join_ali);
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

        assembly_alignments.sort_by_key(|a| a.target_start);

        assemblies.push(Assembly::new(assembly_alignments, confidence_by_id));
    }

    assemblies
}

pub struct AlignmentRange {
    pub ali_id: usize,
    /// The start column of the usage of the alignment relative to the assembly
    pub assembly_col_start: usize,
    /// The end column of the usage of the alignment relative to the assembly
    pub assembly_col_end: usize,
}

///
///
///
///
pub struct Assembly<'a> {
    pub query_id: usize,
    pub strand: Strand,
    pub target_start: usize,
    pub target_end: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub alignments: Vec<&'a Alignment>,
    pub alignment_ranges: Vec<AlignmentRange>,
}

impl<'a> Assembly<'a> {
    fn new(alignments: Vec<&'a Alignment>, confidence_by_id: &HashMap<usize, Vec<f64>>) -> Self {
        let strand = alignments[0].strand;

        let assembly_target_start = alignments.iter().map(|a| a.target_start).min().unwrap();
        let assembly_target_end = alignments.iter().map(|a| a.target_end).max().unwrap();

        let (query_start, query_end) = match strand {
            Strand::Forward => (
                alignments.iter().map(|a| a.query_start).min().unwrap(),
                alignments.iter().map(|a| a.query_end).max().unwrap(),
            ),
            Strand::Reverse => (
                // note: maintain query_start as the larger value here
                alignments.iter().map(|a| a.query_start).max().unwrap(),
                alignments.iter().map(|a| a.query_end).min().unwrap(),
            ),
            Strand::Unset => panic!(),
        };

        let mut alignment_ranges = vec![];

        if alignments.len() == 1 {
            alignment_ranges.push(AlignmentRange {
                ali_id: alignments[0].id,
                assembly_col_start: alignments[0].target_start - assembly_target_start,
                assembly_col_end: alignments[0].target_end - assembly_target_start,
            });
        } else {
            let min_target_start = alignments
                .iter()
                .map(|a| a.target_start)
                .min()
                .expect("empty alignments");

            let max_target_end = alignments
                .iter()
                .map(|a| a.target_end)
                .max()
                .expect("empty alignments");

            let target_len = max_target_end - min_target_start + 1;

            let mut matrix = vec![vec![0.0; alignments.len() + 1]; target_len];
            let mut sources = vec![vec![0usize; alignments.len() + 1]; target_len];

            sources[0].iter_mut().enumerate().for_each(|(i, s)| *s = i);

            (0..target_len).for_each(|col_idx| matrix[col_idx][0] = 1.0);

            alignments
                .iter()
                .enumerate()
                // +1 to get the right row idx (skip state offset)
                .map(|(idx, a)| (idx + 1, a))
                .for_each(|(row_idx, a)| {
                    let conf_vec = confidence_by_id.get(&a.id).expect("alignment ID not found");

                    debug_assert_eq!(a.target_end - a.target_start + 1, conf_vec.len());

                    (a.target_start..=a.target_end)
                        .map(|target_idx| target_idx - min_target_start)
                        .zip(conf_vec)
                        .for_each(|(col_idx, conf)| {
                            matrix[col_idx][row_idx] = *conf;
                            matrix[col_idx][0] = 0.0;
                        });
                });

            matrix.iter_mut().flatten().for_each(|v| *v = v.ln());

            (1..target_len).for_each(|col_idx| {
                let (row_of_max_in_prev_col, &max_score_in_prev_col) = matrix[col_idx - 1]
                    .iter()
                    .enumerate()
                    // skip the skip
                    .skip(1)
                    .max_by(|(_, a), (_, b)| a.partial_cmp(b).unwrap())
                    .expect("mini dp matrix has an empty column");

                let skip_loop_score = matrix[col_idx - 1][0] + (1e55f64).ln() / 30.0;
                let query_to_skip_score = max_score_in_prev_col + (1e55f64).ln() / 2.0;

                if skip_loop_score > query_to_skip_score {
                    matrix[col_idx][0] += skip_loop_score;
                    sources[col_idx][0] = 0;
                } else {
                    matrix[col_idx][0] += query_to_skip_score;
                    sources[col_idx][0] = row_of_max_in_prev_col;
                }

                (1..=alignments.len()).for_each(|row_idx| {
                    let skip_to_query_tuple =
                        (matrix[col_idx - 1][0] + (1e-55f64).ln() / 2.0, 0usize);
                    let loop_tuple = (matrix[col_idx - 1][row_idx], row_idx);
                    let jump_tuple = (
                        max_score_in_prev_col + (1e-55f64).ln(),
                        row_of_max_in_prev_col,
                    );

                    let score_tuples = [skip_to_query_tuple, loop_tuple, jump_tuple];

                    let max_tuple = score_tuples
                        .iter()
                        .max_by(|a, b| a.0.partial_cmp(&b.0).expect("failed to compare floats"))
                        .expect("failed to find max score tuple");

                    matrix[col_idx][row_idx] += max_tuple.0;
                    sources[col_idx][row_idx] = max_tuple.1;
                });
            });

            let last_col = matrix.last().expect("mini dp matrix is empty");

            let (max_row_in_last_col, _) = last_col
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                .unwrap();

            let mut start = matrix.len() - 1;
            let mut prev = max_row_in_last_col;

            let mut trace = vec![];
            (0..target_len - 1).rev().for_each(|col| {
                let source = sources[col + 1][prev];

                if source != prev {
                    trace.push((prev, [col + 1, start]));
                    start = col;
                }
                prev = source;
            });

            trace.push((prev, [0, start]));
            trace.reverse();

            trace
                .iter()
                // filter the skip traces
                .filter(|(i, _)| *i != 0)
                .for_each(|(i, range)| {
                    alignment_ranges.push(AlignmentRange {
                        ali_id: alignments[i - 1].id,
                        assembly_col_start: range[0],
                        assembly_col_end: range[1],
                    });
                });
        };

        Self {
            query_id: alignments[0].query_id,
            strand,
            target_start: assembly_target_start,
            target_end: assembly_target_end,
            query_start,
            query_end,
            alignments,
            alignment_ranges,
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
    pub tandem_repeats: &'a [TandemRepeat],
}

impl<'a> AssemblyGroup<'a> {
    pub fn new(
        group: &ProximityGroup<'a>,
        confidence_avg_by_id: &HashMap<usize, f64>,
        confidence_by_id: &HashMap<usize, Vec<f64>>,
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
                let mut fwd_links = vec![];
                let mut rev_links = vec![];
                if args.assembly_viz {
                    fwd_links = fwd_graph
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

                    rev_links = rev_graph
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
                }

                let mut fwd_assemblies = assembly(&mut fwd_graph, confidence_by_id);

                let mut rev_assemblies = assembly(&mut rev_graph, confidence_by_id);

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

                if args.assembly_viz {
                    if !fwd_ali.is_empty() {
                        let fwd_data = AssemblySodaData::new(
                            &fwd_assemblies,
                            &query_ids,
                            fwd_links,
                            confidence_avg_by_id,
                        );

                        let fwd_path = args.viz_output_path.join(format!("{}-fwd.html", query_id));

                        write_soda_html(
                            &fwd_data,
                            "./fixtures/soda/assembly-template.html",
                            "./fixtures/soda/assembly.js",
                            fwd_path,
                        );
                    }

                    if !rev_ali.is_empty() {
                        let rev_data = AssemblySodaData::new(
                            &rev_assemblies,
                            &query_ids,
                            rev_links,
                            confidence_avg_by_id,
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

                assemblies.append(&mut fwd_assemblies);
                assemblies.append(&mut rev_assemblies);
            });

        assemblies.sort_by_key(|a| a.target_start);

        Self {
            target_id: group.target_id,
            target_start: group.target_start,
            target_end: group.target_end,
            assemblies,
            tandem_repeats: group.tandem_repeats,
        }
    }
}
