use std::{collections::HashSet, hash::Hash};

use itertools::Itertools;

use crate::{
    alignment::{Alignment, Strand},
    chunks::ProximityGroup,
    score_params::ScoreParams,
    AnnotationArgs,
};

/// The direction of an `Edge` in terms of where
/// `&Alignment` B (value) is in relation to `&Alignment` A (key)
/// in the coordinate space of the chromosome
#[derive(Hash, Eq, PartialEq, Clone, Copy, Debug)]
pub enum Direction {
    Left,
    Right,
}

#[derive(Hash, Eq, PartialEq, Clone, Copy, Debug)]
pub enum LinkType {
    Forward,
    Reverse,
    FRInversion,
    RFInversion,
}

#[derive(Clone, Copy, Debug)]
pub struct Edge {
    pub edge_to: usize,
    pub weight: f64,
    pub direction: Direction,
    pub link_type: LinkType,
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

fn piecewise_linear_cost(
    neg_start: f64,
    pos_start: f64,
    neg_slope: f64,
    pos_slope: f64,
    value: f64,
) -> f64 {
    if value < neg_start {
        (value - neg_start).abs() * neg_slope
    } else if value > pos_slope {
        (value - pos_start).abs() * pos_slope
    } else {
        0.0
    }
}

fn get_link_cost(
    annotation_args: &AnnotationArgs,
    score_params: &ScoreParams,
    consensus_gap: f64,
    target_gap: f64,
) -> f64 {
    // TODO: Replace hard-coded values with cli params....
    let min_value = score_params.query_loop_score;
    let value_range = (score_params.query_loop_score - score_params.query_jump_score).abs();
    let lambda = -value_range * (0.1 / annotation_args.target_join_distance as f64).abs();
    let alpha = -value_range * (1.0 / (annotation_args.consensus_join_overlap - 10) as f64).abs();
    let beta = -value_range * (0.5 / (annotation_args.consensus_join_distance - 10) as f64).abs();

    min_value + piecewise_linear_cost(-10.0, 10.0, alpha, beta, consensus_gap) + lambda * target_gap
}

fn link_assemblies(
    graph: &mut [HashSet<Edge>],
    alignments: &[(usize, &Alignment)],
    score_params: &ScoreParams,
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

                let a_length = a.query_end.abs_diff(a.query_start);
                let b_length = b.query_end.abs_diff(b.query_start);
                let min_length = a_length.min(b_length);

                // Query bounds are reversed for reverse sequences, so the start is actually greater than the end (Ex. start: 1510 -> end: 105)
                let (consensus_distance, link_type) = match (a.strand, b.strand) {
                    (Strand::Forward, Strand::Forward) => (
                        b.query_start as isize - a.query_end as isize,
                        LinkType::Forward,
                    ),
                    (Strand::Reverse, Strand::Reverse) => (
                        a.query_end as isize - b.query_start as isize,
                        LinkType::Reverse,
                    ),
                    (Strand::Forward, Strand::Reverse) => (
                        b.query_end as isize - a.query_end as isize,
                        LinkType::FRInversion,
                    ),
                    (Strand::Reverse, Strand::Forward) => (
                        a.query_end as isize - b.query_end as isize,
                        LinkType::RFInversion,
                    ),
                    _ => panic!("Invalid strand types!"),
                };

                let within_target_distance_threshold = match link_type {
                    LinkType::FRInversion | LinkType::RFInversion => {
                        target_distance.abs() < args.inversion_distance
                    }
                    _ => target_distance < args.target_join_distance as isize,
                };

                let consensus_is_colinear = match link_type {
                    LinkType::FRInversion | LinkType::RFInversion => {
                        consensus_distance.abs() < args.inversion_distance
                    }
                    _ => {
                        consensus_distance > -args.consensus_join_overlap
                            && consensus_distance < args.consensus_join_distance
                    }
                };

                // TODO: Hardcoded, change later...
                let is_significant =
                    min_length >= 10 && -consensus_distance <= ((min_length / 2) as isize);

                let weight = get_link_cost(
                    args,
                    score_params,
                    consensus_distance as f64,
                    target_distance as f64,
                );

                let forward_count: usize = graph[a_idx]
                    .iter()
                    .map(|v| matches!(v.direction, Direction::Right) as usize)
                    .sum();
                let not_reached_forward_limit = forward_count < args.max_forward_links;

                if within_target_distance_threshold
                    && consensus_is_colinear
                    && is_significant
                    && not_reached_forward_limit
                {
                    graph[a_idx].insert(Edge {
                        edge_to: b_idx,
                        weight,
                        direction: Direction::Right,
                        link_type,
                    });
                    graph[b_idx].insert(Edge {
                        edge_to: a_idx,
                        weight,
                        direction: Direction::Left,
                        link_type,
                    });
                }
            });
        });
}

/// Represents graph of compatable alignments on the genome.
/// For each alignment, stores all alignments from the same query in front of it.
pub struct AssemblyGraph {
    pub link_graph: Vec<HashSet<Edge>>,
}

impl AssemblyGraph {
    pub fn new(
        group: &ProximityGroup,
        score_params: &ScoreParams,
        annotation_args: &AnnotationArgs,
    ) -> Self {
        let mut query_ids: Vec<usize> = group
            .alignments
            .iter()
            .map(|a| a.query_id)
            .unique()
            .collect();

        query_ids.sort();

        let mut link_graph: Vec<HashSet<Edge>> = vec![HashSet::new(); group.alignments.len()];

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
            .for_each(|(_query_id, alignments)| {
                link_assemblies(
                    &mut link_graph,
                    &alignments.collect_vec(),
                    score_params,
                    annotation_args,
                );
            });

        println!("{:#?}", link_graph);

        Self { link_graph }
    }
}
