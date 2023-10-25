use std::{fs, path::PathBuf};

use itertools::Itertools;

use crate::{
    alignment::{Alignment, AlignmentData, Strand, VecMap},
    alphabet::UTF8_TO_DIGITAL_NUCLEOTIDE,
    chunks::ProximityGroup,
    confidence::confidence,
    matrix::{Matrix, MatrixDef},
    results::Annotation,
    score_params::ScoreParams,
    segments::Segments,
    support::windowed_confidence_slow,
    viterbi::{traceback, viterbi},
    viz::{write_soda_html, AuroraSodaData},
    windowed_scores::windowed_score,
    Args, BACKGROUND_WINDOW_SIZE, SCORE_WINDOW_SIZE,
};

pub fn run_assembly_pipeline(
    group: &ProximityGroup,
    alignment_data: &AlignmentData,
    region_idx: usize,
    args: &Args,
) {
    let mut viz_path = PathBuf::from(format!(
        "./viz/{}-{}-{}-{}",
        group.target_id, region_idx, group.target_start, group.target_end
    ));

    fs::create_dir_all(&viz_path).unwrap();
    viz_path = viz_path.canonicalize().unwrap();

    let matrix_def = MatrixDef::new(group);
    let mut confidence_matrix = Matrix::<f64>::new(&matrix_def);

    windowed_score(
        &mut confidence_matrix,
        group.alignments,
        &alignment_data.substitution_matrices,
        SCORE_WINDOW_SIZE,
        BACKGROUND_WINDOW_SIZE,
    );

    confidence(&mut confidence_matrix);

    let confidence_avg_by_row = windowed_confidence_slow(&mut confidence_matrix);

    let mut query_ids: Vec<usize> = group
        .alignments
        .iter()
        .map(|a| a.query_id)
        .unique()
        .collect();

    query_ids.sort();

    struct AlignmentTuple<'a> {
        row_idx: usize,
        alignment: &'a Alignment,
        confidence: f64,
    }

    #[derive(PartialEq, Clone, Debug)]
    struct Edge {
        node_idx: usize,
        weight: f64,
    }

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
        .for_each(|(query_id, alignments)| {
            // split the forward and reverse stranded alignments
            let (fwd, rev): (Vec<_>, Vec<_>) = alignments
                .into_iter()
                .partition(|a| a.alignment.strand == Strand::Forward);

            if fwd.is_empty() {
                return;
            }

            fwd.iter().zip(fwd.iter().skip(1)).for_each(|(a, b)| {
                debug_assert!(a.alignment.target_start <= b.alignment.target_start);
            });

            let mut graph: Vec<Vec<Edge>> = vec![vec![]; fwd.len()];
            println!();
            (0..fwd.len()).for_each(|a_idx| {
                (a_idx + 1..fwd.len()).for_each(|b_idx| {
                    let a = &fwd[a_idx];
                    let b = &fwd[b_idx];

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
                            node_idx: b_idx,
                            weight: consensus_distance.abs() as f64,
                        });
                        graph[b_idx].push(Edge {
                            node_idx: a_idx,
                            weight: consensus_distance.abs() as f64,
                        });
                    }
                });
            });

            graph
                .iter()
                .enumerate()
                .for_each(|(i, v)| println!("{i} ({:.4}): {v:?}", fwd[i].confidence));
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
            let mut assemblies: VecMap<Vec<usize>> = VecMap::new();

            let mut assembly: Vec<usize>;

            // while we have alignments left to place
            while !remaining_indices.is_empty() {
                // start an assembly with the highest remaining confidence
                assembly = vec![remaining_indices.pop().unwrap()];
                loop {
                    // the alignment we started with is the
                    // one that's going to maintain edges
                    let root_idx = assembly[0];
                    let root_edges = &graph[root_idx];

                    // find the edge with the minimum weight
                    let join_edge = match root_edges.iter().min_by(|a, b| {
                        a.weight
                            .partial_cmp(&b.weight)
                            .expect("failed to compare edge weights")
                    }) {
                        Some(edge) => edge,
                        // if we don't have an edge,
                        // we're done with this assembly
                        None => break,
                    };

                    let join_idx = join_edge.node_idx;

                    // add this alignment to the assembly
                    assembly.push(join_idx);
                    // and remove it from further consideration

                    // println!("{root_idx} {join_idx}, {remaining_indices:?}");
                    // println!("{edges:?}");
                    // println!();

                    match remaining_indices.iter().position(|&idx| idx == join_idx) {
                        Some(pos) => remaining_indices.remove(pos),
                        None => panic!("failed to remove index from remaining"),
                    };

                    // now take it's edges
                    let join_edges = &graph[join_idx];

                    let root_edge_indices = root_edges.iter().map(|e| e.node_idx).collect_vec();

                    let compatible_indices = join_edges
                        .iter()
                        .filter(|e| root_edge_indices.contains(&e.node_idx))
                        .map(|e| e.node_idx)
                        .collect_vec();

                    let root_compatible_edges = root_edges
                        .iter()
                        .filter(|e| compatible_indices.contains(&e.node_idx))
                        .collect_vec();

                    let join_compatible_edges = join_edges
                        .iter()
                        .filter(|e| compatible_indices.contains(&e.node_idx))
                        .collect_vec();

                    debug_assert_eq!(root_compatible_edges.len(), join_compatible_edges.len());

                    let new_root_edges = root_compatible_edges
                        .iter()
                        .zip(join_compatible_edges)
                        .map(|(r, j)| Edge {
                            node_idx: r.node_idx,
                            weight: r.weight.min(j.weight),
                        })
                        .collect_vec();

                    graph[root_idx] = new_root_edges;

                    if !compatible_indices.is_empty() {
                        println!("{}", compatible_indices.len());
                    }

                    graph[join_idx] = vec![];
                }
            }

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
                .collect::<serde_json::Value>();

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
                .collect::<serde_json::Value>();

            let links = fwd
                .iter()
                .enumerate()
                .flat_map(|(a_idx, a_t)| {
                    graph[a_idx]
                        .iter()
                        .map(|edge| &fwd[edge.node_idx])
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
                .collect::<Vec<_>>();

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

pub fn run_pipeline(
    group: &ProximityGroup,
    alignment_data: &AlignmentData,
    region_idx: usize,
    args: &Args,
) {
    let score_params = ScoreParams::new(
        group.alignments.len(),
        args.query_jump_probability,
        args.num_skip_loops_eq_to_jump,
    );

    let matrix_def = MatrixDef::new(group);
    let mut confidence_matrix = Matrix::<f64>::new(&matrix_def);
    let mut viterbi_matrix = Matrix::<f64>::new(&matrix_def);
    let mut sources_matrix = Matrix::<usize>::new(&matrix_def);

    windowed_score(
        &mut confidence_matrix,
        group.alignments,
        &alignment_data.substitution_matrices,
        SCORE_WINDOW_SIZE,
        BACKGROUND_WINDOW_SIZE,
    );

    confidence(&mut confidence_matrix);

    windowed_confidence_slow(&mut confidence_matrix);

    // this removes any initial dead space between alignments

    // let mut active_cols = confidence_matrix.initial_active_cols();
    let mut active_cols = (0..confidence_matrix.num_cols()).collect::<Vec<_>>();

    let mut last_num_cols = active_cols.len();
    let mut results: Vec<Annotation> = vec![];

    let mut trace_strings = vec![];
    let mut fragment_strings = vec![];
    let mut segment_strings = vec![];

    let mut join_id_cnt = 0usize;

    while !active_cols.is_empty() {
        viterbi(
            &confidence_matrix,
            &mut viterbi_matrix,
            &mut sources_matrix,
            &active_cols,
            &score_params,
            alignment_data.query_name_map.size(),
            args.consensus_join_distance,
        );

        // let data = hyper_traceback(&sources_matrix, group.alignments, alignment_data);
        // write_soda_html(
        //     &data,
        //     "./fixtures/soda/trace-template.html",
        //     "./fixtures/soda/trace.js",
        //     "./trace.html",
        // );

        let trace = traceback(&viterbi_matrix, &sources_matrix, &active_cols, join_id_cnt);
        join_id_cnt = trace
            .iter()
            .map(|t| t.join_id)
            .max()
            .expect("trace is empty");

        let mut segments = Segments::new(&trace, &confidence_matrix, &active_cols);
        segments.create_links(
            args.consensus_join_distance,
            args.target_join_distance,
            args.num_skip_loops_eq_to_jump,
            args.min_fragment_length,
        );
        segments.process_links();

        let (mut iteration_results, cols) = segments.get_results_and_active_cols(
            &alignment_data.query_name_map,
            alignment_data.target_name_map.get(group.target_id),
            group.target_start,
            region_idx,
        );

        results.append(&mut iteration_results);
        active_cols = cols;

        if args.viz {
            trace_strings.push(segments.trace_soda_string(&alignment_data.query_name_map));
            fragment_strings.push(segments.fragments_soda_string(&alignment_data.query_name_map));
            segment_strings.push(segments.segments_soda_string());
        }

        if active_cols.len() == last_num_cols {
            panic!();
        }

        last_num_cols = active_cols.len();
    }

    results.sort_by_key(|r| (r.join_id, r.target_start));
    results.sort_by_key(|r| r.target_start);
    results.retain(|r| r.query_name != "skip");

    Annotation::write(&results, &mut std::io::stdout());

    if args.viz {
        let data = AuroraSodaData::new(
            group,
            &alignment_data.query_name_map,
            &results,
            trace_strings,
            fragment_strings,
            segment_strings,
            args,
        );

        let out_path = args.viz_output_path.join(format!(
            "{}-{}-{}.html",
            alignment_data.target_name_map.get(group.target_id),
            matrix_def.target_start,
            matrix_def.target_start + matrix_def.num_cols
        ));

        write_soda_html(
            &data,
            "./fixtures/soda/template.html",
            "./fixtures/soda/viz.js",
            out_path,
        );
    }
}

pub trait StrSliceExt {
    fn to_digital_nucleotides(self) -> Vec<u8>;
}

impl StrSliceExt for &str {
    fn to_digital_nucleotides(self) -> Vec<u8> {
        self.as_bytes()
            .iter()
            .map(|byte| {
                *UTF8_TO_DIGITAL_NUCLEOTIDE
                    .get(byte)
                    .unwrap_or_else(|| panic!("unknown byte: {byte}"))
            })
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use crate::alignment::{Alignment, Strand};

    use super::{run_pipeline, StrSliceExt};

    #[test]
    fn test_pipeline() {
        let ali = [
            Alignment {
                query_id: 1,
                target_seq: "AAAAAAAAAAGGGGGGGGGG".to_digital_nucleotides(),
                query_seq: "AAAAAAAAAAGGGGGGGGGG".to_digital_nucleotides(),
                target_start: 1,
                target_end: 20,
                query_start: 1,
                query_end: 20,
                strand: Strand::Forward,
                substitution_matrix_id: 0,
            },
            Alignment {
                query_id: 2,
                target_seq: "AAAAAAAAAAGGGGGGGGGG".to_digital_nucleotides(),
                query_seq: "AAAAAAAAAAGGGGGGGGGG".to_digital_nucleotides(),
                target_start: 21,
                target_end: 40,
                query_start: 1,
                query_end: 20,
                strand: Strand::Forward,
                substitution_matrix_id: 0,
            },
        ];
    }
}
