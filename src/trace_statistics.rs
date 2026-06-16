use std::fmt::Debug;

use itertools::izip;

use crate::{
    alignment::AlignmentData,
    join_estimation::{JoinEstimator, JoinStatisticsCollector},
    pipeline::NaiveTraceResults,
    segments::{InitialSegments, Segment},
};

#[derive(Debug)]
pub struct RegionStatistics {
    pub total_bases: usize,
    pub unexplained_bases: Vec<usize>,
}

#[derive(Debug, Clone)]
pub struct QueryStatistics<T: JoinEstimator> {
    pub occurances: usize,
    pub coverage: usize,
    pub target_span: usize,
    pub estimator: T,
}

#[derive(Debug)]
pub struct TraceStatistics<T: JoinEstimator> {
    #[allow(dead_code)]
    pub total_bases: usize,
    pub query_statistics: Vec<QueryStatistics<T>>,
    pub region_statistics: Vec<RegionStatistics>,
}

pub enum OccuranceCountingMode {
    Segments,
    #[allow(dead_code)]
    Trace,
}

pub fn calculate_region_statistics(segments: &InitialSegments) -> RegionStatistics {
    let mut region_stat = RegionStatistics {
        total_bases: 0,
        unexplained_bases: Vec::with_capacity(segments.len()),
    };

    let mut unexplained_bases_up_to: usize = 0;
    let mut prior_segment: Option<&Segment> = None;

    for seg in segments.view_segments() {
        if let Some(prior_segment) = prior_segment {
            // If a skip block was the prior block, add it's bases as unexplained.
            if prior_segment.blocks.len() == 1 && prior_segment.blocks[0].row_idx == 0 {
                unexplained_bases_up_to += seg.end_col - seg.start_col + 1;
            }
            unexplained_bases_up_to += seg.start_col - prior_segment.end_col - 1;
            region_stat.total_bases += seg.start_col - prior_segment.end_col - 1;
        }
        region_stat.total_bases += seg.end_col - seg.start_col + 1;
        region_stat.unexplained_bases.push(unexplained_bases_up_to);

        prior_segment = Some(seg);
    }

    region_stat
}

pub fn trace_statistics<S: JoinStatisticsCollector + Debug + Into<E>, E: JoinEstimator>(
    naive_traces: &[NaiveTraceResults<S>],
    alignment_data: &AlignmentData,
    count_mode: OccuranceCountingMode,
) -> TraceStatistics<E> {
    // Asumption... All regions are sorted, no gaps. At least 1 region expected...
    debug_assert!(naive_traces.first().map(|v| v.region_index) == Some(0));
    debug_assert!(naive_traces
        .iter()
        .zip(naive_traces.iter().skip(1))
        .all(|(v1, v2)| v1.region_index + 1 == v2.region_index));

    assert!(naive_traces
        .iter()
        .zip(naive_traces.iter().skip(1))
        .all(|(v1, v2)| v1.region_index + 1 == v2.region_index && v1.target_end < v2.target_start));

    let mut query_stats: Vec<QueryStatistics<E>> = vec![
        QueryStatistics {
            occurances: 0,
            coverage: 0,
            target_span: 0,
            estimator: E::default(),
        };
        alignment_data.query_name_map.size()
    ];

    let mut query_span: Vec<Option<(usize, usize)>> =
        vec![None; alignment_data.query_name_map.size()];

    let mut all_region_stats: Vec<RegionStatistics> = Vec::with_capacity(naive_traces.len());
    // We combine stats for all families to use as a prior (psuedo-count, single sample) for all stats...
    let mut all_family_stats: S = S::new();

    for trace_results in naive_traces.iter() {
        for (_query_id, stats) in trace_results.query_join_statistics.iter() {
            all_family_stats = all_family_stats.combine(stats);
        }

        match count_mode {
            OccuranceCountingMode::Segments => {
                for seg in trace_results.segments.view_segments().iter() {
                    for blk in seg.blocks.iter() {
                        if let Some(query_id) = blk.query_id {
                            query_stats[query_id].occurances += 1;
                            query_stats[query_id].coverage += blk.col_end - blk.col_start + 1;
                            query_span[query_id] = match query_span[query_id] {
                                None => Some((
                                    trace_results.target_start + blk.col_start,
                                    trace_results.target_start + blk.col_end,
                                )),
                                Some((start, end)) => Some((
                                    start.min(trace_results.target_start + blk.col_start),
                                    end.max(trace_results.target_start + blk.col_end),
                                )),
                            }
                        }
                    }
                }
            }
            OccuranceCountingMode::Trace => {
                for trace_blk in trace_results.trace_segments.iter() {
                    query_stats[trace_blk.query_id].occurances += 1;
                    query_stats[trace_blk.query_id].coverage +=
                        trace_blk.col_end - trace_blk.col_start + 1;

                    query_span[trace_blk.query_id] = match query_span[trace_blk.query_id] {
                        None => Some((
                            trace_results.target_start + trace_blk.col_start,
                            trace_results.target_start + trace_blk.col_end,
                        )),
                        Some((start, end)) => Some((
                            start.min(trace_results.target_start + trace_blk.col_start),
                            end.max(trace_results.target_start + trace_blk.col_end),
                        )),
                    }
                }
            }
        }

        all_region_stats.push(calculate_region_statistics(&trace_results.segments));
    }

    // Calculate join statistics for all families using combined prior as a starting point...
    let mut all_join_stats: Vec<S> =
        vec![S::new_from_prior(&all_family_stats, 1); alignment_data.query_name_map.size()];

    for trace_results in naive_traces.iter() {
        for (query_id, stats) in trace_results.query_join_statistics.iter() {
            all_join_stats[*query_id] = all_join_stats[*query_id].combine(stats);
        }
    }

    for (query_info, query_span, join_stat) in
        izip!(query_stats.iter_mut(), query_span.iter(), all_join_stats)
    {
        if let Some((start, end)) = query_span {
            query_info.target_span = end - start + 1;
        }
        query_info.estimator = join_stat.into();
    }

    TraceStatistics {
        total_bases: all_region_stats.iter().map(|v| v.total_bases).sum(),
        query_statistics: query_stats,
        region_statistics: all_region_stats,
    }
}
