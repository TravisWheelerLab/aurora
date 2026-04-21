use crate::{
    alignment::AlignmentData,
    pipeline::NaiveTraceResults,
    segments::SegmentView,
    statistics::{Distribution, ExponentialEstimator},
};

#[derive(Debug)]
pub struct RegionStatistics {
    pub total_bases: usize,
    pub unexplained_bases: Vec<usize>,
}

#[derive(Debug, Clone)]
pub struct QueryStatistics<T: Distribution> {
    pub occurances: usize,
    pub coverage: usize,
    pub target_span: usize,
    pub distribution: T,
}

#[derive(Debug)]
pub struct TraceStatistics<T: Distribution> {
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

pub fn trace_statistics(
    naive_traces: &[NaiveTraceResults],
    alignment_data: &AlignmentData,
    count_mode: OccuranceCountingMode,
) -> TraceStatistics<ExponentialEstimator> {
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

    let mut query_stats = vec![
        QueryStatistics {
            occurances: 0,
            coverage: 0,
            target_span: 0,
            distribution: ExponentialEstimator::unit(),
        };
        alignment_data.query_name_map.size()
    ];

    let mut query_span: Vec<Option<(usize, usize)>> =
        vec![None; alignment_data.query_name_map.size()];

    let mut all_region_stats: Vec<RegionStatistics> = Vec::with_capacity(naive_traces.len());

    for trace_results in naive_traces.iter() {
        match count_mode {
            OccuranceCountingMode::Segments => {
                for seg in trace_results.segments.iter_segments() {
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

        let mut region_stat = RegionStatistics {
            total_bases: 0,
            unexplained_bases: Vec::with_capacity(trace_results.segments.len()),
        };

        let mut unexplained_bases_up_to: usize = 0;
        let mut prior_segment: Option<SegmentView> = None;

        for seg in trace_results.segments.iter_segments() {
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

        all_region_stats.push(region_stat);
    }

    for (query_info, query_span) in query_stats.iter_mut().zip(query_span.iter()) {
        if let Some((start, end)) = query_span {
            query_info.target_span = end - start + 1;
            // We subtract 1 because were looking at distances between each occurance as a sample value.
            query_info.distribution = ExponentialEstimator::new(
                query_info.target_span as f64 / query_info.occurances.saturating_sub(1) as f64,
                query_info.occurances.saturating_sub(1),
            );
        }
    }

    TraceStatistics {
        total_bases: all_region_stats.iter().map(|v| v.total_bases).sum(),
        query_statistics: query_stats,
        region_statistics: all_region_stats,
    }
}
