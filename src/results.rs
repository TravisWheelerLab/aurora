use crate::{
    alignment::{self, AlignmentData, Strand},
    viterbi::{TraceSegment, TraceStep2},
};

pub struct Annotation {
    pub target_name: String,
    pub target_start: usize,
    pub target_end: usize,
    pub query_name: String,
    pub query_start: usize,
    pub query_end: usize,
    pub strand: Strand,
    pub confidence: f64,
    pub join_id: usize,
    pub region_id: usize,
}

#[derive(Default)]
pub struct LineWidths {
    target_name_width: usize,
    target_start_width: usize,
    target_end_width: usize,
    query_name_width: usize,
    query_start_width: usize,
    query_end_width: usize,
    join_id_width: usize,
}

impl Annotation {
    pub fn new(segment: &TraceSegment) -> Self {
        // let target_start = start_step.col_idx + group_target_start;
        // let target_end = end_step.col_idx + group_target_start;
        // let target_length = (target_end - target_start + 1) as f64;

        // Self {
        //     target_name: target_name.to_string(),
        //     target_start,
        //     target_end,
        //     query_name: query_name.to_string(),
        //     query_start: start_step.consensus_pos,
        //     query_end: end_step.consensus_pos,
        //     strand: start_step.strand,
        //     confidence: confidence_sum / target_length,
        //     join_id: end_step.join_id,
        //     region_id,
        // }

        Self {
            target_name: todo!(),
            target_start: todo!(),
            target_end: todo!(),
            query_name: todo!(),
            query_start: todo!(),
            query_end: todo!(),
            strand: todo!(),
            confidence: todo!(),
            join_id: todo!(),
            region_id: todo!(),
        }
    }

    pub fn line(&self, widths: &LineWidths) -> String {
        format!(
            "{:w0$} {:w1$} {:w2$} {:w3$} {:w4$} {:w5$} {} {:4.3} {:w6$} {}",
            self.target_name,
            self.target_start,
            self.target_end,
            self.query_name,
            self.query_start,
            self.query_end,
            self.strand,
            self.confidence,
            self.join_id,
            self.region_id,
            w0 = widths.target_name_width,
            w1 = widths.target_start_width,
            w2 = widths.target_end_width,
            w3 = widths.query_name_width,
            w4 = widths.query_start_width,
            w5 = widths.query_end_width,
            w6 = widths.join_id_width,
        )
    }

    pub fn write(results: &Vec<Annotation>, out: &mut impl std::io::Write) {
        let mut widths = LineWidths::default();

        for result in results {
            widths.target_name_width = widths.target_name_width.max(result.target_name.len());
            widths.target_start_width = widths
                .target_start_width
                .max(result.target_start.to_string().len());
            widths.target_end_width = widths
                .target_end_width
                .max(result.target_end.to_string().len());
            widths.query_name_width = widths.query_name_width.max(result.query_name.len());
            widths.query_start_width = widths
                .query_start_width
                .max(result.query_start.to_string().len());
            widths.query_end_width = widths
                .query_end_width
                .max(result.query_end.to_string().len());
            widths.join_id_width = widths.join_id_width.max(result.join_id.to_string().len())
        }

        for result in results {
            writeln!(out, "{}", result.line(&widths)).expect("failed to write result line");
        }
    }
}

impl std::fmt::Display for Annotation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.line(&LineWidths::default()))
    }
}
