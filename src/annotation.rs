use itertools::Itertools;

use crate::alignment::Strand;

#[derive(Clone)]
pub struct SimpleAnnotation {
    pub target_start: usize,
    pub target_end: usize,
    pub query_id: usize,
    pub query_name: String,
    pub query_start: usize,
    pub query_end: usize,
    pub strand: Strand,
}

#[derive(Clone)]
pub struct AmbiguousAnnotation {
    pub target_name: String,
    pub annotations: Vec<SimpleAnnotation>,
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

fn get_strings(simple_annotations: &[SimpleAnnotation]) -> [String; 6] {
    [
        simple_annotations.iter().map(|v| v.target_start).join(","),
        simple_annotations.iter().map(|v| v.target_end).join(","),
        simple_annotations
            .iter()
            .map(|v| v.query_name.clone())
            .join(","),
        simple_annotations.iter().map(|v| v.query_start).join(","),
        simple_annotations.iter().map(|v| v.query_end).join(","),
        simple_annotations.iter().map(|v| v.strand).join(","),
    ]
}

impl AmbiguousAnnotation {
    pub fn line(&self, widths: &LineWidths) -> String {
        let [ts, te, qn, qs, qe, strand] = get_strings(&self.annotations);

        format!(
            "{:w0$} {:w1$} {:w2$} {:w3$} {:w4$} {:w5$} {} {:4.3} {:w6$} {}",
            self.target_name,
            ts,
            te,
            qn,
            qs,
            qe,
            strand,
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

    pub fn write(results: &Vec<AmbiguousAnnotation>, out: &mut impl std::io::Write) {
        let mut widths = LineWidths::default();

        for result in results {
            let [ts, te, qn, qs, qe, ..] = get_strings(&result.annotations);

            widths.target_name_width = widths.target_name_width.max(result.target_name.len());
            widths.target_start_width = widths.target_start_width.max(ts.len());
            widths.target_end_width = widths.target_end_width.max(te.len());
            widths.query_name_width = widths.query_name_width.max(qn.len());
            widths.query_start_width = widths.query_start_width.max(qs.to_string().len());
            widths.query_end_width = widths.query_end_width.max(qe.to_string().len());
            widths.join_id_width = widths.join_id_width.max(result.join_id.to_string().len())
        }

        for result in results {
            writeln!(out, "{}", result.line(&widths)).expect("failed to write result line");
        }
    }
}

impl std::fmt::Display for AmbiguousAnnotation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.line(&LineWidths::default()))
    }
}
