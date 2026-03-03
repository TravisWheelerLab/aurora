use std::fmt::{format, Debug, Display};

use crate::alignment::Strand;
use itertools::Itertools;

#[derive(Clone)]
pub struct SimpleAnnotation {
    pub target_start: usize,
    pub target_end: usize,
    pub query_id: usize,
    pub query_name: String,
    pub query_start: usize,
    pub query_end: usize,
    pub strand: Strand,
    pub kimura80: f64,
}

#[derive(Clone)]
pub struct AmbiguousAnnotation {
    pub target_name: String,
    pub annotations: Vec<SimpleAnnotation>,
    pub confidence: f64,
    pub join_id: usize,
    pub region_id: usize,
}

#[allow(dead_code)]
pub struct ConcreteAnnotation {
    pub target_name: String,
    pub target_start: usize,
    pub target_end: usize,
    pub query_id: usize,
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
    ambiguous_count_width: usize,
    target_name_width: usize,
    target_start_width: usize,
    target_end_width: usize,
    query_name_width: usize,
    query_start_width: usize,
    query_end_width: usize,
    join_id_width: usize,
    kimura_80_width: usize,
    strand_width: usize,
}

fn get_mutli_option_string<B: Display + PartialEq, F>(
    simple_annotations: &[SimpleAnnotation],
    prop: F,
    simplify: bool,
) -> String
where
    F: Fn(&SimpleAnnotation) -> B,
{
    let prop_ref = &prop;

    if !simplify {
        if let Result::Ok(val) = simple_annotations.iter().map(prop_ref).all_equal_value() {
            val.to_string()
        } else {
            simple_annotations.iter().map(prop_ref).join(",")
        }
    } else {
        // Only allow 1 element...
        simple_annotations.iter().take(1).map(prop_ref).join(",")
    }
}

fn get_mutli_option_string_with_format<B: Display + PartialEq, F, G>(
    simple_annotations: &[SimpleAnnotation],
    prop: F,
    simplify: bool,
    formatter: G,
) -> String
where
    F: Fn(&SimpleAnnotation) -> B,
    G: Fn(&B) -> String,
{
    let prop_ref = &prop;
    let get_and_format = |v| formatter(&prop_ref(v));

    if !simplify {
        if let Result::Ok(val) = simple_annotations.iter().map(prop_ref).all_equal_value() {
            formatter(&val)
        } else {
            simple_annotations.iter().map(get_and_format).join(",")
        }
    } else {
        // Only allow 1 element...
        simple_annotations
            .iter()
            .take(1)
            .map(get_and_format)
            .join(",")
    }
}

fn get_strings(simple_annotations: &[SimpleAnnotation], simplify: bool) -> [String; 7] {
    [
        get_mutli_option_string(simple_annotations, |v| v.target_start, simplify),
        get_mutli_option_string(simple_annotations, |v| v.target_end, simplify),
        get_mutli_option_string(simple_annotations, |v| v.query_name.clone(), simplify),
        get_mutli_option_string(simple_annotations, |v| v.query_start, simplify),
        get_mutli_option_string(simple_annotations, |v| v.query_end, simplify),
        get_mutli_option_string(simple_annotations, |v| v.strand, simplify),
        get_mutli_option_string_with_format(
            simple_annotations,
            |v| v.kimura80,
            simplify,
            |v| format!("{:4.3}", v),
        ),
    ]
}

impl AmbiguousAnnotation {
    pub fn line(&self, widths: &LineWidths, simplify: bool) -> String {
        let [ts, te, qn, qs, qe, strand, k80] = get_strings(&self.annotations, simplify);
        format!(
            "{:w0$} {:w1$} {:w2$} {:w3$} {:w4$} {:w5$} {:w6$} {:w7$} {:4.3} {:w8$} {:w9$} {}",
            self.annotations.len(),
            self.target_name,
            ts,
            te,
            qn,
            qs,
            qe,
            strand,
            self.confidence,
            k80,
            self.join_id,
            self.region_id,
            w0 = widths.ambiguous_count_width,
            w1 = widths.target_name_width,
            w2 = widths.target_start_width,
            w3 = widths.target_end_width,
            w4 = widths.query_name_width,
            w5 = widths.query_start_width,
            w6 = widths.query_end_width,
            w7 = widths.strand_width,
            w8 = widths.kimura_80_width,
            w9 = widths.join_id_width,
        )
    }

    pub fn write(
        results: &[AmbiguousAnnotation],
        out: &mut impl std::io::Write,
        simplified: bool,
    ) -> std::io::Result<()> {
        let mut widths = LineWidths::default();

        for result in results {
            let [ts, te, qn, qs, qe, strand, k80] = get_strings(&result.annotations, simplified);

            widths.ambiguous_count_width = widths
                .ambiguous_count_width
                .max(result.annotations.len().to_string().len());
            widths.target_name_width = widths.target_name_width.max(result.target_name.len());
            widths.target_start_width = widths.target_start_width.max(ts.len());
            widths.target_end_width = widths.target_end_width.max(te.len());
            widths.query_name_width = widths.query_name_width.max(qn.len());
            widths.query_start_width = widths.query_start_width.max(qs.len());
            widths.query_end_width = widths.query_end_width.max(qe.len());
            widths.join_id_width = widths.join_id_width.max(result.join_id.to_string().len());
            widths.kimura_80_width = widths.kimura_80_width.max(k80.len());
            widths.strand_width = widths.strand_width.max(strand.len());
        }

        for result in results {
            writeln!(out, "{}", result.line(&widths, simplified))?;
        }

        Ok(())
    }

    pub fn get_target_bounds(&self) -> (usize, usize) {
        (
            self.annotations
                .iter()
                .map(|a| a.target_start)
                .min()
                .unwrap_or(0),
            self.annotations
                .iter()
                .map(|a| a.target_end)
                .max()
                .unwrap_or(0),
        )
    }
}

pub enum AnnotationConversionError {
    NotConcrete,
}

impl TryFrom<&AmbiguousAnnotation> for ConcreteAnnotation {
    type Error = AnnotationConversionError;

    fn try_from(value: &AmbiguousAnnotation) -> Result<Self, Self::Error> {
        if value.annotations.len() != 1 {
            Result::Err(AnnotationConversionError::NotConcrete)
        } else {
            let inner_val = &value.annotations[0];

            Result::Ok(ConcreteAnnotation {
                target_name: value.target_name.clone(),
                target_start: inner_val.target_start,
                target_end: inner_val.target_end,
                query_id: inner_val.query_id,
                query_name: inner_val.query_name.clone(),
                query_start: inner_val.query_start,
                query_end: inner_val.query_end,
                strand: inner_val.strand,
                confidence: value.confidence,
                join_id: value.join_id,
                region_id: value.region_id,
            })
        }
    }
}

impl std::fmt::Display for AmbiguousAnnotation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.line(&LineWidths::default(), false))
    }
}
