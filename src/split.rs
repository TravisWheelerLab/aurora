use crate::viterbi::TraceSegment;

#[derive(Debug)]
pub struct MatrixRange {
    pub col_start: usize,
    pub col_end: usize,
}

///
///
///
///
pub struct SplitResults {
    pub trace_ambiguous: Vec<TraceSegment>,
    pub trace_conclusive: Vec<TraceSegment>,
    pub resolved_assembly_rows: Vec<usize>,
    pub unresolved_assembly_rows: Vec<usize>,
    pub competed_assembly_rows: Vec<usize>,
    pub inactive_col_ranges: Vec<MatrixRange>,
}
