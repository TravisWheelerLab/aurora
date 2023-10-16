pub struct ScoreParams {
    /// T_m from the paper
    pub query_jump_score: f64,
    /// T_ms from the paper
    pub query_to_skip_score: f64,
    /// T_s from the paper
    pub query_loop_score: f64,
    /// T_ss from the paper
    pub skip_loop_score: f64,
}

impl ScoreParams {
    pub fn new(
        num_alignments: usize,
        query_jump_probability: f64,
        num_skip_loops_eq_to_jump: usize,
    ) -> Self {
        let query_jump_probability_per = query_jump_probability / num_alignments as f64;

        let query_jump_score = query_jump_probability_per.ln();

        // jumping to the skip state and then jumping back to a query sequence
        // should be the same cost as jumping between query sequences
        let query_to_skip_score =
            (query_jump_score / 2.0) + query_jump_score / num_skip_loops_eq_to_jump as f64;

        // staying in the same query sequence is essentially free (this should
        // not be zero mathematically, but it becomes zero due to floating point
        // arithmetic)
        // TODO: if we never want to parameterize this, we should remove it entirely
        let query_loop_score = 0.0;

        // staying in the skip state for some amount of positions (default = 30)
        // should be the same cost as jumping between query sequences
        let skip_loop_score = query_jump_score / num_skip_loops_eq_to_jump as f64;

        Self {
            query_jump_score,
            query_to_skip_score,
            query_loop_score,
            skip_loop_score,
        }
    }
}
