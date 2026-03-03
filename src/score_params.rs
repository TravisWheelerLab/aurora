#[derive(Debug)]
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

pub fn approximate_ideal_skip_state_score(
    num_skip_loops_match_jump: f64,
    query_jump_penalty_nats: f64,
    penalty_shift: f64,
) -> f64 {
    -(query_jump_penalty_nats / num_skip_loops_match_jump) + penalty_shift
}

fn fast_select(a: f64, b: f64, switch: bool) -> f64 {
    /*let sw = (!switch as u64).wrapping_sub(1);
    f64::from_bits((a.to_bits() & sw) | (b.to_bits() & !sw))*/
    if switch {
        a
    } else {
        b
    }
}

impl ScoreParams {
    pub fn new(
        num_alignments: usize,
        query_jump_penalty_nats: f64,
        num_skip_loops_eq_to_jump: usize,
    ) -> Self {
        let query_jump_score = query_jump_penalty_nats - (num_alignments as f64).ln();
        // jumping to the skip state and then jumping back to a query sequence
        // should be the same cost as jumping between query sequences
        // note for 2nd term: add a tiny bit of penalty to make entering and leaving asymmetrical
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

    /// Compute the log-probability (base e) of transitioning between 2 states, given the following:
    ///  - Is one of the two states a skip state?
    ///  - Is the prior state a different row than the current state (different alignment).
    pub fn transition(&self, is_skip: bool, prior_is_different: bool) -> f64 {
        fast_select(
            fast_select(self.query_to_skip_score, self.query_jump_score, is_skip),
            fast_select(self.skip_loop_score, self.query_loop_score, is_skip),
            prior_is_different,
        )
    }
}
