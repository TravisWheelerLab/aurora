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

fn fast_select(a: f64, b: f64, switch: bool) -> f64 {
    /*let sw = (!switch as u64).wrapping_sub(1);
    f64::from_bits((a.to_bits() & sw) | (b.to_bits() & !sw))*/
    if switch {
        a
    } else {
        b
    }
}

const NEG_LN2: f64 = -0.6931471806;

// See https://cran.uib.no/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
fn ln_1m_exp(v: f64) -> f64 {
    let v_neg = -v.abs();
    if v_neg < NEG_LN2 {
        (-v_neg.exp()).ln_1p()
    } else {
        (-v_neg.exp_m1()).ln()
    }
}

impl ScoreParams {
    pub fn new(
        num_alignments: usize,
        query_jump_penalty_nats: f64,
        skip_penalty_multiplier: f64,
    ) -> Self {
        let query_jump_score = query_jump_penalty_nats - (num_alignments.max(1) as f64).ln();
        // jumping to the skip state and then jumping back to a query sequence
        // should be the same cost as jumping between query sequences
        let query_to_skip_score = skip_penalty_multiplier * query_jump_score;

        // staying in the same query sequence is essentially free (this should
        // not be zero mathematically, but it becomes zero due to floating point
        // arithmetic)
        let query_loop_score = ln_1m_exp(query_jump_score);

        // staying in the skip state for some amount of positions is the same as staying in a query sequence...
        let skip_loop_score = query_loop_score; // 0.0;

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

    pub fn join_transition(&self, prior_is_skip: bool, join_probability: f64) -> f64 {
        let true_join_cost = fast_select(
            self.query_to_skip_score + (self.query_loop_score - self.query_jump_score).abs(),
            self.query_loop_score,
            prior_is_skip,
        );
        let false_join_cost = fast_select(
            self.query_to_skip_score,
            self.query_jump_score,
            prior_is_skip,
        );
        true_join_cost * join_probability + (1.0 - join_probability) * false_join_cost
    }
}
