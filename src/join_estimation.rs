use std::fmt::Debug;

use crate::{
    assembly::block_target_distance,
    segments::Block,
    statistics::{ln_add_exp, Distribution, ExponentialEstimator, HalfT},
};

pub trait JoinEstimator: Clone + Default + Debug {
    fn predict(&self, first_block: &Block, second_block: &Block, log_space: bool) -> f64;
}

pub trait JoinStatisticsCollector: Clone + Debug {
    fn new() -> Self;
    fn new_from_prior(bayesian_prior: &Self) -> Self;
    fn combine(&self, other: &Self) -> Self;
    fn add(&mut self, first_block: &Block, second_block: &Block, neighbors: bool, joinable: bool);
}

#[derive(Debug, Clone, Default)]
pub struct BayesianJoinEstimator {
    target_distance_join: ExponentialEstimator,
    target_distance_nojoin: ExponentialEstimator,
    divergence_join: HalfT,
    divergence_nojoin: HalfT,
    join_prior: f64,
}

impl JoinEstimator for BayesianJoinEstimator {
    fn predict(&self, first_block: &Block, second_block: &Block, log_space: bool) -> f64 {
        let target_dist = block_target_distance(first_block, second_block) as f64;
        // Absolute value as t-dist is symmetric and we want to get prob in tail, also, we know the mean is 0...
        let divergence_diff = (second_block.kimura80 - first_block.kimura80).abs();

        let join_score = self.join_prior.ln()
            + self.target_distance_join.logpdf(target_dist)
            + self.divergence_join.logpdf(divergence_diff);
        let nojoin_score = (-self.join_prior).ln_1p()
            + self.target_distance_nojoin.logpdf(target_dist)
            + self.divergence_nojoin.logpdf(divergence_diff);

        let score_norm = ln_add_exp(join_score, nojoin_score);
        let score = join_score - score_norm;

        if log_space {
            score
        } else {
            score.exp()
        }
    }
}

impl From<BayesianJoinStatistics> for BayesianJoinEstimator {
    fn from(value: BayesianJoinStatistics) -> Self {
        Self::from(&value)
    }
}

impl From<&BayesianJoinStatistics> for BayesianJoinEstimator {
    fn from(statistics: &BayesianJoinStatistics) -> Self {
        let join_psuedo_count = statistics.joinable_count.max(1);
        let nojoin_psuedo_count = statistics.unjoinable_count.max(1);

        let join_td_mean =
            (statistics.joinable_target_distance_sum as f64 / join_psuedo_count as f64).max(1.0);
        let nojoin_td_mean = (statistics.unjoinable_target_distance_sum as f64
            / nojoin_psuedo_count as f64)
            .max(join_td_mean);

        // Divergence distributions should have a mean of 0, so we assume that...
        let join_div_mean =
            (statistics.joinable_divergence_sum / join_psuedo_count as f64).max(1.0);
        let nojoin_div_mean =
            (statistics.unjoinable_divergence_sum / nojoin_psuedo_count as f64).max(join_div_mean);

        Self {
            target_distance_join: ExponentialEstimator::new(join_td_mean, join_psuedo_count),
            target_distance_nojoin: ExponentialEstimator::new(nojoin_td_mean, nojoin_psuedo_count),
            divergence_join: HalfT::from_sample_mean(join_div_mean, join_psuedo_count),
            divergence_nojoin: HalfT::from_sample_mean(nojoin_div_mean, nojoin_psuedo_count),
            // We take sqrt since we count all pairs, not just neighbors.
            join_prior: (join_psuedo_count as f64
                / (nojoin_psuedo_count + join_psuedo_count) as f64)
                .sqrt(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct BayesianJoinStatistics {
    joinable_target_distance_sum: usize,
    unjoinable_target_distance_sum: usize,
    joinable_divergence_sum: f64,
    unjoinable_divergence_sum: f64,
    joinable_count: usize,
    unjoinable_count: usize,
}

impl JoinStatisticsCollector for BayesianJoinStatistics {
    fn new() -> Self {
        Self {
            joinable_target_distance_sum: 0,
            unjoinable_target_distance_sum: 0,
            joinable_divergence_sum: 0.0,
            unjoinable_divergence_sum: 0.0,
            joinable_count: 0,
            unjoinable_count: 0,
        }
    }

    fn new_from_prior(bayesian_prior: &Self) -> Self {
        let join_psuedo_count = bayesian_prior.joinable_count.max(1);
        let nojoin_psuedo_count = bayesian_prior.unjoinable_count.max(1);

        Self {
            joinable_target_distance_sum: bayesian_prior.joinable_target_distance_sum
                / join_psuedo_count,
            unjoinable_target_distance_sum: bayesian_prior.unjoinable_target_distance_sum
                / nojoin_psuedo_count,
            joinable_divergence_sum: bayesian_prior.joinable_divergence_sum
                / join_psuedo_count as f64,
            unjoinable_divergence_sum: bayesian_prior.unjoinable_divergence_sum
                / nojoin_psuedo_count as f64,
            joinable_count: 1,
            unjoinable_count: 1,
        }
    }

    fn add(&mut self, first_block: &Block, second_block: &Block, _neighbors: bool, joinable: bool) {
        let target_dist = block_target_distance(first_block, second_block).abs() as usize;
        let divergence_diff = (second_block.kimura80 - first_block.kimura80).abs();

        if joinable {
            self.joinable_target_distance_sum += target_dist;
            self.joinable_divergence_sum += divergence_diff;
            self.joinable_count += 1;
        } else {
            self.unjoinable_target_distance_sum += target_dist;
            self.unjoinable_divergence_sum += divergence_diff;
            self.unjoinable_count += 1;
        }
    }

    fn combine(&self, other: &Self) -> Self {
        Self {
            joinable_target_distance_sum: self.joinable_target_distance_sum
                + other.joinable_target_distance_sum,
            unjoinable_target_distance_sum: self.unjoinable_target_distance_sum
                + other.unjoinable_target_distance_sum,
            joinable_divergence_sum: self.joinable_divergence_sum + other.joinable_divergence_sum,
            unjoinable_divergence_sum: self.unjoinable_divergence_sum
                + other.unjoinable_divergence_sum,
            joinable_count: self.joinable_count + other.joinable_count,
            unjoinable_count: self.unjoinable_count + other.unjoinable_count,
        }
    }
}
