use std::{fmt::Debug, ops};

use crate::{
    assembly::{block_consensus_distance, block_length_on_query, block_target_distance, LinkType},
    segments::Block,
    statistics::{ln_add_exp, Distribution, ExponentialEstimator, Frechet, HalfT, Laplace},
};

pub trait JoinEstimator: Clone + Default + Debug {
    fn predict(&self, first_block: &Block, second_block: &Block, log_space: bool) -> f64;
}

pub trait JoinStatisticsCollector: Clone + Debug {
    fn new() -> Self;
    fn new_from_prior(bayesian_prior: &Self, pseudo_count: usize) -> Self;
    fn combine(&self, other: &Self) -> Self;
    fn add(&mut self, first_block: &Block, second_block: &Block, neighbors: bool, joinable: bool);
}

#[derive(Debug, Clone, Default)]
pub struct BayesianJoinEstimator {
    target_distance_join: ExponentialEstimator,
    target_distance_nojoin: ExponentialEstimator,
    divergence_join: HalfT,
    divergence_nojoin: HalfT,
    consensus_distance_join: Frechet,
    consensus_distance_nojoin: Laplace,
    join_prior: f64,
}

impl JoinEstimator for BayesianJoinEstimator {
    fn predict(&self, first_block: &Block, second_block: &Block, log_space: bool) -> f64 {
        let target_dist = block_target_distance(first_block, second_block) as f64;
        // Absolute value as t-dist is symmetric and we want to get prob in tail, also, we know the mean is 0...
        let divergence_diff = (second_block.kimura80 - first_block.kimura80).abs();
        let (consensus_dist, _join_type) = block_consensus_distance(first_block, second_block);
        let rel_con_dist = consensus_dist as f64
            / block_length_on_query(first_block).max(block_length_on_query(second_block)) as f64;

        /*
        println!("{:#?}", self);
        println!(
            "{} {} {}",
            rel_con_dist,
            self.consensus_distance_join.pdf(rel_con_dist),
            self.consensus_distance_nojoin.pdf(rel_con_dist)
        );
        println!(
            "{} {} {}",
            target_dist,
            self.target_distance_join.pdf(target_dist),
            self.target_distance_nojoin.pdf(target_dist)
        );
        println!(
            "{} {} {}",
            divergence_diff,
            self.divergence_join.pdf(divergence_diff),
            self.divergence_nojoin.pdf(divergence_diff)
        );*/

        let join_score = self.join_prior.ln()
            + self.target_distance_join.logpdf(target_dist)
            + self.divergence_join.logpdf(divergence_diff);
        //+ self.consensus_distance_join.logpdf(rel_con_dist);
        let nojoin_score = (-self.join_prior).ln_1p()
            + self.target_distance_nojoin.logpdf(target_dist)
            + self.divergence_nojoin.logpdf(divergence_diff);
        //+ self.consensus_distance_nojoin.logpdf(rel_con_dist);

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

#[derive(Debug, Clone, Copy)]
struct MomentEstimator {
    sum_square: f64,
    sum: f64,
    samples: usize,
}

impl MomentEstimator {
    fn new() -> Self {
        Self {
            sum_square: 0.0,
            sum: 0.0,
            samples: 0,
        }
    }

    fn to_psuedo_count(&self, count: usize) -> Self {
        Self {
            sum_square: (self.sum_square / self.samples.max(1) as f64) * count as f64,
            sum: (self.sum / self.samples.max(1) as f64) * count as f64,
            samples: count,
        }
    }

    fn mean(&self) -> f64 {
        self.sum / self.samples.max(1) as f64
    }

    fn variance(&self) -> f64 {
        // TODO: Use shifted data alg for more accuracy...
        (self.sum_square - (self.sum * self.sum) / self.samples.max(1) as f64)
            / (self.samples.max(2) as f64 - 1.0)
    }

    fn standard_deviation(&self) -> f64 {
        self.variance().sqrt()
    }

    fn samples(&self) -> usize {
        self.samples
    }
}

impl Default for MomentEstimator {
    fn default() -> Self {
        Self::new()
    }
}

impl ops::Add<MomentEstimator> for MomentEstimator {
    type Output = MomentEstimator;

    fn add(self, rhs: MomentEstimator) -> Self::Output {
        Self {
            sum_square: self.sum_square + rhs.sum_square,
            sum: self.sum + rhs.sum,
            samples: self.samples + rhs.samples,
        }
    }
}

impl ops::AddAssign<MomentEstimator> for MomentEstimator {
    fn add_assign(&mut self, rhs: MomentEstimator) {
        self.sum_square += rhs.sum_square;
        self.sum += rhs.sum;
        self.samples += rhs.samples;
    }
}

impl ops::AddAssign<f64> for MomentEstimator {
    fn add_assign(&mut self, rhs: f64) {
        self.sum_square += rhs * rhs;
        self.sum += rhs;
        self.samples += 1;
    }
}

impl From<MomentEstimator> for ExponentialEstimator {
    fn from(value: MomentEstimator) -> Self {
        Self::new(value.mean(), value.samples().max(1))
    }
}

impl From<MomentEstimator> for HalfT {
    fn from(value: MomentEstimator) -> Self {
        Self::from_sample_mean(value.mean(), value.samples().max(1))
    }
}

impl From<MomentEstimator> for Laplace {
    fn from(value: MomentEstimator) -> Self {
        Self::from_moments(value.mean(), value.standard_deviation())
    }
}

impl From<&BayesianJoinStatistics> for BayesianJoinEstimator {
    fn from(statistics: &BayesianJoinStatistics) -> Self {
        Self {
            target_distance_join: statistics.joinable_target_distance.into(),
            target_distance_nojoin: statistics.unjoinable_target_distance.into(),
            divergence_join: statistics.joinable_divergence.into(),
            divergence_nojoin: statistics.unjoinable_divergence.into(),
            consensus_distance_join: Frechet::from_log_moments(
                statistics.joinable_consensus_log.mean(),
                statistics.joinable_consensus_log.standard_deviation(),
                -0.05,
            ),
            consensus_distance_nojoin: statistics.unjoinable_consensus.into(),
            // We take sqrt since we count all pairs, not just neighbors.
            join_prior: (statistics.joinable_target_distance.samples() as f64
                / (statistics.joinable_target_distance.samples()
                    + statistics.unjoinable_target_distance.samples())
                .max(1) as f64)
                .sqrt(),
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct BayesianJoinStatistics {
    joinable_target_distance: MomentEstimator,
    unjoinable_target_distance: MomentEstimator,
    joinable_divergence: MomentEstimator,
    unjoinable_divergence: MomentEstimator,
    joinable_consensus_log: MomentEstimator,
    unjoinable_consensus: MomentEstimator,
}

impl JoinStatisticsCollector for BayesianJoinStatistics {
    fn new() -> Self {
        Self::default()
    }

    fn new_from_prior(bayesian_prior: &Self, pseudo_count: usize) -> Self {
        Self {
            joinable_target_distance: bayesian_prior
                .joinable_target_distance
                .to_psuedo_count(pseudo_count),
            unjoinable_target_distance: bayesian_prior
                .unjoinable_target_distance
                .to_psuedo_count(pseudo_count),
            joinable_divergence: bayesian_prior
                .joinable_divergence
                .to_psuedo_count(pseudo_count),
            unjoinable_divergence: bayesian_prior
                .unjoinable_divergence
                .to_psuedo_count(pseudo_count),
            joinable_consensus_log: bayesian_prior
                .joinable_consensus_log
                .to_psuedo_count(pseudo_count),
            unjoinable_consensus: bayesian_prior
                .unjoinable_consensus
                .to_psuedo_count(pseudo_count),
        }
    }

    fn add(&mut self, first_block: &Block, second_block: &Block, _neighbors: bool, joinable: bool) {
        let target_dist = block_target_distance(first_block, second_block).abs() as usize;
        let divergence_diff = (second_block.kimura80 - first_block.kimura80).abs();
        let (consensus_dist, join_type) = block_consensus_distance(first_block, second_block);
        let rel_con_dist = consensus_dist as f64
            / block_length_on_query(first_block).max(block_length_on_query(second_block)) as f64;

        if joinable {
            self.joinable_target_distance += target_dist as f64;
            self.joinable_divergence += divergence_diff;
            if matches!(join_type, LinkType::Forward | LinkType::Reverse) {
                //println!("CDist: {}", rel_con_dist);
                self.joinable_consensus_log += (rel_con_dist + 0.05).max(1e-50).ln();
            }
        } else {
            self.unjoinable_target_distance += target_dist as f64;
            self.unjoinable_divergence += divergence_diff;
            if matches!(join_type, LinkType::Forward | LinkType::Reverse) {
                self.unjoinable_consensus += rel_con_dist;
            }
        }
    }

    fn combine(&self, other: &Self) -> Self {
        Self {
            joinable_target_distance: self.joinable_target_distance
                + other.joinable_target_distance,
            unjoinable_target_distance: self.unjoinable_target_distance
                + other.unjoinable_target_distance,
            joinable_divergence: self.joinable_divergence + other.joinable_divergence,
            unjoinable_divergence: self.unjoinable_divergence + other.unjoinable_divergence,
            joinable_consensus_log: self.joinable_consensus_log + other.joinable_consensus_log,
            unjoinable_consensus: self.unjoinable_consensus + other.unjoinable_consensus,
        }
    }
}
