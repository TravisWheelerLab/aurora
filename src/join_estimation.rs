use crate::{
    assembly::block_target_distance,
    segments::Block,
    statistics::{Distribution, ExponentialEstimator, HalfT},
};

pub trait JoinEstimator<T: JoinStatistics> {
    fn predict(&self, first_block: &Block, second_block: &Block, log_space: bool) -> f64;
    fn from_statistics(statistics: T) -> Self;
}

pub trait JoinStatistics {
    fn new() -> Self;
    fn combine(&self, other: &Self) -> Self;
    fn add(&mut self, first_block: &Block, second_block: &Block, neighbors: bool, joinable: bool);
}

pub struct BayesianJoinEstimator {
    target_distance_join: ExponentialEstimator,
    target_distance_background: ExponentialEstimator,
    divergence_join: HalfT,
    divergence_background: HalfT,
}

impl JoinEstimator<BayesianJoinStatistics> for BayesianJoinEstimator {
    fn from_statistics(statistics: BayesianJoinStatistics) -> Self {
        let join_td_mean =
            statistics.joinable_target_distance_sum as f64 / statistics.joinable_count as f64;
        let all_td_mean = statistics.all_target_distance_sum as f64 / statistics.all_count as f64;

        // Divergence distributions should have a mean of 0, so we assume that...
        let join_div_std = statistics.join_divergence_square_sum / statistics.joinable_count as f64;
        let all_div_std = statistics.divergence_square_sum / statistics.all_count as f64;

        Self {
            target_distance_join: ExponentialEstimator::new(
                join_td_mean,
                statistics.joinable_count,
            ),
            target_distance_background: ExponentialEstimator::new(
                all_td_mean,
                statistics.all_count,
            ),
            divergence_join: HalfT::new(join_div_std, statistics.joinable_count),
            divergence_background: HalfT::new(all_div_std, statistics.all_count),
        }
    }

    fn predict(&self, first_block: &Block, second_block: &Block, log_space: bool) -> f64 {
        let prior_acc: f64 = 0.95; // Accuracy of the prior estimator of joins...
        let target_dist = block_target_distance(first_block, second_block) as f64;
        // Absolute value as t-dist is symmetric and we want to get prob in tail, also, we know the mean is 0...
        let divergence_diff = (second_block.kimura80 - first_block.kimura80).abs();

        let target_likelihood = self.target_distance_join.logccdf(target_dist)
            - self.target_distance_background.logccdf(target_dist);
        let diverg_likelihood = self.divergence_join.logccdf(divergence_diff)
            - self.divergence_background.logccdf(divergence_diff);

        let score = target_likelihood + diverg_likelihood + prior_acc.ln();

        if log_space {
            score
        } else {
            score.exp()
        }
    }
}

pub struct BayesianJoinStatistics {
    joinable_target_distance_sum: usize,
    all_target_distance_sum: usize,
    divergence_sum: f64,
    divergence_square_sum: f64,
    join_divergence_sum: f64,
    join_divergence_square_sum: f64,
    joinable_count: usize,
    all_count: usize,
}

impl JoinStatistics for BayesianJoinStatistics {
    fn new() -> Self {
        Self {
            joinable_target_distance_sum: 0,
            all_target_distance_sum: 0,
            divergence_sum: 0.0,
            divergence_square_sum: 0.0,
            join_divergence_sum: 0.0,
            join_divergence_square_sum: 0.0,
            joinable_count: 0,
            all_count: 0,
        }
    }

    fn add(&mut self, first_block: &Block, second_block: &Block, neighbors: bool, joinable: bool) {
        let target_dist = block_target_distance(first_block, second_block).abs() as usize;
        let divergence_diff = second_block.kimura80 - first_block.kimura80;

        if joinable {
            self.joinable_target_distance_sum += target_dist;
            self.join_divergence_sum += divergence_diff;
            self.join_divergence_square_sum += divergence_diff * divergence_diff;
            self.joinable_count += 1;
        }

        self.all_target_distance_sum += target_dist;
        self.divergence_sum += divergence_diff;
        self.divergence_square_sum += divergence_diff * divergence_diff;
        self.all_count += 1;
    }

    fn combine(&self, other: &Self) -> Self {
        Self {
            joinable_target_distance_sum: self.joinable_target_distance_sum
                + other.joinable_target_distance_sum,
            all_target_distance_sum: self.all_target_distance_sum + other.all_target_distance_sum,
            divergence_sum: self.divergence_sum + other.divergence_sum,
            divergence_square_sum: self.divergence_square_sum + other.divergence_square_sum,
            join_divergence_sum: self.join_divergence_sum + other.join_divergence_sum,
            join_divergence_square_sum: self.join_divergence_square_sum
                + other.join_divergence_square_sum,
            joinable_count: self.joinable_count + other.joinable_count,
            all_count: self.all_count + other.all_count,
        }
    }
}
