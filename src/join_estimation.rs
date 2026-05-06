use crate::{
    segments::Block,
    statistics::{ExponentialEstimator, StudentsT},
};

trait JoinEstimator {
    fn predict(&self, first_block: &Block, second_block: &Block) -> f64;
}

trait JoinStatistics<T: JoinEstimator> {
    fn new() -> Self;
    fn combine(&self, other: &Self) -> Self;
    fn add(&self, first_block: &Block, second_block: &Block, neighbors: bool, joinable: bool);
    fn to_estimator(&self) -> T;
}

struct BayesianJoinEstimator {
    target_distance_join: ExponentialEstimator,
    target_distance_background: ExponentialEstimator,
    divergence_join: StudentsT,
    divergence_background: StudentsT,
}

struct BayesianJoinStatistics {
    joinable_target_distance_sum: usize,
    all_target_distance_sum: usize,
    divergence_sum: f64,
    divergence_square_sum: f64,
    join_divergence_sum: f64,
    join_divergence_square_sum: f64,
    divergence_offset: f64,
    joinable_count: usize,
    all_count: usize,
}
