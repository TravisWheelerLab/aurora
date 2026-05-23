use std::cmp::Ordering;

use crate::segments::MergeIterator;
use itertools::izip;
// Implementation of P2 estimator.
// See "The P2 Algorithm for Dynamic Statistical Computing Calculation of Quantiles and Histograms Without Storing Observations"
// at https://www.cse.wustl.edu/~jain/papers/ftp/psqr.pdf
//
// We replace the P2 interpolation with PCHIP instead (See paper A Method for Constructing Local Monotone Piecewise Cubic Interpolants by F. N. Fritsch and J. Butland, or https://doi.org/10.1137/0905021)

struct QuantileEstimatorData<'a> {
    ranks: &'a [usize],
    values: &'a [f64],
    targets: &'a [f64],
    observations: &'a usize,
}

struct MutableQuantileEstimatorData<'a> {
    ranks: &'a mut [usize],
    values: &'a mut [f64],
    targets: &'a [f64],
    observations: &'a mut usize,
}

fn cubic_hermite_spline(x0: f64, y0: f64, x1: f64, y1: f64, m0: f64, m1: f64, x: f64) -> f64 {
    let t = (x - x0) / (x1 - x0);
    let ms0 = (x1 - x0) * m0;
    let ms1 = (x1 - x0) * m1;

    let h0 = (2.0 * t - 3.0) * t * t + 1.0;
    let h1 = ((t - 1.0) * t + 1.0) * t;
    let h2 = (2.0 * t + 3.0) * t * t;
    let h3 = (t - 1.0) * t * t;

    h0 * y0 + h1 * ms0 + h2 * y1 + h3 * ms1
}

fn pchip_point_derivative(dx0: f64, dy0: f64, dx1: f64, dy1: f64) -> f64 {
    if dy0 * dy1 > 0.0 {
        let alpha = (1.0 / 3.0) * (1.0 + dx1 / (dx0 + dx1));
        dy0 * dy1 / (alpha * dy1 + (1.0 - alpha) * dy0)
    } else {
        0.0
    }
}

fn pchip_prediction(ranks: &[f64; 4], values: &[f64; 4], x: f64) -> f64 {
    let m0 = pchip_point_derivative(
        ranks[1] - ranks[0],
        values[1] - values[0],
        ranks[2] - ranks[1],
        values[2] - values[1],
    );
    let m1 = pchip_point_derivative(
        ranks[2] - ranks[1],
        values[2] - values[1],
        ranks[3] - ranks[2],
        values[3] - values[2],
    );

    cubic_hermite_spline(ranks[1], values[1], ranks[2], values[2], m0, m1, x)
}

fn debug_check_valid_estimator(
    ranks: &[usize],
    values: &[f64],
    targets: &[f64],
    observations: usize,
) {
    debug_assert!(ranks.len() > 2);
    debug_assert!(ranks.len() == values.len() && ranks.len() == targets.len());
    debug_assert!(values.is_sorted() && ranks.is_sorted() && targets.is_sorted());
    debug_assert!(targets.first() == Some(&0.0) && targets.last() == Some(&1.0));
    debug_assert!(observations >= ranks.len());
    debug_assert!(ranks.first() == Some(&0) && ranks.last() == Some(&observations));
}

fn debug_check_uninitialized_estimator(ranks: &[usize], values: &[f64], targets: &[f64]) {
    debug_assert!(ranks.len() > 2);
    debug_assert!(ranks.len() == values.len() && ranks.len() == targets.len());
    debug_assert!(targets.is_sorted());
    debug_assert!(targets.first() == Some(&0.0) && targets.last() == Some(&1.0));
}

fn _add_sample_to_estimator(data: MutableQuantileEstimatorData, sample: f64) {
    let MutableQuantileEstimatorData {
        ranks,
        values,
        targets,
        observations,
    } = data;
    debug_check_valid_estimator(ranks, values, targets, *observations);
    // Find where sample falls within distribution...
    let p = values.partition_point(|&v| v <= sample);
    let bound_p = p.min(values.len() - 1);

    // Update extremes...
    if bound_p == 0 {
        values[bound_p] = values[bound_p].min(sample);
    } else if bound_p == (values.len() - 1) {
        values[bound_p] = values[bound_p].max(sample);
    }

    // Increment ranks of markers above newly inserted sample...
    for i in (bound_p + 1)..ranks.len() {
        ranks[i] = ranks[1] + 1;
    }

    // Adjust inner markers to within 1 of their target quantile using p2 formula...
    for i in 1..(values.len() - 1) {
        let target_rank = (targets[i] * (*observations) as f64) as usize;
        let true_rank = ranks[i];
        if true_rank.abs_diff(target_rank) > 1 {
            let new_rank: usize = target_rank.clamp(
                ranks[i - 1].saturating_add(1),
                ranks[i + 1].saturating_sub(1),
            );
            if new_rank == true_rank {
                continue;
            }

            let idx_shift = if new_rank > target_rank { 1 } else { 0 };
            let indexes = [
                (i + idx_shift).saturating_sub(2),
                (i + idx_shift).saturating_sub(1),
                (i + idx_shift),
                (i + idx_shift).saturating_add(1).min(ranks.len() - 1),
            ];

            ranks[i] = new_rank;
            values[i] = pchip_prediction(
                &indexes.map(|i| ranks[i] as f64),
                &indexes.map(|i| values[i]),
                new_rank as f64,
            );
        }
    }

    *observations += 1;
}

fn _merge_estimators(
    q1: QuantileEstimatorData,
    q2: QuantileEstimatorData,
    mut new_estimator: MutableQuantileEstimatorData,
) {
    debug_check_valid_estimator(q1.ranks, q1.values, q1.targets, *q1.observations);
    debug_check_valid_estimator(q2.ranks, q2.values, q2.targets, *q2.observations);
    debug_check_uninitialized_estimator(
        new_estimator.ranks,
        new_estimator.values,
        new_estimator.targets,
    );

    assert!(new_estimator.targets.len() <= (*q1.observations + *q2.observations));

    fn get_at(a: &QuantileEstimatorData, i: usize) -> (usize, f64) {
        (a.ranks[i], a.values[i])
    }

    fn set_at(a: &mut MutableQuantileEstimatorData, i: usize, data: (usize, f64)) {
        a.ranks[i] = data.0;
        a.values[i] = data.1;
    }

    // Initialize the min/max quantiles...
    if q1.values[0] <= q2.values[0] {
        set_at(&mut new_estimator, 0, get_at(&q1, 0));
    } else {
        set_at(&mut new_estimator, 0, get_at(&q2, 0));
    }

    let new_est_len = new_estimator.ranks.len();
    if q1.values[q1.values.len() - 1] >= q2.values[q2.values.len() - 1] {
        set_at(
            &mut new_estimator,
            new_est_len - 1,
            get_at(&q1, q1.ranks.len() - 1),
        );
    } else {
        set_at(
            &mut new_estimator,
            new_est_len - 1,
            get_at(&q2, q2.ranks.len() - 1),
        );
    }

    // May eventually replace with algorithm that doesn't use extra memory...
    // Calculate a "merged" quantiles by linearly iterpolating ranks based on the values we see...
    let mut dual_est_quants: Vec<(f64, f64)> = Vec::with_capacity(q1.ranks.len() + q2.ranks.len());

    let mut q1_prior: Option<(usize, f64)> = None;
    let mut q2_prior: Option<(usize, f64)> = None;

    let mut q1_idx = 0;
    let mut q2_idx = 0;

    loop {
        let q1_past_end = q1_idx >= q1.values.len();
        let q2_past_end = q2_idx >= q2.values.len();

        if q1_past_end && q2_past_end {
            break;
        } else if q1_past_end {
            let next = get_at(&q2, q2_idx);
            dual_est_quants.push(((next.0 + q1_prior.map(|v| v.0).unwrap_or(0)) as f64, next.1));
            q2_prior = Some(next);
            q2_idx += 1;
        } else if q2_past_end {
            let next = get_at(&q1, q1_idx);
            dual_est_quants.push(((next.0 + q2_prior.map(|v| v.0).unwrap_or(0)) as f64, next.1));
            q1_prior = Some(next);
            q1_idx += 1;
        } else if q1.values[q1_idx] <= q2.values[q2_idx] {
            let other_next = get_at(&q2, q2_idx);
            let next = get_at(&q1, q1_idx);
            let w = q2_prior
                .map(|other_prior| (next.1 - other_prior.1) / (other_next.1 - other_prior.1))
                .unwrap_or(0.0);
            let other_rank_est = q2_prior
                .map(|other_prior| other_prior.0 as f64 * (1.0 - w) + other_next.0 as f64 * w)
                .unwrap_or(0.0);
            dual_est_quants.push((next.0 as f64 + other_rank_est, next.1));
            q1_prior = Some(next);
            q1_idx += 1;
        } else {
            let other_next = get_at(&q1, q1_idx);
            let next = get_at(&q2, q2_idx);
            let w = q1_prior
                .map(|other_prior| (next.1 - other_prior.1) / (other_next.1 - other_prior.1))
                .unwrap_or(0.0);
            let other_rank_est = q1_prior
                .map(|other_prior| other_prior.0 as f64 * (1.0 - w) + other_next.0 as f64 * w)
                .unwrap_or(0.0);
            dual_est_quants.push((next.0 as f64 + other_rank_est, next.1));
            q2_prior = Some(next);
            q2_idx += 1;
        }
    }

    // New number of observations is the sum of both...
    *(new_estimator.observations) = *(q1.observations) + *(q2.observations);

    // Solve all inner quantiles using traditional interpolation...
    let mut index_between = 0;

    for ti in 1..new_estimator.targets.len() - 1 {
        // Calculate new rank...
        let target = new_estimator.targets[ti];
        let approx_obs_rank = ((target * *(new_estimator.observations) as f64) as usize).clamp(
            1 + ti,
            *(new_estimator.observations) - (new_estimator.targets.len() - (ti + 1)),
        );

        // Find where it lands in cdf...
        while index_between < dual_est_quants.len()
            && (approx_obs_rank as f64) < dual_est_quants[index_between].0
        {
            index_between += 1;
        }

        // Get pchip estimate for the value...
        let indexes = [
            index_between.saturating_sub(2),
            index_between.saturating_sub(1),
            index_between,
            index_between
                .saturating_add(1)
                .min(dual_est_quants.len() - 1),
        ];

        new_estimator.ranks[ti] = approx_obs_rank;
        new_estimator.values[ti] = pchip_prediction(
            &indexes.map(|i| dual_est_quants[i].0),
            &indexes.map(|i| dual_est_quants[i].1),
            approx_obs_rank as f64,
        )
    }
}

pub trait QuantileEstimator {
    fn update(&mut self, sample: f64);
    fn update_all(&mut self, samples: &[f64]) {
        for &s in samples.iter() {
            self.update(s);
        }
    }
    fn combine(&self, other: &Self) -> Self;
}

#[derive(Clone)]
struct FixedSizeQuantileEstimator<const N: usize> {
    values: [f64; N],
    ranks: [usize; N],
    targets: [f64; N],
    observations: usize,
}

impl<const N: usize> FixedSizeQuantileEstimator<N> {
    pub fn new(targets: &[f64; N]) -> Self {
        Self {
            values: [0.0; N],
            ranks: [0; N],
            targets: targets.clone(),
            observations: 0,
        }
    }

    fn _as_data(&self) -> QuantileEstimatorData<'_> {
        QuantileEstimatorData {
            ranks: &self.ranks,
            values: &self.values,
            targets: &self.targets,
            observations: &self.observations,
        }
    }

    fn _as_mut_data(&mut self) -> MutableQuantileEstimatorData<'_> {
        MutableQuantileEstimatorData {
            ranks: &mut self.ranks,
            values: &mut self.values,
            targets: &self.targets,
            observations: &mut self.observations,
        }
    }

    fn _is_initialized(&self) -> bool {
        self.observations >= N
    }
}

impl<const N: usize> QuantileEstimator for FixedSizeQuantileEstimator<N> {
    fn update(&mut self, sample: f64) {
        match (self.observations + 1).cmp(&self.values.len()) {
            Ordering::Less => {
                self.values[self.observations] = sample;
                self.observations += 1;
            }
            Ordering::Equal => {
                self.values[self.observations] = sample;
                self.values.sort_by(|a, b| a.total_cmp(b));
                for i in 0..self.ranks.len() {
                    self.ranks[i] = i / (self.ranks.len() - 1);
                }
                self.observations += 1;
            }
            Ordering::Greater => {
                _add_sample_to_estimator(self._as_mut_data(), sample);
            }
        }
    }

    fn combine(&self, other: &Self) -> Self {
        match (self._is_initialized(), other._is_initialized()) {
            (true, true) => {
                let mut new_quant_est = Self::new(&self.targets);

                _merge_estimators(
                    self._as_data(),
                    other._as_data(),
                    new_quant_est._as_mut_data(),
                );

                new_quant_est
            }
            (true, false) | (false, false) => {
                let mut new_quant_est = self.clone();
                new_quant_est.update_all(&other.values[..other.observations]);
                new_quant_est
            }
            (false, true) => {
                let mut new_quant_est = other.clone();
                new_quant_est.update_all(&self.values[..self.observations]);
                new_quant_est
            }
        }
    }
}
