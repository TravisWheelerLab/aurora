use std::cmp::Ordering;

use crate::statistics::Distribution;
use itertools::Itertools;
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
    let h1 = ((t - 2.0) * t + 1.0) * t;
    let h2 = (-2.0 * t + 3.0) * t * t;
    let h3 = (t - 1.0) * t * t;

    h0 * y0 + h1 * ms0 + h2 * y1 + h3 * ms1
}

fn pchip_point_derivative(dx0: f64, dy0: f64, dx1: f64, dy1: f64) -> f64 {
    let s0 = if dx0 != 0.0 { dy0 / dx0 } else { 0.0 };
    let s1 = if dx1 != 0.0 { dy1 / dx1 } else { 0.0 };
    if s0 * s1 > 0.0 {
        let alpha = (1.0 / 3.0) * (1.0 + dx1 / (dx0 + dx1));
        s0 * s1 / (alpha * s1 + (1.0 - alpha) * s0)
    } else {
        0.0
    }
}

fn pchip_prediction(x_points: &[f64; 4], y_points: &[f64; 4], x: f64) -> f64 {
    debug_assert!(x_points.is_sorted() && x <= x_points[2] && x >= x_points[1]);
    let m0 = pchip_point_derivative(
        x_points[1] - x_points[0],
        y_points[1] - y_points[0],
        x_points[2] - x_points[1],
        y_points[2] - y_points[1],
    );
    let m1 = pchip_point_derivative(
        x_points[2] - x_points[1],
        y_points[2] - y_points[1],
        x_points[3] - x_points[2],
        y_points[3] - y_points[2],
    );

    cubic_hermite_spline(
        x_points[1],
        y_points[1],
        x_points[2],
        y_points[2],
        m0,
        m1,
        x,
    )
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
    debug_assert!(ranks.first() == Some(&0) && ranks.last() == Some(&(observations - 1)));
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
    let p = values.partition_point(|&v| v < sample);
    let bound_p = p.min(values.len() - 1);

    // Update extremes...
    if bound_p == 0 {
        values[bound_p] = values[bound_p].min(sample);
    } else if bound_p == (values.len() - 1) {
        values[bound_p] = values[bound_p].max(sample);
    }

    // Increment ranks of markers above newly inserted sample...
    for i in bound_p.max(1)..ranks.len() {
        ranks[i] = ranks[i] + 1;
    }

    // Adjust inner markers to within 1 of their target quantile using p2 formula...
    for i in 1..(values.len() - 1) {
        // Observations hasn't been incremented yet, don't need to subtract 1...
        let target_rank = (targets[i] * (*observations) as f64) as usize;
        let current_rank = ranks[i];
        if current_rank.abs_diff(target_rank) > 1 {
            //println!("{:?}, {}, {}", ranks, i, ranks[i]);
            let new_rank: usize = target_rank.clamp(
                ranks[i - 1].saturating_add(1),
                ranks[i + 1].saturating_sub(1),
            );
            if new_rank == current_rank {
                continue;
            }

            let idx_shift = if new_rank > current_rank { 1 } else { 0 };

            let indexes = [
                (i + idx_shift).saturating_sub(2),
                (i + idx_shift).saturating_sub(1),
                (i + idx_shift),
                (i + idx_shift).saturating_add(1).min(ranks.len() - 1),
            ];

            values[i] = pchip_prediction(
                &indexes.map(|i| ranks[i] as f64),
                &indexes.map(|i| values[i]),
                new_rank as f64,
            );
            ranks[i] = new_rank;
        }
    }

    *observations += 1;
}

fn _merge_estimators(
    q1: QuantileEstimatorData,
    q2: QuantileEstimatorData,
    new_estimator: MutableQuantileEstimatorData,
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
            dual_est_quants.push((
                (next.0 + q1_prior.map(|v| v.0 + 1).unwrap_or(0)) as f64,
                next.1,
            ));
            q2_prior = Some(next);
            q2_idx += 1;
        } else if q2_past_end {
            let next = get_at(&q1, q1_idx);
            dual_est_quants.push((
                (next.0 + q2_prior.map(|v| v.0 + 1).unwrap_or(0)) as f64,
                next.1,
            ));
            q1_prior = Some(next);
            q1_idx += 1;
        } else if q1.values[q1_idx] <= q2.values[q2_idx] {
            let other_next = get_at(&q2, q2_idx);
            let next = get_at(&q1, q1_idx);
            let w = q2_prior
                .map(|other_prior| (next.1 - other_prior.1) / (other_next.1 - other_prior.1))
                .unwrap_or(0.0);
            let other_rank_est = q2_prior
                .map(|other_prior| {
                    (other_prior.0 + 1) as f64 * (1.0 - w) + (other_next.0 + 1) as f64 * w
                })
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
                .map(|other_prior| {
                    (other_prior.0 + 1) as f64 * (1.0 - w) + (other_next.0 + 1) as f64 * w
                })
                .unwrap_or(0.0);
            dual_est_quants.push((next.0 as f64 + other_rank_est, next.1));
            q2_prior = Some(next);
            q2_idx += 1;
        }
    }

    // New number of observations is the sum of both...
    *(new_estimator.observations) = *(q1.observations) + *(q2.observations);
    let rank_range = *(new_estimator.observations) - 1;

    // Solve all inner quantiles using traditional interpolation...
    let mut index_between = 0;

    new_estimator.ranks.first_mut().map(|r| *r = 0);
    new_estimator.ranks.last_mut().map(|r| *r = rank_range);
    new_estimator
        .values
        .first_mut()
        .map(|v| *v = dual_est_quants[0].1);
    new_estimator
        .values
        .last_mut()
        .map(|v| *v = dual_est_quants[dual_est_quants.len() - 1].1);

    for ti in 1..new_estimator.targets.len() - 1 {
        // Calculate new rank...
        let target = new_estimator.targets[ti];
        let approx_obs_rank = ((target * rank_range as f64) as usize)
            .clamp(ti, rank_range - (new_estimator.targets.len() - (ti + 1)));

        // Find where it lands in cdf...
        while index_between < dual_est_quants.len()
            && (approx_obs_rank as f64) > dual_est_quants[index_between].0
        {
            index_between += 1;
        }

        // Get pchip estimate for the value...
        let indexes = [
            index_between.saturating_sub(2),
            index_between.saturating_sub(1),
            index_between.min(dual_est_quants.len() - 1),
            index_between
                .saturating_add(1)
                .min(dual_est_quants.len() - 1),
        ];

        new_estimator.values[ti] = pchip_prediction(
            &indexes.map(|i| dual_est_quants[i].0),
            &indexes.map(|i| dual_est_quants[i].1),
            approx_obs_rank as f64,
        );
        new_estimator.ranks[ti] = approx_obs_rank;
    }
}

trait PrimativeCast<T> {
    fn as_(&self) -> T;
}

impl PrimativeCast<f64> for f64 {
    #[inline]
    fn as_(&self) -> f64 {
        *self
    }
}

impl PrimativeCast<f64> for usize {
    #[inline]
    fn as_(&self) -> f64 {
        *self as f64
    }
}

fn _interpolated_value_prediction<
    I: PartialOrd + PrimativeCast<f64> + Copy,
    O: PartialOrd + PrimativeCast<f64> + Copy,
>(
    xs: &[I],
    ys: &[O],
    x: f64,
    lower_val: f64,
    upper_val: f64,
    not_enough_data_value: f64,
) -> f64 {
    debug_assert!(xs.is_sorted());
    debug_assert!(xs.len() == ys.len());

    if xs.len() < 1 {
        return not_enough_data_value;
    }

    let idx = xs.partition_point(|&v| v.as_() < x);
    if idx > xs.len() {
        upper_val
    } else if idx == 0 {
        lower_val
    } else {
        let indexes = [
            idx.saturating_sub(2),
            idx.saturating_sub(1),
            idx,
            idx.saturating_add(1).min(xs.len() - 1),
        ];

        pchip_prediction(
            &indexes.map(|i| xs[i].as_()),
            &indexes.map(|i| ys[i].as_()),
            x,
        )
    }
}

pub trait QuantileEstimator: Distribution {
    fn from_prior(prior: &Self, count: usize) -> Self;
    fn update(&mut self, sample: f64);
    fn update_all(&mut self, samples: &[f64]) {
        for &s in samples.iter() {
            self.update(s);
        }
    }
    fn combine(&self, other: &Self) -> Self;
    #[allow(dead_code)]
    fn samples(&self) -> usize;
}

trait SimpleQuantileEstimatorRepresentation: Clone {
    fn new_like(other: &Self) -> Self;
    fn _data(&self) -> QuantileEstimatorData<'_>;
    fn _mut_data(&mut self) -> MutableQuantileEstimatorData<'_>;
    fn _is_initialized(&self) -> bool {
        let data = self._data();
        *data.observations >= data.ranks.len()
    }
}

impl<Q: SimpleQuantileEstimatorRepresentation> QuantileEstimator for Q {
    fn samples(&self) -> usize {
        *self._data().observations
    }

    fn from_prior(prior: &Self, count_per_entry: usize) -> Self {
        let prior_data = prior._data();
        let mut new_self = Self::new_like(prior);
        let new_data = new_self._mut_data();

        let new_observations = count_per_entry.max(1) * prior_data.ranks.len();

        for i in 0..new_data.targets.len() {
            let closest_rank = ((new_data.targets[i] * (new_observations - 1) as f64) as usize)
                .clamp(
                    i,
                    (new_observations - 1) - (new_data.targets.len() - (i + 1)),
                );

            new_data.ranks[i] = closest_rank;
            new_data.values[i] = prior.ppf(closest_rank as f64 / (new_observations - 1) as f64)
        }
        *new_data.observations = new_observations;

        new_self
    }

    fn update(&mut self, sample: f64) {
        let data = self._mut_data();

        match (*data.observations + 1).cmp(&data.values.len()) {
            Ordering::Less => {
                data.values[*data.observations] = sample;
                *data.observations += 1;
            }
            Ordering::Equal => {
                data.values[*data.observations] = sample;
                data.values.sort_by(|a, b| a.total_cmp(b));
                for i in 0..data.ranks.len() {
                    data.ranks[i] = i;
                }
                *data.observations += 1;
            }
            Ordering::Greater => {
                _add_sample_to_estimator(data, sample);
            }
        }
    }

    fn combine(&self, other: &Self) -> Self {
        match (self._is_initialized(), other._is_initialized()) {
            (true, true) => {
                let mut new_quant_est = Self::new_like(&self);

                _merge_estimators(self._data(), other._data(), new_quant_est._mut_data());

                new_quant_est
            }
            (true, false) | (false, false) => {
                let other_data = other._data();
                let mut new_quant_est = self.clone();
                new_quant_est.update_all(&other_data.values[..*other_data.observations]);
                new_quant_est
            }
            (false, true) => {
                let self_data = self._data();
                let mut new_quant_est = other.clone();
                new_quant_est.update_all(&self_data.values[..*self_data.observations]);
                new_quant_est
            }
        }
    }
}

impl<Q: SimpleQuantileEstimatorRepresentation> Distribution for Q {
    fn cdf(&self, x: f64) -> f64 {
        let data = self._data();
        if self._is_initialized() {
            _interpolated_value_prediction(
                data.values,
                data.ranks,
                x,
                0.0,
                (*data.observations - 1) as f64,
                0.0,
            ) / (*data.observations - 1).max(1) as f64
        } else {
            let xs_sorted = data.values[..*data.observations]
                .iter()
                .copied()
                .sorted_by(|a, b| a.total_cmp(b))
                .collect_vec();
            let ys = (0..*data.observations).collect_vec();
            _interpolated_value_prediction(
                &xs_sorted,
                &ys,
                x,
                0.0,
                (*data.observations - 1) as f64,
                0.0,
            )
        }
    }

    fn logcdf(&self, x: f64) -> f64 {
        self.cdf(x).ln()
    }

    fn ccdf(&self, x: f64) -> f64 {
        1.0 - self.cdf(x)
    }

    fn logccdf(&self, x: f64) -> f64 {
        (-self.cdf(x)).ln_1p()
    }

    fn ppf(&self, p: f64) -> f64 {
        let data = self._data();
        let est_rank = p.clamp(0.0, 1.0) * (*data.observations - 1) as f64;

        let data = self._data();
        let (min_val, max_val) = self.support();

        if self._is_initialized() {
            _interpolated_value_prediction(
                data.ranks,
                data.values,
                est_rank,
                min_val,
                max_val,
                0.0_f64.clamp(min_val, max_val),
            )
        } else {
            let ys_sorted = data.values[..*data.observations]
                .iter()
                .copied()
                .sorted_by(|a, b| a.total_cmp(b))
                .collect_vec();
            let xs = (0..*data.observations).collect_vec();
            _interpolated_value_prediction(
                &xs,
                &ys_sorted,
                est_rank,
                min_val,
                max_val,
                0.0_f64.clamp(min_val, max_val),
            )
        }
    }

    fn pdf(&self, _x: f64) -> f64 {
        // Will have to calculate derivatives, cache normalization factor (such that area under curve is 1)...
        // May be worth splitting out into different class to allow pre-processing this stuff...
        // TODO: Would prefer quintic splines for this... Allows us to avoid normalization...
        panic!("Currently not supported!");
    }

    fn logpdf(&self, x: f64) -> f64 {
        self.pdf(x).ln()
    }

    fn support(&self) -> (f64, f64) {
        let data = self._data();

        if self._is_initialized() {
            (
                *data.values.first().unwrap_or(&f64::NEG_INFINITY),
                *data.values.last().unwrap_or(&f64::INFINITY),
            )
        } else {
            data.values[..*data.observations]
                .iter()
                .copied()
                .minmax()
                .into_option()
                .unwrap_or((f64::NEG_INFINITY, f64::INFINITY))
        }
    }
}

#[derive(Clone, Debug)]
pub struct FixedSizeQuantileEstimator<const N: usize> {
    values: [f64; N],
    ranks: [usize; N],
    targets: [f64; N],
    observations: usize,
}

impl<const N: usize> FixedSizeQuantileEstimator<N> {
    pub fn new(targets: &[f64; N]) -> Self {
        assert!(
            targets.is_sorted() && targets.first() == Some(&0.0) && targets.last() == Some(&1.0)
        );
        Self {
            values: [0.0; N],
            ranks: [0; N],
            targets: targets.clone(),
            observations: 0,
        }
    }
}

impl<const N: usize> SimpleQuantileEstimatorRepresentation for FixedSizeQuantileEstimator<N> {
    fn new_like(other: &Self) -> Self {
        Self::new(&other.targets)
    }

    fn _data(&self) -> QuantileEstimatorData<'_> {
        QuantileEstimatorData {
            ranks: &self.ranks,
            values: &self.values,
            targets: &self.targets,
            observations: &self.observations,
        }
    }

    fn _mut_data(&mut self) -> MutableQuantileEstimatorData<'_> {
        MutableQuantileEstimatorData {
            ranks: &mut self.ranks,
            values: &mut self.values,
            targets: &self.targets,
            observations: &mut self.observations,
        }
    }
}

#[derive(Clone, Debug)]
pub struct VectorQuantileEstimator {
    values: Vec<f64>,
    ranks: Vec<usize>,
    targets: Vec<f64>,
    observations: usize,
}

impl VectorQuantileEstimator {
    pub fn new(targets: &[f64]) -> Self {
        assert!(
            targets.is_sorted() && targets.first() == Some(&0.0) && targets.last() == Some(&1.0)
        );
        Self {
            values: vec![0.0; targets.len()],
            ranks: (0..targets.len()).collect_vec(),
            targets: Vec::from(targets),
            observations: 0,
        }
    }
}

impl SimpleQuantileEstimatorRepresentation for VectorQuantileEstimator {
    fn new_like(other: &Self) -> Self {
        Self::new(&other.targets)
    }

    fn _data(&self) -> QuantileEstimatorData<'_> {
        QuantileEstimatorData {
            ranks: &self.ranks,
            values: &self.values,
            targets: &self.targets,
            observations: &self.observations,
        }
    }

    fn _mut_data(&mut self) -> MutableQuantileEstimatorData<'_> {
        MutableQuantileEstimatorData {
            ranks: &mut self.ranks,
            values: &mut self.values,
            targets: &self.targets,
            observations: &mut self.observations,
        }
    }
}

pub mod custom_quantile_estimator {
    use super::*;
    use std::f64::consts::E;

    macro_rules! replace_expr {
        ($_t:tt,$sub:expr) => {
            $sub
        };
    }

    macro_rules! count_exprs {
        ($($val:expr),+) => {<[()]>::len(&[$(replace_expr!($val,())),+])};
    }

    macro_rules! implement_fixed_quantile_estimator {
        ($name:ident[$($val:expr),+]) => {
            #[derive(Clone, Debug)]
            pub struct $name {
                values: [f64; Self::COUNT],
                ranks: [usize; Self::COUNT],
                observations: usize,
            }

            impl $name {
                const TARGETS: [f64; count_exprs!($($val),+) + 2] = [0.0, $($val),+, 1.0];
                const COUNT: usize = Self::TARGETS.len();

                pub fn new() -> Self {
                    Self {
                        values: [0.0; _],
                        ranks: [0; _],
                        observations: 0
                    }
                }
            }

            impl Default for $name {
                fn default() -> Self {
                    Self::new()
                }
            }

            impl SimpleQuantileEstimatorRepresentation for $name {
                fn new_like(_other: &Self) -> Self {
                    Self::default()
                }
                fn _data(&self) -> QuantileEstimatorData<'_> {
                    QuantileEstimatorData {
                        ranks: &self.ranks,
                        values: &self.values,
                        targets: &Self::TARGETS,
                        observations: &self.observations,
                    }
                }
                fn _mut_data(&mut self) -> MutableQuantileEstimatorData<'_> {
                    MutableQuantileEstimatorData {
                        ranks: &mut self.ranks,
                        values: &mut self.values,
                        targets: &Self::TARGETS,
                        observations: &mut self.observations,
                    }
                }
            }
        };
    }

    implement_fixed_quantile_estimator!(FrechetQuant[0.5 / E, 0.25, 1.0 / E, 0.5, 0.5 + 1.0 / 2.0 * E, 0.75]);
}

#[cfg(test)]
mod test {
    use crate::{
        p2estimator::{FixedSizeQuantileEstimator, QuantileEstimator, VectorQuantileEstimator},
        statistics::{linspace, Distribution, Exponential},
    };
    use itertools::Itertools;
    use rand::{rngs::Xoshiro256PlusPlus, RngExt, SeedableRng};

    fn is_close(a: f64, b: f64) -> bool {
        let rel_tol = 1e-9;
        let abs_tol = 0.0;
        if a == b {
            true
        } else {
            (a - b).abs() <= (rel_tol * (a.abs()).max(b.abs())).max(abs_tol)
        }
    }

    #[test]
    fn quantiles_on_exponential_dist() {
        let expon = Exponential::new(1.0);
        let mut estimator =
            FixedSizeQuantileEstimator::new(&[0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]);

        let mut rng = Xoshiro256PlusPlus::seed_from_u64(12345654321);

        for _ in 0..10_000 {
            let sample = expon.ppf(rng.random());

            estimator.update(sample);
        }

        assert!(estimator.samples() == 10_000);

        for val in linspace(0.0, 0.90, 90) {
            assert!((estimator.ppf(val) - expon.ppf(val)).abs() <= 0.04);
            let dist_val = expon.ppf(val);
            assert!((estimator.cdf(dist_val) - expon.cdf(dist_val)).abs() <= 0.04);

            // Basic probability distribution checks...
            assert!(is_close(
                estimator.ccdf(dist_val),
                1.0 - estimator.cdf(dist_val)
            ));

            assert!(is_close(
                estimator.logcdf(dist_val),
                estimator.cdf(dist_val).ln()
            ));
            assert!(is_close(
                estimator.logccdf(dist_val),
                estimator.ccdf(dist_val).ln()
            ));
        }
    }

    #[test]
    fn test_quantile_merging() {
        let expon = Exponential::new(1.0);
        let mut merged_estimator =
            VectorQuantileEstimator::new(&linspace(0.0, 1.0, 10).collect_vec());

        let mut rng = Xoshiro256PlusPlus::seed_from_u64(12345654321);

        for _ in 0..100 {
            let targets: Vec<f64> = linspace(0.0, 1.0, rng.random_range(5..15)).collect();
            let mut estimator = VectorQuantileEstimator::new(&targets);

            for _ in 0..100 {
                estimator.update(expon.ppf(rng.random()));
            }

            merged_estimator = merged_estimator.combine(&estimator);
        }

        assert!(merged_estimator.samples() == 10_000);

        for val in linspace(0.0, 0.75, 75) {
            assert!((merged_estimator.ppf(val) - expon.ppf(val)).abs() <= 0.1);
            let dist_val = expon.ppf(val);
            assert!((merged_estimator.cdf(dist_val) - expon.cdf(dist_val)).abs() <= 0.04)
        }
    }
}
