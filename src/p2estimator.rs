use std::{cmp::Ordering, ops::Neg};

/// Implementation of P2 estimator.
/// See "The P2 Algorithm for Dynamic Statistical Computing Calculation of Quantiles and Histograms Without Storing Observations"
/// at https://www.cse.wustl.edu/~jain/papers/ftp/psqr.pdf
///
/// We replace the P2 interpolation with PCHIP instead (See paper A Method for Constructing Local Monotone Piecewise Cubic Interpolants by F. N. Fritsch and J. Butland, or https://doi.org/10.1137/0905021)
use num_traits::{float::TotalOrder, Float, FromPrimitive, Num, Unsigned};

struct P2HistogramPoint<F: Float, I: Unsigned> {
    value: F,
    rank: I,
    target: F,
}

fn get_sign<A: Num + PartialOrd + Neg<Output = A>, B: Num + PartialOrd + Neg<Output = B>>(
    val: A,
) -> B {
    if val >= A::zero() {
        B::one()
    } else {
        -B::one()
    }
}

fn inc_or_dec<A: Num + PartialOrd, B: Num + PartialOrd + Neg<Output = B>>(val: A, delta: B) -> A {
    match delta.partial_cmp(&B::zero()) {
        Some(Ordering::Less) => val - A::one(),
        Some(Ordering::Greater) => val + A::one(),
        _ => val,
    }
}

fn cubic_hermite_spline<F: Float + FromPrimitive>(
    x0: F,
    y0: F,
    x1: F,
    y1: F,
    m0: F,
    m1: F,
    x: F,
) -> F {
    let t = (x - x0) / (x1 - x0);
    let ms0 = (x1 - x0) * m0;
    let ms1 = (x1 - x0) * m1;

    let _1 = F::one();
    let _2 = F::from_i32(2).unwrap();
    let _3 = F::from_i32(3).unwrap();

    let h0: F = (_2 * t - _3) * t * t + _1;
    let h1 = ((t - _1) * t + _1) * t;
    let h2 = (_2 * t + _3) * t * t;
    let h3 = (t - _1) * t * t;

    h0 * y0 + h1 * ms0 + h2 * y1 + h3 * ms1
}

fn pchip_point_derivative<F: Float + FromPrimitive>(dx0: F, dy0: F, dx1: F, dy1: F) -> F {
    if dy0 * dy1 > F::zero() {
        let _1 = F::one();
        let one_third = _1 / F::from_i32(3).unwrap();
        let alpha = one_third * (_1 + dx1 / (dx0 + dx1));
        dy0 * dy1 / (alpha * dy1 + (_1 - alpha) * dy0)
    } else {
        F::zero()
    }
}

fn secant_diff<F: Float, I: Unsigned + Copy + Ord + Into<F>>(
    point0: Option<&P2HistogramPoint<F, I>>,
    point1: Option<&P2HistogramPoint<F, I>>,
) -> (F, F) {
    if let (Some(p0), Some(p1)) = (point0, point1) {
        ((p1.rank - p0.rank).into(), p1.value - p0.value)
    } else {
        // Assume slope at endpoints of CDF is 0...
        (F::zero(), F::zero())
    }
}

fn pchip_prediction<F: Float + FromPrimitive, I: Unsigned + Copy + Ord + Into<F>>(
    point0: Option<&P2HistogramPoint<F, I>>,
    point1: &P2HistogramPoint<F, I>,
    point2: &P2HistogramPoint<F, I>,
    point3: Option<&P2HistogramPoint<F, I>>,
    x: F,
) -> F {
    let s0 = secant_diff(point0, Some(point1));
    let s1 = secant_diff(Some(point1), Some(point2));
    let s2 = secant_diff(Some(point2), point3);
    let m0 = pchip_point_derivative(s0.0, s0.1, s1.0, s1.1);
    let m1 = pchip_point_derivative(s1.0, s1.1, s2.0, s2.1);

    cubic_hermite_spline(
        point1.rank.into(),
        point1.value,
        point2.rank.into(),
        point2.value,
        m0,
        m1,
        x,
    )
}

struct QuantileEstimator<'a, F: Float, I: Unsigned + Copy + Ord + Into<F>> {
    observations: I,
    points: &'a mut [P2HistogramPoint<F, I>],
}

impl<
        'a,
        F: Float + TotalOrder + Into<usize> + FromPrimitive,
        I: Unsigned + Copy + Ord + Into<F> + From<usize> + Into<usize>,
    > QuantileEstimator<'a, F, I>
{
    fn _standard_update(&mut self, sample: F) {
        // Find where sample falls within distribution...
        let p = self.points.partition_point(|v| v.value <= sample);
        let bound_p = p.min(self.points.len() - 1);

        // Update extremes...
        if bound_p == 0 {
            self.points[bound_p].value = self.points[bound_p].value.min(sample);
        } else if bound_p == (self.points.len() - 1) {
            self.points[bound_p].value = self.points[bound_p].value.max(sample);
        }

        // Increment ranks of markers above newly inserted sample...
        for i in (bound_p + 1)..self.points.len() {
            self.points[i].rank = self.points[i].rank + I::one();
        }

        // Adjust inner markers to within 1 of their target quantile using p2 formula...
        for i in 1..(self.points.len() - 1) {
            let target_rank: F = (self.points[i].target * self.observations.into()).floor();
            let true_rank: F = self.points[i].rank.into();
            let true_rank_int: usize = self.points[i].rank.into();
            if (true_rank - target_rank).abs() > F::one() {
                let target_rank_int: usize = target_rank.into();
                let lower_rank: usize = self.points[i - 1].rank.into();
                let upper_rank: usize = self.points[i + 1].rank.into();
                let new_rank: usize = target_rank_int
                    .clamp(lower_rank.saturating_add(1), upper_rank.saturating_sub(1));
                if new_rank == true_rank_int {
                    continue;
                }

                let shift = if new_rank > target_rank_int { 1 } else { 0 };

                self.points[i].rank = new_rank.into();
                self.points[i].value = pchip_prediction(
                    if i + shift > 2 {
                        Some(&self.points[i + shift - 2])
                    } else {
                        None
                    },
                    &self.points[i - shift - 1],
                    &self.points[i + shift],
                    if i + shift + 1 < self.points.len() {
                        Some(&self.points[i + shift + 1])
                    } else {
                        None
                    },
                    F::from_usize(new_rank).unwrap(),
                );
            }
        }

        self.observations = self.observations + I::one();
    }

    fn _pre_init_update(&mut self, sample: F) {
        let nxt_idx: usize = self.observations.into();
        self.points[nxt_idx].value = sample;
        self.observations = self.observations + I::one();
    }

    fn _initialize(&mut self) {
        panic!("Fix!");
        // TODO: Fix this...
        self.points.sort_by(|a, b| a.value.total_cmp(&b.value));
        self.points.iter_mut().enumerate().for_each(|(i, p)| {
            p.rank = i.into();
        });
    }

    pub fn update(&mut self, sample: F) {
        let obs: usize = self.observations.into();
        match (obs + 1).cmp(&self.points.len()) {
            Ordering::Less => self._pre_init_update(sample),
            Ordering::Equal => {
                self._pre_init_update(sample);
                self._initialize();
            }
            Ordering::Greater => {
                self._standard_update(sample);
            }
        }
    }

    pub fn is_initialized(&self) -> bool {
        let obs: usize = self.observations.into();
        obs >= self.points.len()
    }

    fn combine(&mut self, other: &P2HistogramData<F, I>) {
        // TODO: Need to think about how to do this efficiently while maintaining accuracy...
        panic!("Not implemented!")
    }
}
