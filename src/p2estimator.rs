use std::{cmp::Ordering, ops::Neg};

/// Implementation of P2 estimator.
/// See "The P2 Algorithm for Dynamic Statistical Computing Calculation of Quantiles and Histograms Without Storing Observations"
/// at https://www.cse.wustl.edu/~jain/papers/ftp/psqr.pdf
use num_traits::{float::TotalOrder, Float, Num, Unsigned};

struct P2HistogramPoint<F: Float, I: Unsigned> {
    value: F,
    rank: I,
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

fn linear_prediction<F: Float, I: Unsigned + Copy + Ord + Into<F>>(
    points: &[P2HistogramPoint<F, I>; 3],
    d: isize,
) -> F {
    let n: [F; 3] = points.each_ref().map(|v| v.rank.into());
    let q: [F; 3] = points.each_ref().map(|v| v.value);
    let d_f: F = get_sign(d);
    let d_off = (1 + d) as usize;

    q[1] + d_f * ((q[d_off] - q[1]) / (n[d_off] - n[1]))
}

fn parabolic_prediction<F: Float, I: Unsigned + Copy + Ord + Into<F>>(
    points: &[P2HistogramPoint<F, I>; 3],
    d: isize,
) -> F {
    let n: [F; 3] = points.each_ref().map(|v| v.rank.into());
    let q: [F; 3] = points.each_ref().map(|v| v.value);
    let d: F = get_sign(d);

    let left = (n[1] - n[0] + d) * ((q[2] - q[1]) / (n[2] - n[1]));
    let right = (n[2] - n[1] - d) * ((q[1] - q[0]) / (n[1] - n[0]));
    q[1] + (d / (n[2] - n[0])) * (left + right)
}

fn _p2update<F: Float, I: Unsigned + Copy + Ord + Into<F> + From<usize>>(
    points: &mut [P2HistogramPoint<F, I>],
    center_index: usize,
    observations: I,
    total_points: I,
) {
    // Actual rank desired for the given quantile...
    let ci: I = center_index.into();
    let rank_proposal: F =
        (ci * (observations - I::one())).into() / (total_points - I::one()).into();
    let d: F = rank_proposal - points[1].rank.into();

    if d >= F::one() && (points[2].rank - points[1].rank) > I::one()
        || (d <= -F::one()) && points[1].rank - points[0].rank > I::one()
    {
        let d: isize = get_sign(d);
        let mut p_est = parabolic_prediction(
            (&points[center_index - 1..center_index + 1])
                .as_array()
                .unwrap(),
            d,
        );
        if p_est <= points[center_index - 1].value || p_est >= points[center_index + 1].value {
            p_est = linear_prediction(
                (&points[center_index - 1..center_index + 1])
                    .as_array()
                    .unwrap(),
                d,
            );
        }

        points[center_index].value = p_est;
        points[center_index].rank = inc_or_dec(points[center_index].rank, d);
    }
}

struct P2HistogramData<'a, F: Float, I: Unsigned + Copy + Ord + Into<F>> {
    observations: I,
    points: &'a mut [P2HistogramPoint<F, I>],
}

impl<'a, F: Float + TotalOrder, I: Unsigned + Copy + Ord + Into<F> + From<usize> + Into<usize>>
    P2HistogramData<'a, F, I>
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
            _p2update(self.points, i, self.observations, self.points.len().into());
        }

        self.observations = self.observations + I::one();
    }

    fn _pre_init_update(&mut self, sample: F) {
        let nxt_idx: usize = self.observations.into();
        self.points[nxt_idx].value = sample;
        self.observations = self.observations + I::one();
    }

    fn _initialize(&mut self) {
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
