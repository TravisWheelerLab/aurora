/// Implementation of P2 estimator.
/// See "The P2 Algorithm for Dynamic Statistical Computing Calculation of Quantiles and Histograms Without Storing Observations"
/// at https://www.cse.wustl.edu/~jain/papers/ftp/psqr.pdf
use num_traits::{AsPrimitive, Float, Unsigned};

struct P2HistogramPoint<F: Float, I: Unsigned> {
    value: F,
    rank: I,
}

fn linear_prediction<F: Float + From<isize>, I: Unsigned + Copy + Ord + Into<F>>(
    points: &[P2HistogramPoint<F, I>; 3],
    d: isize,
) -> F {
    let n: [F; 3] = points.each_ref().map(|v| v.rank.into());
    let q: [F; 3] = points.each_ref().map(|v| v.value);
    let d_f: F = d.into();
    let d_off = 1 + d as usize;

    q[1] + d_f * ((q[d_off] - q[1]) / (n[d_off] - n[1]))
}

fn parabolic_prediction<F: Float + From<isize>, I: Unsigned + Copy + Ord + Into<F>>(
    points: &[P2HistogramPoint<F, I>; 3],
    d: isize,
) -> F {
    let n: [F; 3] = points.each_ref().map(|v| v.rank.into());
    let q: [F; 3] = points.each_ref().map(|v| v.value);
    let d = d.into();

    let left = (n[1] - n[0] + d) * ((q[2] - q[1]) / (n[2] - n[1]));
    let right = (n[2] - n[1] - d) * ((q[1] - q[0]) / (n[1] - n[0]));
    q[1] + (d / (n[2] - n[0])) * (left + right)
}

fn p2update<F: Float + AsPrimitive<isize>, I: Unsigned + Ord + Into<F>>(
    points: &[P2HistogramPoint<F, I>],
    center_index: I,
    observations: I,
    total_points: I,
    proposal: F,
) -> (F, I) {
    // Actual rank desired for the given quantile...
    let rank_proposal: F =
        (center_index * (observations - I::one())).into() / (total_points - I::one()).into();
    let d: F = rank_proposal - points[1].rank.into();

    if d >= F::one() && (points[2].rank - points[1].rank) > I::one()
        || (d <= -F::one()) && points[1].rank - points[0].rank > I::one()
    {
        let d: isize = if d >= F::zero() { 1 } else { -1 };
        let n: [F; 3];
        for i in 0..3 {
            n[i] = points[i].rank.into();
        }

        let est = points[1].value + (d / (points[2].rank - points[0].rank).into());
        est
    } else {
        (points[1].value, I::zero())
    }
}

struct P2HistogramData<'a, F: Float, I: Unsigned + Into<F>> {
    observations: I,
    points: &'a mut [P2HistogramPoint<F, I>],
}
