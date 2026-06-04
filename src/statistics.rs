use core::f64;
use puruspe::{beta, betai, invbetai};
use std::fmt::Debug;

pub fn ln_add_exp(a: f64, b: f64) -> f64 {
    let max = a.max(b);
    let min = a.min(b);
    // TODO: Possibly use more stable ln_1p_exp at https://github.com/JuliaStats/LogExpFunctions.jl/files/8218470/log1pexp.pdf (Implemented at https://github.com/JuliaStats/LogExpFunctions.jl/blob/master/src/basicfuns.jl#L263)
    max + (min - max).exp().ln_1p()
}

// TODO: Support for generic floating types...
#[allow(dead_code)]
pub trait Distribution {
    fn pdf(&self, x: f64) -> f64;
    fn cdf(&self, x: f64) -> f64;
    fn ppf(&self, p: f64) -> f64;
    fn support(&self) -> (f64, f64);
    fn ccdf(&self, x: f64) -> f64;
    fn logpdf(&self, x: f64) -> f64;
    fn logcdf(&self, x: f64) -> f64;
    fn logccdf(&self, x: f64) -> f64;
}

#[allow(dead_code)]
pub trait ParameterizedDistribution: Distribution + Debug + Default + Clone {
    fn unit() -> Self {
        Self::default()
    }
}

#[derive(Clone, Debug)]
pub struct Exponential {
    lambda: f64,
}

impl ParameterizedDistribution for Exponential {}

impl Exponential {
    pub fn new(lambda: f64) -> Self {
        Self { lambda }
    }

    pub fn from_scale(beta: f64) -> Self {
        Self::new(1.0 / beta)
    }
}

impl Default for Exponential {
    fn default() -> Self {
        Self::new(1.0)
    }
}

impl Distribution for Exponential {
    fn pdf(&self, x: f64) -> f64 {
        self.lambda * (-self.lambda * x).exp()
    }

    fn cdf(&self, x: f64) -> f64 {
        1.0 - (-self.lambda * x).exp()
    }

    fn ppf(&self, p: f64) -> f64 {
        -(1.0 - p).ln() / self.lambda
    }

    fn ccdf(&self, x: f64) -> f64 {
        (-self.lambda * x).exp()
    }

    fn logpdf(&self, x: f64) -> f64 {
        self.lambda.ln() - self.lambda * x
    }

    fn logcdf(&self, x: f64) -> f64 {
        (-(-self.lambda * x).exp()).ln_1p()
    }

    fn logccdf(&self, x: f64) -> f64 {
        -self.lambda * x
    }

    fn support(&self) -> (f64, f64) {
        (0.0, f64::INFINITY)
    }
}

#[derive(Debug, Clone)]
pub struct ExponentialEstimator {
    sample_mean: f64,
    degrees_of_freedom: usize,
}

impl ParameterizedDistribution for ExponentialEstimator {}

impl ExponentialEstimator {
    pub fn new(sample_mean: f64, sample_size: usize) -> Self {
        Self {
            sample_mean: sample_mean,
            degrees_of_freedom: sample_size,
        }
    }
}

impl From<ExponentialEstimator> for Exponential {
    fn from(value: ExponentialEstimator) -> Self {
        Self::from_scale(value.sample_mean)
    }
}

impl Default for ExponentialEstimator {
    fn default() -> Self {
        Self {
            sample_mean: 1.0,
            degrees_of_freedom: 1,
        }
    }
}

impl Distribution for ExponentialEstimator {
    fn logpdf(&self, x: f64) -> f64 {
        let n = self.degrees_of_freedom as f64;
        let sm = self.sample_mean;
        ((n + 1.0) * n.ln() + n * sm.ln()) - ((n + 1.0) * (n * sm + x).ln())
    }

    fn pdf(&self, x: f64) -> f64 {
        self.logpdf(x).exp()
    }

    fn logccdf(&self, x: f64) -> f64 {
        let n = self.degrees_of_freedom as f64;
        let sm = self.sample_mean;
        n * ((n * sm).ln() - (n * sm + x).ln())
    }

    fn logcdf(&self, x: f64) -> f64 {
        (-self.logccdf(x).exp()).ln_1p()
    }

    fn cdf(&self, x: f64) -> f64 {
        -(self.logccdf(x).exp_m1())
    }

    fn ccdf(&self, x: f64) -> f64 {
        self.logccdf(x).exp()
    }

    fn ppf(&self, p: f64) -> f64 {
        let n = self.degrees_of_freedom as f64;
        let sm = self.sample_mean;
        (n * sm) * ((1.0 - p).powf(-1.0 / n) - 1.0)
    }

    fn support(&self) -> (f64, f64) {
        (0.0, f64::INFINITY)
    }
}

#[derive(Debug, Clone)]
pub struct HalfT {
    standard_deviation: f64,
    degrees_of_freedom: usize,
}

impl ParameterizedDistribution for HalfT {}

impl HalfT {
    #[allow(dead_code)]
    pub fn new(standard_deviation: f64, degrees_of_freedom: usize) -> Self {
        Self {
            standard_deviation,
            degrees_of_freedom,
        }
    }

    pub fn from_sample_mean(mean: f64, degrees_of_freedom: usize) -> Self {
        Self {
            standard_deviation: mean * (2.0 / f64::consts::PI).sqrt(),
            degrees_of_freedom,
        }
    }
}

impl Default for HalfT {
    fn default() -> Self {
        Self {
            standard_deviation: 1.0,
            degrees_of_freedom: 1,
        }
    }
}

impl Distribution for HalfT {
    fn logpdf(&self, x: f64) -> f64 {
        let v = self.degrees_of_freedom as f64;
        let s = self.standard_deviation;
        let z = x / s;
        if z >= 0.0 {
            let norm = (2.0_f64).ln() - (0.5 * v.ln() + beta(0.5, 0.5 * v).ln() + s.ln());
            norm - 0.5 * (v + 1.0) * ((z * z) / v).ln_1p()
        } else {
            0.0
        }
    }

    fn pdf(&self, x: f64) -> f64 {
        self.logpdf(x).exp()
    }

    fn cdf(&self, x: f64) -> f64 {
        let v = self.degrees_of_freedom as f64;
        let z = x / self.standard_deviation;
        if z >= 0.0 {
            1.0 - betai(0.5 * v, 0.5, v / (z * z + v))
        } else {
            0.0
        }
    }

    fn logcdf(&self, x: f64) -> f64 {
        self.cdf(x).ln()
    }

    fn ppf(&self, p: f64) -> f64 {
        let v = self.degrees_of_freedom as f64;
        let inv_out = invbetai(1.0 - p, 0.5 * v, 0.5);
        let x_unit = (v / inv_out - v).sqrt();
        x_unit * self.standard_deviation
    }

    fn ccdf(&self, x: f64) -> f64 {
        1.0 - self.cdf(x)
    }

    fn logccdf(&self, x: f64) -> f64 {
        self.ccdf(x).ln()
    }

    fn support(&self) -> (f64, f64) {
        (0.0, f64::INFINITY)
    }
}

#[derive(Debug, Clone)]
pub struct Lomax {
    alpha: f64,
    lambda: f64,
}

impl Lomax {
    pub fn new(alpha: f64, lambda: f64) -> Self {
        Self { alpha, lambda }
    }
}

impl ParameterizedDistribution for Lomax {}

impl Default for Lomax {
    fn default() -> Self {
        Self::new(1.0, 1.0)
    }
}

impl Distribution for Lomax {
    fn logpdf(&self, x: f64) -> f64 {
        let a = self.alpha;
        let y = self.lambda;
        (a / y).ln() - (a + 1.0) * (1.0 + x / y).ln()
    }

    fn pdf(&self, x: f64) -> f64 {
        self.logpdf(x).exp()
    }

    fn logccdf(&self, x: f64) -> f64 {
        let a = self.alpha;
        let y = self.lambda;
        -a * (1.0 + x / y).ln()
    }

    fn ccdf(&self, x: f64) -> f64 {
        self.logccdf(x).exp()
    }

    fn cdf(&self, x: f64) -> f64 {
        -(self.logccdf(x).exp_m1())
    }

    fn logcdf(&self, x: f64) -> f64 {
        (-self.ccdf(x)).ln_1p()
    }

    fn ppf(&self, p: f64) -> f64 {
        let a = self.alpha;
        let y = self.lambda;
        y * ((1.0 - p).powf(-1.0 / a) - 1.0)
    }

    fn support(&self) -> (f64, f64) {
        (0.0, f64::INFINITY)
    }
}

#[derive(Debug, Clone)]
pub struct AssymetricLaplace {
    mode: f64,
    scale: f64,
    mode_quantile: f64,
}

impl ParameterizedDistribution for AssymetricLaplace {}

impl AssymetricLaplace {
    pub fn new(mode: f64, scale: f64, mode_quantile: f64) -> Self {
        Self {
            mode,
            scale,
            mode_quantile,
        }
    }

    pub fn from_exponential_halves(mode: f64, negative_mean: f64, positive_mean: f64) -> Self {
        let nm = negative_mean.max(1e-8);
        let pm = positive_mean.max(1e-8);
        Self::new(mode, (nm * pm) / (nm + pm), 1.0 / (pm / nm + 1.0))
    }

    pub fn symmetric_from_moments(mean: f64, standard_deviation: f64) -> Self {
        Self::new(mean, standard_deviation / (8.0_f64.sqrt()), 0.5)
    }
}

impl Default for AssymetricLaplace {
    fn default() -> Self {
        Self::symmetric_from_moments(0.0, 1.0)
    }
}

impl Distribution for AssymetricLaplace {
    fn logpdf(&self, x: f64) -> f64 {
        let m = self.mode;
        let l = self.scale;
        let p = self.mode_quantile;
        let exp_comp = if x <= m {
            ((1.0 - p) / l) * (x - m)
        } else {
            -(p / l) * (x - m)
        };

        ((p * (1.0 - p)) / l).ln() + exp_comp
    }

    fn pdf(&self, x: f64) -> f64 {
        self.logpdf(x).exp()
    }

    fn cdf(&self, x: f64) -> f64 {
        let m = self.mode;
        let l = self.scale;
        let p = self.mode_quantile;
        if x <= m {
            p * (((1.0 - p) / l) * (x - m)).exp()
        } else {
            1.0 - (1.0 - p) * (-(p / l) * (x - m)).exp()
        }
    }

    fn logcdf(&self, x: f64) -> f64 {
        self.cdf(x).ln()
    }

    fn ccdf(&self, x: f64) -> f64 {
        1.0 - self.cdf(x)
    }

    fn logccdf(&self, x: f64) -> f64 {
        self.ccdf(x).ln()
    }

    fn ppf(&self, p: f64) -> f64 {
        let m = self.mode;
        let l = self.scale;
        let pm = self.mode_quantile;
        if p <= pm {
            m + (l / (1.0 - pm)) * (p / pm).ln()
        } else {
            m - (l / pm) * ((1.0 - p) / (1.0 - pm)).ln()
        }
    }

    fn support(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use std::fmt::Debug;

    pub fn linspace(start: f64, stop: f64, steps: usize) -> impl Iterator<Item = f64> {
        (0..steps)
            .map(move |n| n as f64 / (steps as f64 - 1.0))
            .map(move |n| start * (1.0 - n) + stop * n)
    }

    // Add debug trait to allow for printout...
    pub trait TestDistribution: Distribution + Debug {}
    impl<T: Distribution + Debug> TestDistribution for T {}

    fn as_box<T: TestDistribution + 'static>(d: T) -> Box<dyn TestDistribution> {
        Box::new(d)
    }

    use super::{Exponential, ParameterizedDistribution};

    fn get_dists() -> [Box<dyn TestDistribution>; 5] {
        [
            as_box(Exponential::unit()),
            as_box(ExponentialEstimator::unit()),
            as_box(HalfT::unit()),
            as_box(AssymetricLaplace::unit()),
            as_box(Lomax::unit()),
        ]
    }

    fn is_close(a: f64, b: f64) -> bool {
        let rel_tol = 1e-9;
        let abs_tol = 0.0;
        (a - b).abs() <= (rel_tol * (a.abs()).max(b.abs())).max(abs_tol)
    }

    #[test]
    fn basic_distribution_propery_checks() {
        for dist in get_dists() {
            println!("Testing distribution: {:?}", dist);
            let (mut low, mut high) = dist.support();

            assert!(dist.cdf(low) == 0.0);
            assert!(dist.cdf(high) == 1.0);

            if high == f64::INFINITY {
                high = 5.0;
            }
            if low == f64::NEG_INFINITY {
                low = -5.0;
            }

            for x in linspace(low, high, 100) {
                // Basic properties...
                // println!("{x} -> {} vs {}", dist.tpdf(x), dist.tlogpdf(x).exp());
                assert!(is_close(dist.pdf(x), dist.logpdf(x).exp()));
                assert!(is_close(dist.cdf(x), dist.logcdf(x).exp()));
                assert!(is_close(dist.ccdf(x), dist.logccdf(x).exp()));
                assert!(is_close(dist.ccdf(x), 1.0 - dist.cdf(x)));
                // println!("{x} -> {}", dist.tppf(dist.tcdf(x)));
                assert!(is_close(dist.ppf(dist.cdf(x)), x));
            }
        }
    }
}
