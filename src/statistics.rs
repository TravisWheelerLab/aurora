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
pub trait Distribution: Clone + Debug + Default {
    fn pdf(&self, x: f64) -> f64;
    fn cdf(&self, x: f64) -> f64;
    fn ppf(&self, p: f64) -> f64;
    fn support(&self) -> (f64, f64);
    fn ccdf(&self, x: f64) -> f64;
    fn logpdf(&self, x: f64) -> f64;
    fn logcdf(&self, x: f64) -> f64;
    fn logccdf(&self, x: f64) -> f64;

    fn unit() -> Self {
        Self::default()
    }
}

#[derive(Clone, Debug)]
pub struct Exponential {
    lambda: f64,
}

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
        self.cdf(x).ln()
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

#[cfg(test)]
mod test {
    use crate::statistics::{ExponentialEstimator, HalfT};
    use std::fmt::Debug;

    pub trait TestDistribution: Debug {
        fn tpdf(&self, x: f64) -> f64;
        fn tcdf(&self, x: f64) -> f64;
        fn tppf(&self, p: f64) -> f64;
        fn tsupport(&self) -> (f64, f64);
        fn tccdf(&self, x: f64) -> f64;
        fn tlogpdf(&self, x: f64) -> f64;
        fn tlogcdf(&self, x: f64) -> f64;
        fn tlogccdf(&self, x: f64) -> f64;
    }

    impl<T: Distribution> TestDistribution for T {
        fn tpdf(&self, x: f64) -> f64 {
            self.pdf(x)
        }
        fn tcdf(&self, x: f64) -> f64 {
            self.cdf(x)
        }
        fn tppf(&self, p: f64) -> f64 {
            self.ppf(p)
        }
        fn tsupport(&self) -> (f64, f64) {
            self.support()
        }
        fn tccdf(&self, x: f64) -> f64 {
            self.ccdf(x)
        }
        fn tlogpdf(&self, x: f64) -> f64 {
            self.logpdf(x)
        }
        fn tlogcdf(&self, x: f64) -> f64 {
            self.logcdf(x)
        }
        fn tlogccdf(&self, x: f64) -> f64 {
            self.logccdf(x)
        }
    }

    fn as_box<T: Distribution + 'static>(d: T) -> Box<dyn TestDistribution> {
        Box::new(d)
    }

    use super::{Distribution, Exponential};

    fn get_dists() -> [Box<dyn TestDistribution>; 3] {
        [
            as_box(Exponential::unit()),
            as_box(ExponentialEstimator::unit()),
            as_box(HalfT::unit()),
        ]
    }

    fn is_close(a: f64, b: f64) -> bool {
        let rel_tol = 1e-9;
        let abs_tol = 0.0;
        (a - b).abs() <= (rel_tol * (a.abs()).max(b.abs())).max(abs_tol)
    }

    fn linspace(start: f64, stop: f64, steps: usize) -> impl Iterator<Item = f64> {
        (0..steps)
            .map(move |n| n as f64 / (steps as f64 - 1.0))
            .map(move |n| start * (1.0 - n) + stop * n)
    }

    #[test]
    fn basic_distribution_propery_checks() {
        for dist in get_dists() {
            println!("Testing distribution: {:?}", dist);
            let (mut low, mut high) = dist.tsupport();

            if high == f64::INFINITY {
                high = 5.0;
            }
            if low == f64::INFINITY {
                low = -5.0;
            }

            for x in linspace(low, high, 100) {
                // Basic properties...
                assert!(is_close(dist.tpdf(x), dist.tlogpdf(x).exp()));
                assert!(is_close(dist.tcdf(x), dist.tlogcdf(x).exp()));
                assert!(is_close(dist.tccdf(x), dist.tlogccdf(x).exp()));
                assert!(is_close(dist.tccdf(x), 1.0 - dist.tcdf(x)));
                assert!(is_close(dist.tppf(dist.tcdf(x)), x));
            }
        }
    }

    #[test]
    fn test_exponential_distribution() {
        let _dist = Exponential::unit();
    }
}
