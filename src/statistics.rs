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
pub trait Distribution: Clone {
    fn pdf(&self, x: f64) -> f64;
    fn cdf(&self, x: f64) -> f64;
    fn ppf(&self, p: f64) -> f64;
    fn support(&self) -> (f64, f64);
    fn ccdf(&self, x: f64) -> f64;
    fn logpdf(&self, x: f64) -> f64;
    fn logcdf(&self, x: f64) -> f64;
    fn logccdf(&self, x: f64) -> f64;
}

pub trait ParameterizedDistribution: Distribution + Debug + Default {
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
pub struct Frechet {
    alpha: f64,
    scale: f64,
    minimum: f64,
}

impl ParameterizedDistribution for Frechet {}

impl Frechet {
    pub fn new(alpha: f64, scale: f64, minimum: f64) -> Self {
        Self {
            alpha,
            scale,
            minimum,
        }
    }
}

impl Default for Frechet {
    fn default() -> Self {
        Self {
            alpha: 1.0,
            scale: 1.0,
            minimum: 0.0,
        }
    }
}

impl Distribution for Frechet {
    fn logpdf(&self, x: f64) -> f64 {
        let a = self.alpha;
        let s = self.scale;
        let m = self.minimum;
        if x > m {
            (a / s).ln() + -(a + 1.0) * ((x - m) / s).ln() + -((x - m) / s).powf(-a)
        } else {
            f64::NEG_INFINITY
        }
    }

    fn pdf(&self, x: f64) -> f64 {
        self.logpdf(x).exp()
    }

    fn cdf(&self, x: f64) -> f64 {
        self.logcdf(x).exp()
    }

    fn logcdf(&self, x: f64) -> f64 {
        let a = self.alpha;
        let s = self.scale;
        let m = self.minimum;
        if x > m {
            -((x - m) / s).powf(-a)
        } else {
            f64::NEG_INFINITY
        }
    }

    fn ppf(&self, p: f64) -> f64 {
        let a = self.alpha;
        let s = self.scale;
        let m = self.minimum;
        if p >= 1.0 {
            f64::INFINITY
        } else if p <= 0.0 {
            m
        } else {
            m + s * (-p.min(1.0).ln()).powf(1.0 / -a)
        }
    }

    fn ccdf(&self, x: f64) -> f64 {
        1.0 - self.cdf(x)
    }

    fn logccdf(&self, x: f64) -> f64 {
        self.ccdf(x).ln()
    }

    fn support(&self) -> (f64, f64) {
        (self.minimum, f64::INFINITY)
    }
}

#[derive(Debug, Clone)]
pub struct Gumbel {
    location: f64,
    scale: f64,
}

impl Gumbel {
    pub fn new(location: f64, scale: f64) -> Self {
        Self { location, scale }
    }
}

impl Default for Gumbel {
    fn default() -> Self {
        Self::new(0.0, 1.0)
    }
}

impl ParameterizedDistribution for Gumbel {}

impl Distribution for Gumbel {
    fn logpdf(&self, x: f64) -> f64 {
        let mu = self.location;
        let beta = self.scale;
        let z = (x - mu) / beta;
        (1.0 / beta).ln() - (z + (-z).exp())
    }

    fn pdf(&self, x: f64) -> f64 {
        self.logpdf(x).exp()
    }

    fn cdf(&self, x: f64) -> f64 {
        self.logcdf(x).exp()
    }

    fn logcdf(&self, x: f64) -> f64 {
        let mu = self.location;
        let beta = self.scale;
        let z = (x - mu) / beta;
        -((-z).exp())
    }

    fn ppf(&self, p: f64) -> f64 {
        let mu = self.location;
        let beta = self.scale;
        mu - beta * (-p.ln()).ln()
    }

    fn ccdf(&self, x: f64) -> f64 {
        1.0 - self.cdf(x)
    }

    fn logccdf(&self, x: f64) -> f64 {
        self.ccdf(x).ln()
    }

    fn support(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
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
pub struct Laplace {
    mean: f64,
    scale: f64,
}

impl ParameterizedDistribution for Laplace {}

impl Laplace {
    pub fn new(mean: f64, scale: f64) -> Self {
        Self { mean, scale }
    }

    pub fn from_moments(mean: f64, standard_deviation: f64) -> Self {
        Self {
            mean,
            scale: standard_deviation / f64::consts::SQRT_2,
        }
    }
}

impl Default for Laplace {
    fn default() -> Self {
        Self {
            mean: 0.0,
            scale: 1.0,
        }
    }
}

impl Distribution for Laplace {
    fn logpdf(&self, x: f64) -> f64 {
        let mu = self.mean;
        let b = self.scale;
        (0.5 / b).ln() + -((x - mu).abs() / b)
    }

    fn pdf(&self, x: f64) -> f64 {
        self.logpdf(x).exp()
    }

    fn cdf(&self, x: f64) -> f64 {
        let mu = self.mean;
        let b = self.scale;
        0.5 + 0.5 * (x - mu).signum() * (1.0 - (-(x - mu).abs() / b).exp())
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
        let mu = self.mean;
        let b = self.scale;
        let p = p.clamp(0.0, 1.0);
        mu - b * (p - 0.5).signum() * (1.0 - 2.0 * (p - 0.5).abs()).ln()
    }

    fn support(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }
}

pub fn linspace(start: f64, stop: f64, steps: usize) -> impl Iterator<Item = f64> {
    (0..steps)
        .map(move |n| n as f64 / (steps as f64 - 1.0))
        .map(move |n| start * (1.0 - n) + stop * n)
}

#[cfg(test)]
mod test {
    use super::*;
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

    impl<T: ParameterizedDistribution> TestDistribution for T {
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

    fn as_box<T: ParameterizedDistribution + 'static>(d: T) -> Box<dyn TestDistribution> {
        Box::new(d)
    }

    use super::{Exponential, ParameterizedDistribution};

    fn get_dists() -> [Box<dyn TestDistribution>; 7] {
        [
            as_box(Exponential::unit()),
            as_box(ExponentialEstimator::unit()),
            as_box(HalfT::unit()),
            as_box(Frechet::unit()),
            as_box(Laplace::unit()),
            as_box(Gumbel::unit()),
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
            let (mut low, mut high) = dist.tsupport();

            assert!(dist.tcdf(low) == 0.0);
            assert!(dist.tcdf(high) == 1.0);

            if high == f64::INFINITY {
                high = 5.0;
            }
            if low == f64::NEG_INFINITY {
                low = -5.0;
            }

            for x in linspace(low, high, 100) {
                // Basic properties...
                // println!("{x} -> {} vs {}", dist.tpdf(x), dist.tlogpdf(x).exp());
                assert!(is_close(dist.tpdf(x), dist.tlogpdf(x).exp()));
                assert!(is_close(dist.tcdf(x), dist.tlogcdf(x).exp()));
                assert!(is_close(dist.tccdf(x), dist.tlogccdf(x).exp()));
                assert!(is_close(dist.tccdf(x), 1.0 - dist.tcdf(x)));
                // println!("{x} -> {}", dist.tppf(dist.tcdf(x)));
                assert!(is_close(dist.tppf(dist.tcdf(x)), x));
            }
        }
    }
}
