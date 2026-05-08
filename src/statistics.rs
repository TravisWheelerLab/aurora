use puruspe::{beta, betai, invbetai};
use std::{f64, fmt::Debug};

#[allow(dead_code)]
pub trait Distribution: Clone + Debug {
    fn unit() -> Self;
    fn pdf(&self, x: f64) -> f64;
    fn cdf(&self, x: f64) -> f64;
    fn ppf(&self, p: f64) -> f64;
    fn support(&self) -> (f64, f64);

    fn ccdf(&self, x: f64) -> f64 {
        1.0 - self.cdf(x)
    }
    fn logpdf(&self, x: f64) -> f64 {
        self.pdf(x).ln()
    }
    fn logcdf(&self, x: f64) -> f64 {
        self.cdf(x).ln()
    }
    fn logccdf(&self, x: f64) -> f64 {
        self.ccdf(x).ln()
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

impl Distribution for Exponential {
    fn unit() -> Self {
        Self::new(1.0)
    }

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

impl Distribution for ExponentialEstimator {
    fn unit() -> Self {
        Self {
            sample_mean: 1.0,
            degrees_of_freedom: 1,
        }
    }

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
pub struct StudentsT {
    mean: f64,
    standard_deviation: f64,
    degrees_of_freedom: usize,
}

impl StudentsT {
    pub fn new(mean: f64, standard_deviation: f64, degrees_of_freedom: usize) -> Self {
        Self {
            mean,
            standard_deviation,
            degrees_of_freedom,
        }
    }
}

impl Distribution for StudentsT {
    fn unit() -> Self {
        Self {
            mean: 0.0,
            standard_deviation: 1.0,
            degrees_of_freedom: 1,
        }
    }

    fn pdf(&self, x: f64) -> f64 {
        let v = self.degrees_of_freedom as f64;
        let s = self.standard_deviation;
        let z = (x - self.mean) / s;
        1.0 / (v.sqrt() * beta(0.5, 0.5 * v) * s) * (1.0 + (z * z) / v).powf(-0.5 * (v + 1.0))
    }

    fn cdf(&self, x: f64) -> f64 {
        let v = self.degrees_of_freedom as f64;
        let z = (x - self.mean) / self.standard_deviation;
        let beta_comp = betai(0.5 * v, 0.5, v / (z * z + v));
        if z > 0.0 {
            1.0 - 0.5 * beta_comp
        } else {
            0.5 * beta_comp
        }
    }

    fn ppf(&self, p: f64) -> f64 {
        let v = self.degrees_of_freedom as f64;
        let p_in = if p <= 0.5 { 2.0 * p } else { 2.0 * (1.0 - p) };
        let inv_out = invbetai(p_in, 0.5 * v, 0.5);
        let x_unit = (v / inv_out - v).sqrt();
        x_unit * self.standard_deviation + self.mean
    }

    fn support(&self) -> (f64, f64) {
        (f64::NEG_INFINITY, f64::INFINITY)
    }
}

#[derive(Debug, Clone)]
pub struct HalfT {
    standard_deviation: f64,
    degrees_of_freedom: usize,
}

impl HalfT {
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

impl Distribution for HalfT {
    fn unit() -> Self {
        Self {
            standard_deviation: 1.0,
            degrees_of_freedom: 1,
        }
    }

    fn pdf(&self, x: f64) -> f64 {
        let v = self.degrees_of_freedom as f64;
        let s = self.standard_deviation;
        let z = x / s;
        2.0 / (v.sqrt() * beta(0.5, 0.5 * v) * s) * (1.0 + (z * z) / v).powf(-0.5 * (v + 1.0))
    }

    fn cdf(&self, x: f64) -> f64 {
        let v = self.degrees_of_freedom as f64;
        let z = x / self.standard_deviation;
        let beta_comp = betai(0.5 * v, 0.5, v / (z * z + v));
        if z > 0.0 {
            1.0 - 0.5 * beta_comp
        } else {
            0.5 * beta_comp
        }
    }

    fn ppf(&self, p: f64) -> f64 {
        let v = self.degrees_of_freedom as f64;
        let p_in = if p <= 0.5 { 2.0 * p } else { 2.0 * (1.0 - p) };
        let inv_out = invbetai(p_in, 0.5 * v, 0.5);
        let x_unit = (v / inv_out - v).sqrt();
        x_unit * self.standard_deviation
    }

    fn support(&self) -> (f64, f64) {
        (0.0, f64::INFINITY)
    }
}
