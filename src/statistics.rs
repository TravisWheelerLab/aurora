#[allow(dead_code)]
pub trait Distribution: Clone {
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
        (1.0 - self.ccdf(x)).ln()
    }
}

#[derive(Clone)]
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
