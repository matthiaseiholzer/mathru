// Beta distribution
use crate::{algebra::abstr::Real, special, special::beta, statistics::distrib::Continuous};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::clone::Clone;

/// Beta distribution
///
/// Fore more information:
/// <a href="https://en.wikipedia.org/wiki/Beta_distribution">https://en.wikipedia.org/wiki/Beta_distribution</a>
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Copy, Debug)]
pub struct Beta<T>
{
    p: T,
    q: T,
}

impl<T> Beta<T> where T: Real
{
    /// Create a probability distribution
    ///
    /// # Arguments
    ///
    /// * `p`: &alpha; > 0.0
    /// * `q`: &beta; > 0.0
    ///
    /// # Panics
    /// if p <= 0.0 || q <= 0.0
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Beta, Continuous};
    ///
    /// let distrib: Beta<f64> = Beta::new(0.2, 0.3);
    /// ```
    pub fn new(p: T, q: T) -> Beta<T>
    {
        if p < T::zero() || q <= T::zero()
        {
            panic!()
        }
        Beta { p, q }
    }
}

impl<T> Continuous<T> for Beta<T> where T: Real + beta::Beta
{
    /// Probability density function
    ///
    /// # Arguments
    ///
    /// * `x` &isin &#2115;
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Beta, Continuous};
    ///
    /// let distrib: Beta<f64> = Beta::new(0.2, 0.3);
    /// let x: f64 = 0.5;
    /// let p: f64 = distrib.pdf(x);
    /// ```
    fn pdf(self: &Self, x: T) -> T
    {
        if x < T::zero()
        {
            panic!();
        }
        if x > T::one()
        {
            return T::one();
        }
        return x.pow(self.p - T::one()) * (T::one() - x).pow(self.q - T::one())
               / special::beta::beta(self.p, self.q);
    }

    /// Cumulative distribution function
    ///
    /// # Arguments
    ///
    /// * `x`
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Beta, Continuous};
    ///
    /// let distrib: Beta<f64> = Beta::new(0.3, 0.2);
    /// let x: f64 = 0.4;
    /// let p: f64 = distrib.cdf(x);
    /// ```
    fn cdf(self: &Self, x: T) -> T
    {
        beta::beta_inc_reg(x, self.p, self.q)
    }

    /// Quantile function of inverse cdf
    fn quantile(self: &Self, _p: T) -> T
    {
        unimplemented!();
    }

    /// Expected value
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Beta, Continuous};
    ///
    /// let distrib: Beta<f64> = Beta::new(0.2, 0.3);
    /// let mean: f64 = distrib.mean();
    /// ```
    fn mean(self: &Self) -> T
    {
        self.p / (self.p + self.q)
    }

    /// Variance
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Beta, Continuous};
    ///
    /// let distrib: Beta<f64> = Beta::new(0.2, 0.3);
    /// let var: f64 = distrib.variance();
    /// ```
    fn variance(self: &Self) -> T
    {
        self.p * self.q / ((self.p + self.q + T::one()) * (self.p + self.q).pow(T::from_f64(2.0)))
    }

    /// Skewness
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Beta, Continuous};
    ///
    /// let distrib: Beta<f64> = Beta::new(0.2, 0.3);
    /// let skewness: f64 = distrib.skewness();
    /// ```
    fn skewness(self: &Self) -> T
    {
        return T::from_f64(2.0) * (self.q - self.p) * (self.p + self.q + T::one()).sqrt()
               / ((self.p + self.q + T::from_f64(2.0)) * (self.q * self.p).sqrt());
    }

    fn median(self: &Self) -> T
    {
        unimplemented!();
    }

    fn entropy(self: &Self) -> T
    {
        unimplemented!();
    }
}
