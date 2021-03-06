//! Poisson distribution
use crate::{
    algebra::abstr::Real,
    special,
    special::gamma::Gamma,
    statistics::{combins, distrib::Discrete},
};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::clone::Clone;

/// Poisson distribution
///
/// Fore more information:
/// <a href="https://en.wikipedia.org/wiki/Poisson_distribution">https://en.wikipedia.org/wiki/Poisson_distribution</a>
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Copy, Debug)]
pub struct Poisson<T>
{
    gamma: T,
}

impl<T> Poisson<T> where T: Real
{
    /// Creates a probability distribution
    ///
    /// # Arguments
    ///
    /// * `gamma` gamma > 0.0
    ///
    /// # Panics
    ///
    /// if gamma <= 0.0
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::Poisson;
    ///
    /// let distrib: Poisson<f64> = Poisson::new(&0.2);
    /// ```
    pub fn new(gamma: &T) -> Poisson<T>
    {
        if *gamma <= T::zero()
        {
            panic!();
        }

        Poisson { gamma: *gamma }
    }
}

impl<T> Discrete<T, u32, u32> for Poisson<T>
    where T: Real + Gamma
{
    /// Probability mass function
    ///
    /// # Arguments
    ///
    /// * `x` Random variable x &isin; &#x2115;
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Discrete, Poisson};
    ///
    /// let distrib: Poisson<f64> = Poisson::new(&0.2);
    /// let x: u32 = 5;
    /// let p: f64 = distrib.pmf(x);
    /// ```
    fn pmf<'a>(self: &'a Self, x: u32) -> T
    {
        let k_fact: T = T::from_u64(combins::factorial(x));
        self.gamma.pow(T::from_u32(x)) * (-self.gamma).exp() / k_fact
    }

    /// Cumulative distribution function of the Bernoulli distribution
    ///
    /// # Arguments
    ///
    /// * `x` Random variable x &isin; &#x2115;
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Discrete, Poisson};
    ///
    /// let distrib: Poisson<f64> = Poisson::new(&0.2);
    /// let x: u32 = 4;
    /// let p: f64 = distrib.cdf(x);
    /// ```
    fn cdf<'a>(self: &'a Self, x: u32) -> T
    {
        return special::gamma::gamma_ur(T::from_u32(x + 1), self.gamma);
    }

    /// Expected value
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Discrete, Poisson};
    ///
    /// let distrib: Poisson<f64> = Poisson::new(&0.2);
    /// let mean: f64 = distrib.mean();
    /// ```
    fn mean<'a>(self: &'a Self) -> T
    {
        return self.gamma;
    }

    /// Variance
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::statistics::distrib::{Discrete, Poisson};
    ///
    /// let distrib: Poisson<f64> = Poisson::new(&0.2);
    /// let var: f64 = distrib.variance();
    /// ```
    fn variance<'a>(self: &'a Self) -> T
    {
        return self.gamma;
    }
}
