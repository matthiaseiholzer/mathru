//! Solves an ODE using midpoint method.
use super::{explicit_method::ExplicitMethod, ExplicitODE};
use crate::{
    algebra::{abstr::Real, linear::Vector},
    analysis::differential_equation::ordinary::ButcherFixedStepSize
};
use std::clone::Clone;

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};


/// Solves an ODE using midpoint method.
///
/// <a href="https://en.wikipedia.org/wiki/Midpoint_method">https://en.wikipedia.org/wiki/Midpoint_method</a>
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug)]
pub struct Midpoint<T>
{
    butcher: ButcherFixedStepSize<T>
}

impl<T> Default for Midpoint<T> where T: Real
{
    /// Creates a Euler instance
    fn default() -> Midpoint<T>
    {
        let a: Vec<T> = vec![T::from_f64(0.5)];
        let b: Vec<T> = vec![T::zero(), T::one()];
        let c: Vec<T> = vec![T::from_f64(0.5)];

        return Midpoint {
            butcher: ButcherFixedStepSize::new(a, b, c)
        };
    }
}

impl<T> ExplicitMethod<T> for Midpoint<T> where T: Real
{
    fn do_step<F>(self: &Self, prob: &F, t_n: &T, x_n: &Vector<T>, h: &T) -> Vector<T>
        where F: ExplicitODE<T>
    {
        return self.butcher.do_step(prob, t_n, x_n, h);
    }

    // The midpoint methods is a 2nd order method
    fn order(self: &Self) -> u8
    {
        return 2;
    }
}
