use super::ImplicitODE;
use crate::algebra::{abstr::Real, linear::Vector};

pub trait ImplicitFixedStepSizeMethod<T>
    where T: Real
{
    fn do_step<F>(self: &Self, prob: &F, t_n: &T, x_n: &Vector<T>, h: &T) -> Vector<T>
        where F: ImplicitODE<T>;

    fn order(self: &Self) -> u8;
}
