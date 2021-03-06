use super::ExplicitODE;
use crate::algebra::linear::Vector;

///
pub trait ExplicitMethod<T>
{
    ///
    fn do_step<F>(self: &Self, prob: &F, t_n: &T, x_n: &Vector<T>, h: &T) -> Vector<T>
        where F: ExplicitODE<T>;

    fn order(self: &Self) -> u8;
}

pub trait ExplicitEmbeddedMethod<T>
{
    ///
    fn do_step<F>(self: &Self, prob: &F, t_n: &T, x_n: &Vector<T>, h: &T) -> (Vector<T>, Vector<T>)
        where F: ExplicitODE<T>;

    fn order(self: &Self) -> (u8, u8);
}
