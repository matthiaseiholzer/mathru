use super::ImplicitODE;
use crate::algebra::linear::{Vector};

///
pub trait ImplicitFixedStepSizeMethod<T>
{
    fn do_step<F>(self: &Self, prob: &F, t_n: &T, x_n: &Vector<T>, h: &T) -> Vector<T>
        where F: ImplicitODE<T>;

    fn order(self: &Self) -> u8;
}

//pub trait ExplicitAdaptiveMethod<T>
//{
//    ///
//    fn do_step<F>(self: &Self, prob: &F, t_n: &T, x_n: &Vector<T>, h: &T) -> (Vector<T>, Vector<T>)
//        where F: ExplicitODE<T>;
//
//    fn order(self: &Self) -> (u8, u8);
//}