/// Exponential function and its inverse
///
///<a href="https://en.wikipedia.org/wiki/Exponential_function">https://en.wikipedia.org/wiki/Exponential_function</a>
///

use std::{f32, f64};

pub trait Exponential
{
    /// Euler's number
    fn e() -> Self;

    ///Exponential function
    fn exp(self: Self) -> Self;

    /// Natural logiarithm function
    fn ln(self: Self) -> Self;
}

macro_rules! exponential_impl {
    ($t:ty, $e: expr) => {
        impl Exponential for $t
        {
            ///
            fn e() -> Self
            {
                $e
            }

            ///Exponential function
            fn exp(self: Self) -> Self
            {
                self.exp()
            }

            ///Logarithm function
            fn ln(self: Self) -> Self
            {
                self.ln()
            }
        }
    };
}

exponential_impl!(f32, f32::consts::E);
exponential_impl!(f64, f64::consts::E);