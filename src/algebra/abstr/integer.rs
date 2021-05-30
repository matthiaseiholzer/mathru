use crate::algebra::abstr::Ring;

/// Integer
///
///<a href="https://en.wikipedia.org/wiki/Integer">https://en.wikipedia.org/wiki/Integer</a>
pub trait Integer: Ring + Eq
{
}

macro_rules! impl_integer
{
    ($($t:ty),*) =>
    {
    	$(
        impl Integer for $t
        {

        }
        )*
    }
}

impl_integer!(i8, i16, i32, i64, i128);
