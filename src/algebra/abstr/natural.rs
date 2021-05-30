//use super::Semiring;

/// Natural number
///
///<a href="https://en.wikipedia.org/wiki/Natural_number">https://en.wikipedia.org/wiki/Natural_number</a>
pub trait Natural: Eq
{
}

macro_rules! impl_natural
{
	($($t:ty),*) =>
    {
    	$(
        impl Natural for $t
        {

        }
        )*
    }
}

impl_natural!(u8, u16, u32, u64, u128);
