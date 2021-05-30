/// Power functions
///
///<a href="https://en.wikipedia.org/wiki/Exponentiation#Power_functions">https://en.wikipedia.org/wiki/Exponentiation#Power_functions</a>
pub trait Power
{
    fn pow(self: Self, exp: Self) -> Self;

    fn root(self: Self, root: Self) -> Self;

    fn sqrt(self: Self) -> Self;
}

macro_rules! power_impl {
    ($t:ty) => {
        impl Power for $t
        {
            fn pow(self: Self, exp: Self) -> Self
            {
                return self.powf(exp);
            }

            fn root(self: Self, root: Self) -> Self
            {
                return self.powf(1.0 / root);
            }

            fn sqrt(self: Self) -> Self
            {
                return self.powf(0.5);
            }
        }
    };
}

power_impl!(f32);
power_impl!(f64);

macro_rules! power_impl_integer {
    ($t:ty) => {
        impl Power for $t
        {
            fn pow(self: Self, _exp: Self) -> Self
            {
                unimplemented!();
            }

            fn root(self: Self, _root: Self) -> Self
            {
                unimplemented!();
            }

            fn sqrt(self: Self) -> Self
            {
                unimplemented!();
            }
        }
    };
}

power_impl_integer!(u8);
power_impl_integer!(u16);
power_impl_integer!(u32);
power_impl_integer!(u64);
power_impl_integer!(u128);
power_impl_integer!(usize);
power_impl_integer!(i8);
power_impl_integer!(i16);
power_impl_integer!(i32);
power_impl_integer!(i64);
power_impl_integer!(i128);
power_impl_integer!(isize);