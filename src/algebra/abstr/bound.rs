use std::{f32, f64, i128, i16, i32, i64, i8, u128, u16, u32, u64, u8};

/// Numbers which have upper and lower bounds
pub trait Bound
{
    /// returns the smallest finite number this type can represent
    fn lower_bound() -> Self;
    /// returns the largest finite number this type can represent
    fn upper_bound() -> Self;
}

macro_rules! impl_bound {
    ($t:ty, $min:expr, $max:expr) => {
        impl Bound for $t
        {
            fn lower_bound() -> $t
            {
                $min
            }

            fn upper_bound() -> $t
            {
                $max
            }
        }
    };
}

impl_bound!(u8, u8::MIN, u8::MAX);
impl_bound!(u16, u16::MIN, u16::MAX);
impl_bound!(u32, u32::MIN, u32::MAX);
impl_bound!(u64, u64::MIN, u64::MAX);
impl_bound!(u128, u128::MIN, u128::MAX);
impl_bound!(i8, i8::MIN, i8::MAX);
impl_bound!(i16, i16::MIN, i16::MAX);
impl_bound!(i32, i32::MIN, i32::MAX);
impl_bound!(i64, i64::MIN, i64::MAX);
impl_bound!(i128, i128::MIN, i128::MAX);
impl_bound!(f32, f32::MIN, f32::MAX);
impl_bound!(f64, f64::MIN, f64::MAX);
