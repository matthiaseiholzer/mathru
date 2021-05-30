//! Magma
use super::operator::{Addition, Multiplication, Operator};
use std::{
    ops::{Add, AddAssign, Mul, MulAssign},
};


/// A Magma is a pair $`(\mathbb{M}, \circ)`$, composed by a set $`\mathbb{M}`$
/// and a binary inner operation $`\circ`$:
///
/// ```math
/// \circ: \mathbb{M} \times \mathbb{M} \rightarrow \mathbb{M} ,\\\\ (x, y) \mapsto x \circ y
/// ```
pub trait Magma<O: Operator>: Clone + PartialEq
{
    /// binary operation
    fn operate(self, rhs: Self) -> Self;
}


macro_rules! impl_magma
(
    ($O:ty; $op: ident; $($T:ty),*) =>
    {
        $(
            impl Magma<$O> for $T
            {
                fn operate(self, rhs: Self) -> Self
                {
                    self.$op(rhs)
                }
            }
        )*
    }
);

impl_magma!(Addition; add; u8, u16, u32, u64, u128, i8, i16, i32, i64, i128, f32, f64);
impl_magma!(Multiplication; mul; u8, u16, u32, u64, u128, i8, i16, i32, i64, i128, f32, f64);

/// Syntactic sugar for Magma::<Addition>::operate(a, b)
pub trait MagmaAdd: Magma<Addition> + Add<Self, Output = Self> + AddAssign<Self>
{
}

macro_rules! impl_magmaadd
(
    ($($T:ty),*) =>
    {
        $(
            impl MagmaAdd for $T
            {

            }
        )*
    }
);

impl_magmaadd!(u8, u16, u32, u64, u128, i8, i16, i32, i64, i128, f32, f64);

/// Syntactic sugar for Magma::<Multiplication>::operate(a, b)
pub trait MagmaMul: Magma<Multiplication> + Mul<Self, Output = Self> + MulAssign<Self>
{
}

macro_rules! impl_magmamul
(
    ($($T:ty),*) =>
    {
        $(
            impl MagmaMul for $T
            {

            }
        )*
    }
);

impl_magmamul!(u8, u16, u32, u64, u128, i8, i16, i32, i64, i128, f32, f64);
