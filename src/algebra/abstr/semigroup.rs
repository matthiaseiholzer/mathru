//! Semigroup
use super::{
    magma::{Magma, MagmaAdd, MagmaMul},
    operator::{Addition, Multiplication, Operator},
};

/// A Semigroup is a pair $`(\mathbb{S}, \circ)`$, composed by a set
/// $`\mathbb{S}`$ and a binary inner operation $`\circ`$: # Definition
///
/// ```math
/// \circ: \mathbb{S} \times \mathbb{S} \rightarrow \mathbb{S} , (x, y) \mapsto x \circ y
/// ```
/// and is associative
/// $`x, y, z \in \mathbb{S}`$
/// $`x \circ (y \circ z) = (x \circ y) \circ z`$
pub trait Semigroup<O: Operator + Copy>: Magma<O>
{
    fn is_associative(self: Self, y: Self, z: Self) -> bool
    {
        return self.clone().operate(y.clone()).operate(z.clone()) == self.operate(y.operate(z));
    }
}

macro_rules! impl_semigroup
(
    ($O:ty; $op: ident; $($T:ty),*) =>
    {
        $(
            impl Semigroup<$O> for $T
            {
            }
        )*
    }
);

impl_semigroup!(Addition; add; u8, u16, u32, u64, u128, i8, i16, i32, i64, i128, f32, f64);
impl_semigroup!(Multiplication; mul; u8, u16, u32, u64, u128, i8, i16, i32, i64, i128, f32, f64);

pub trait SemigroupAdd: Semigroup<Addition> + MagmaAdd
{
}

macro_rules! impl_semigroupadd
(
    ($($T:ty),*) =>
    {
        $(
            impl SemigroupAdd for $T
            {

            }
        )*
    }
);

impl_semigroupadd!(u8, u16, u32, u64, u128, i8, i16, i32, i64, i128, f32, f64);

pub trait SemigroupMul: Semigroup<Multiplication> + MagmaMul
{
}

macro_rules! impl_semigroupmul
(
    ($($T:ty),*) =>
    {
        $(
            impl SemigroupMul for $T
            {

            }
        )*
    }
);

impl_semigroupmul!(u8, u16, u32, u64, u128, i8, i16, i32, i64, i128, f32, f64);
