//! quasigroup

use super::{Addition, Identity, Magma, Multiplication, Operator};

/// A quasigroup is a magma which has the divisibility property (or Latin square
/// property). Divisibility is a weak form of right and left invertibility.
///
/// $`\forall a, b \in \mathbb{Q}, \exists! r, l \in \mathbb{Q}`$ such that $`l
/// ∘ a = b`$  and $`a ∘ r = b`$
///
/// The solutions for $`r`$ and $`l`$ are:
///
/// $`r = a \ b`$ and $`l = b / a`$
///
/// where $`\`$ is the left and $`/`$ is th right division.
pub trait Quasigroup<O: Operator>: Magma<O> + Identity<O> + PartialEq
{

}

macro_rules! impl_quasigroup(
    ($T:ty, $($S:ty),*) =>
    {
        $(
        impl Quasigroup<$T> for $S
        {

        }
        )*
    }
);

impl_quasigroup!(Addition, i8, i16, i32, i64, i128, f32, f64);
impl_quasigroup!(Multiplication, f32, f64);
