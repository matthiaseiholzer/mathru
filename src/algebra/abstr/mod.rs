//! Abstract algebra
//!
//! <a href="https://en.wikipedia.org/wiki/Abstract_algebra">https://en.wikipedia.org/wiki/Abstract_algebra</a>

pub use self::{
    abeliangroup::{AbelianGroup, AbelianGroupAdd, AbelianGroupMul},
    abs_diff_eq::{AbsDiff, AbsDiffEq},
    group::{Group, GroupAdd, GroupMul},
    identity::Identity,
    loop_::Loop,
    magma::{Magma, MagmaAdd, MagmaMul},
    monoid::{Monoid, MonoidAdd, MonoidMul, One, Zero},
    operator::{Addition, Multiplication, Operator},
    quasigroup::Quasigroup,
    relative_eq::{Relative, RelativeEq},
    semigroup::{Semigroup, SemigroupAdd, SemigroupMul},
};
pub use self::{
    complex::Complex,
    field::Field,
    integer::Integer,
    natural::Natural,
    real::Real,
    ring::{CommutativeRing, Ring},
    scalar::Scalar,
    sign::Sign,
};
pub use self::polynomial::Polynomial;
#[cfg(feature = "lapack")]
pub use self::scalar::{Blas, Lapack};

pub mod cast;
mod field;
mod ring;
mod scalar;
mod sign;
mod semiring;
mod abeliangroup;
// mod complex;
mod group;
mod identity;
mod integer;
mod loop_;
mod magma;
mod monoid;
mod natural;
mod operator;
mod quasigroup;
mod real;
mod semigroup;
mod abs_diff_eq;
mod relative_eq;

mod polynomial;
#[macro_use]
//pub mod real;
//pub mod integer;
//pub mod natural;
pub mod complex;
