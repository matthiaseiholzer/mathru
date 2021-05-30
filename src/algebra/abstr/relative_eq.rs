use crate::algebra::abstr::abs_diff_eq::AbsDiffEq;

/// The requisite parameters for testing for approximate equality using a
/// relative based comparison.
///
/// This is not normally used directly, rather via the
/// `assert_relative_{eq|ne}!` and `relative_{eq|ne}!` macros.
///
/// # Example
///
/// ```rust
/// use mathru::algebra::abstr::Relative;
///
/// Relative::default().eq(&1.0, &1.0);
/// Relative::default().epsilon(f64::EPSILON).eq(&1.0, &1.0);
/// Relative::default().max_relative(1.0).eq(&1.0, &1.0);
/// Relative::default().epsilon(f64::EPSILON).max_relative(1.0).eq(&1.0, &1.0);
/// Relative::default().max_relative(1.0).epsilon(f64::EPSILON).eq(&1.0, &1.0);
/// ```
pub struct Relative<A, B = A>
    where
        A: RelativeEq<B> + ?Sized,
        B: ?Sized,
{
    /// The tolerance to use when testing values that are close together.
    pub epsilon: A::Epsilon,
    /// The relative tolerance for testing values that are far-apart.
    pub max_relative: A::Epsilon,
}

impl<A, B> Default for Relative<A, B>
    where
        A: RelativeEq<B> + ?Sized,
        B: ?Sized,
{
    fn default() -> Relative<A, B>
    {
        Relative {
            epsilon: A::default_epsilon(),
            max_relative: A::default_max_relative(),
        }
    }
}

impl<A, B> Relative<A, B>
    where
        A: RelativeEq<B> + ?Sized,
        B: ?Sized,
{
    /// Replace the epsilon value with the one specified.
    pub fn epsilon(self, epsilon: A::Epsilon) -> Relative<A, B>
    {
        Relative { epsilon, ..self }
    }

    /// Replace the maximum relative value with the one specified.
    pub fn max_relative(self, max_relative: A::Epsilon) -> Relative<A, B>
    {
        Relative {
            max_relative,
            ..self
        }
    }

    /// Perform the equality comparison
    pub fn eq(self, lhs: &A, rhs: &B) -> bool
    {
        A::relative_eq(lhs, rhs, self.epsilon, self.max_relative)
    }

    /// Perform the inequality comparison
    pub fn ne(self, lhs: &A, rhs: &B) -> bool
    {
        A::relative_ne(lhs, rhs, self.epsilon, self.max_relative)
    }
}


/// Equality comparisons between two numbers using both the absolute difference and
/// relative based comparisons.
pub trait RelativeEq<Rhs = Self>: AbsDiffEq<Rhs>
    where
        Rhs: ?Sized,
{
    /// The default relative tolerance for testing values that are far-apart.
    ///
    /// This is used when no `max_relative` value is supplied to the [`relative_eq`] macro.
    fn default_max_relative() -> Self::Epsilon;

    /// A test for equality that uses a relative comparison if the values are far apart.
    fn relative_eq(
        &self,
        other: &Rhs,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool;

    /// The inverse of [`RelativeEq::relative_eq`].
    fn relative_ne(
        &self,
        other: &Rhs,
        epsilon: Self::Epsilon,
        max_relative: Self::Epsilon,
    ) -> bool {
        !Self::relative_eq(self, other, epsilon, max_relative)
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Base implementations
///////////////////////////////////////////////////////////////////////////////////////////////////

// Implementation based on: [Comparing Floating Point Numbers, 2012 Edition]
// (https://randomascii.wordpress.com/2012/02/25/comparing-floating-point-numbers-2012-edition/)
macro_rules! impl_relative_eq_float {
    ($T:ident, $U:ident) => {
        impl RelativeEq for $T {

            fn default_max_relative() -> $T {
                $T::EPSILON
            }

            fn relative_eq(&self, other: &$T, epsilon: $T, max_relative: $T) -> bool {
                // Handle same infinities
                if self == other {
                    return true;
                }

                // Handle remaining infinities
                if $T::is_infinite(*self) || $T::is_infinite(*other) {
                    return false;
                }

                let abs_diff = $T::abs(self - other);

                // For when the numbers are really close together
                if abs_diff <= epsilon {
                    return true;
                }

                let abs_self = $T::abs(*self);
                let abs_other = $T::abs(*other);

                let largest = if abs_other > abs_self {
                    abs_other
                } else {
                    abs_self
                };

                // Use a relative difference comparison
                abs_diff <= largest * max_relative
            }
        }
    };
}

impl_relative_eq_float!(f32, i32);
impl_relative_eq_float!(f64, i64);

/// Approximate equality using both the absolute difference and relative based comparisons.
#[macro_export]
macro_rules! relative_eq {
    ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)*) => {
        $crate::algebra::abstr::Relative::default()$(.$opt($val))*.eq(&$lhs, &$rhs)
    };
    ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)*,) => {
        $crate::algebra::abstr::Relative::default()$(.$opt($val))*.eq(&$lhs, &$rhs)
    };
}

/// Approximate inequality using both the absolute difference and relative based comparisons.
#[macro_export]
macro_rules! relative_ne {
    ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)*) => {
        $crate::algebra::abstr::Relative::default()$(.$opt($val))*.ne(&$lhs, &$rhs)
    };
    ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)*,) => {
        $crate::algebra::abstr::Relative::default()$(.$opt($val))*.ne(&$lhs, &$rhs)
    };
}

/// An assertion that delegates to [`relative_eq!`], and panics with a helpful error on failure.
#[macro_export(local_inner_macros)]
macro_rules! assert_relative_eq {
    ($given:expr, $expected:expr $(, $opt:ident = $val:expr)*) => {
        __assert_approx!(relative_eq, $given, $expected $(, $opt = $val)*)
    };
    ($given:expr, $expected:expr $(, $opt:ident = $val:expr)*,) => {
        __assert_approx!(relative_eq, $given, $expected $(, $opt = $val)*)
    };
}

/// An assertion that delegates to [`relative_ne!`], and panics with a helpful error on failure.
#[macro_export(local_inner_macros)]
macro_rules! assert_relative_ne {
    ($given:expr, $expected:expr $(, $opt:ident = $val:expr)*) => {
        __assert_approx!(relative_ne, $given, $expected $(, $opt = $val)*)
    };
    ($given:expr, $expected:expr $(, $opt:ident = $val:expr)*,) => {
        __assert_approx!(relative_ne, $given, $expected $(, $opt = $val)*)
    };
}