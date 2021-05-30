use std::{f32, f64};

/// The requisite parameters for testing for approximate equality using a
/// absolute difference based comparison.
///
/// This is not normally used directly, rather via the
/// `assert_abs_diff_{eq|ne}!` and `abs_diff_{eq|ne}!` macros.
///
/// # Example
///
/// ```rust
/// use std::f64;
/// use mathru::algebra::abstr::AbsDiff;
///
/// AbsDiff::default().eq(&1.0, &1.0);
/// AbsDiff::default().epsilon(f64::EPSILON).eq(&1.0, &1.0);
/// ```
pub struct AbsDiff<A, B = A>
    where
        A: AbsDiffEq<B> + ?Sized,
        B: ?Sized,
{
    /// The tolerance to use when testing values that are close together.
    pub epsilon: A::Epsilon,
}

impl<A, B> Default for AbsDiff<A, B>
    where
        A: AbsDiffEq<B> + ?Sized,
        B: ?Sized,
{
    fn default() -> AbsDiff<A, B> {
        AbsDiff {
            epsilon: A::default_epsilon(),
        }
    }
}

impl<A, B> AbsDiff<A, B>
    where
        A: AbsDiffEq<B> + ?Sized,
        B: ?Sized,
{
    /// Replace the epsilon value with the one specified.
    pub fn epsilon(self, epsilon: A::Epsilon) -> AbsDiff<A, B> {
        AbsDiff { epsilon, ..self }
    }

    /// Perform the equality comparison
    pub fn eq(self, lhs: &A, rhs: &B) -> bool
    {
        A::abs_diff_eq(lhs, rhs, self.epsilon)
    }

    /// Perform the inequality comparison
    pub fn ne(self, lhs: &A, rhs: &B) -> bool
    {
        A::abs_diff_ne(lhs, rhs, self.epsilon)
    }
}

/// Equality that is defined using the absolute difference of two numbers.
pub trait AbsDiffEq<Rhs = Self>: PartialEq<Rhs>
    where
        Rhs: ?Sized,
{
    /// Used for specifying relative comparisons.
    type Epsilon;

    /// The default tolerance to use when testing values that are close together.
    ///
    /// This is used when no `epsilon` value is supplied to the [`abs_diff_eq!`], [`relative_eq!`],
    fn default_epsilon() -> Self::Epsilon;

    /// A test for equality that uses the absolute difference to compute the approximate
    /// equality of two numbers.
    fn abs_diff_eq(&self, other: &Rhs, epsilon: Self::Epsilon) -> bool;

    /// The inverse of [`AbsDiffEq::abs_diff_eq`].
    fn abs_diff_ne(&self, other: &Rhs, epsilon: Self::Epsilon) -> bool
    {
        !Self::abs_diff_eq(self, other, epsilon)
    }
}

macro_rules! impl_abs_diff_eq
{
    ($T:ident, $default_epsilon:expr) =>
    {
        impl AbsDiffEq for $T
        {
            type Epsilon = $T;

            fn default_epsilon() -> $T
            {
                $default_epsilon
            }

            fn abs_diff_eq(&self, other: &$T, epsilon: $T) -> bool
            {
                (
                    if self > other
                    {
                        self - other
                    } else {
                        other - self
                    }
                ) <= epsilon
            }
        }
    };
}

impl_abs_diff_eq!(f32, f32::EPSILON);
impl_abs_diff_eq!(f64, f64::EPSILON);

/// Approximate equality using the absolute difference.
#[macro_export]
macro_rules! abs_diff_eq {
    ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)*) => {
        $crate::algebra::abstr::AbsDiff::default()$(.$opt($val))*.eq(&$lhs, &$rhs)
    };
    ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)*,) => {
        $crate::algebra::abstr::AbsDiff::default()$(.$opt($val))*.eq(&$lhs, &$rhs)
    };
}

/// Approximate inequality using the absolute difference.
#[macro_export]
macro_rules! abs_diff_ne {
    ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)*) => {
        $crate::algebra::abstr::AbsDiff::default()$(.$opt($val))*.ne(&$lhs, &$rhs)
    };
    ($lhs:expr, $rhs:expr $(, $opt:ident = $val:expr)*,) => {
        $crate::algebra::abstr::AbsDiff::default()$(.$opt($val))*.ne(&$lhs, &$rhs)
    };
}

#[doc(hidden)]
#[macro_export]
macro_rules! __assert_approx {
    ($eq:ident, $given:expr, $expected:expr) => {{
        let (given, expected) = (&($given), &($expected));

        if !$eq!(*given, *expected) {
            panic!(
"assert_{}!({}, {})
    left  = {:?}
    right = {:?}
",
                stringify!($eq),
                stringify!($given),
                stringify!($expected),
                given, expected,
            );
        }
    }};
    ($eq:ident, $given:expr, $expected:expr, $($opt:ident = $val:expr),+) => {{
        let (given, expected) = (&($given), &($expected));

        if !$eq!(*given, *expected, $($opt = $val),+) {
            panic!(
"assert_{}!({}, {}, {})
    left  = {:?}
    right = {:?}
",
                stringify!($eq),
                stringify!($given),
                stringify!($expected),
                stringify!($($opt = $val),+),
                given, expected,
            );
        }
    }};
}

/// An assertion that delegates to [`abs_diff_eq!`], and panics with a helpful error on failure.
#[macro_export(local_inner_macros)]
macro_rules! assert_abs_diff_eq {
    ($given:expr, $expected:expr $(, $opt:ident = $val:expr)*) => {
        __assert_approx!(abs_diff_eq, $given, $expected $(, $opt = $val)*)
    };
    ($given:expr, $expected:expr $(, $opt:ident = $val:expr)*,) => {
        __assert_approx!(abs_diff_eq, $given, $expected $(, $opt = $val)*)
    };
}

/// An assertion that delegates to [`abs_diff_ne!`], and panics with a helpful error on failure.
#[macro_export(local_inner_macros)]
macro_rules! assert_abs_diff_ne {
    ($given:expr, $expected:expr $(, $opt:ident = $val:expr)*) => {
        __assert_approx!(abs_diff_ne, $given, $expected $(, $opt = $val)*)
    };
    ($given:expr, $expected:expr $(, $opt:ident = $val:expr)*,) => {
        __assert_approx!(abs_diff_ne, $given, $expected $(, $opt = $val)*)
    };
}
