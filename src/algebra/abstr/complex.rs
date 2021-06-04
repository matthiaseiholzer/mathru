use crate::{
    algebra::abstr::{
        cast,
        cast::{FromPrimitive, NumCast, ToPrimitive},
        AbelianGroup, AbelianGroupAdd, AbelianGroupMul, Addition, CommutativeRing,
        Field, Group, GroupAdd, GroupMul, Identity, Loop, Magma, MagmaAdd,
        MagmaMul, Monoid, MonoidAdd, MonoidMul, Multiplication, One, Quasigroup, Real, Ring,
        Semigroup, SemigroupAdd, SemigroupMul, Sign, Zero,
    },
    algebra::abstr::{AbsDiffEq, RelativeEq},
};

use std::{
    cmp::Ordering,
    fmt,
    fmt::Display,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

/// Complex number in cartesian form
#[derive(Debug, Clone, Copy)]
pub struct Complex<T>
{
    /// Real portion of the complex number
    pub re: T,
    /// Imaginary  portion of the complex number
    pub im: T,
}

impl<T> Complex<T>
    where T: Real
{
    /// Create a new complex
    pub fn new(re: T, im: T) -> Complex<T>
    {
        Complex { re, im }
    }
}

/// Returns -self = -Re(self) - i Im(self)
impl<T> Neg for Complex<T> where T: Real
{
    type Output = Complex<T>;

    fn neg(self: Self) -> Complex<T>
    {
        Complex { re: -self.re,
                  im: -self.im }
    }
}

/// Complex numbers can not be compared
///
/// # Panics
///
/// # FIXME
impl<T> PartialOrd for Complex<T>
    where T: Real
{
    fn partial_cmp(self: &Self, other: &Self) -> Option<Ordering>
    {
        if self == other {
            return Some(Ordering::Equal)
        }
        else
        {
            if self.abs().to_f64() > other.abs().to_f64()
            {
                return Some(Ordering::Greater)
            }
            else
            {
                return Some(Ordering::Less)
            }
        }
    }
}

#[cfg(feature = "native")]
impl<T> Complex<T>
    where T: Real
{
    /// Returns the complex conjugate
    /// conj(self) = Re(self) - i Im(self)
    pub fn conj(self: Self) -> Complex<T>
    {
        Complex { re: self.re,
                  im: -self.im }
    }

    /// Returns the argument of the complex number
    pub fn arg(self: Self) -> Self
    {
        Complex { re: self.im.arctan2(self.re),
                  im: T::zero() }
    }
}

#[cfg(feature = "lapack")]
impl<T> Complex<T>
    where T: Real
{
    /// Returns the complex conjugate
    /// conj(self) = Re(self) - i Im(self)
    pub fn conj(self: Self) -> Complex<T>
    {
        Complex { re: self.re,
            im: -self.im }
    }

    /// Returns the argument of the complex number
    pub fn arg(self: Self) -> Self
    {
        Complex { re: self.im.arctan2(self.re),
            im: T::zero() }
    }
}

impl<T> Sign for Complex<T>
    where T: Real
{
    fn sign(self: &Self) -> Self
    {
        unimplemented!()
    }

    /// Absolute value of the complex number
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{
    ///     algebra::abstr::{Complex, cast::ToPrimitive, Sign},
    /// };
    ///
    /// let a: Complex<f64> = Complex::new(1.0, 2.0);
    /// let refer: f64 = (5.0_f64).sqrt();
    /// assert_eq!(refer, a.abs().to_f64());
    /// ```
    fn abs(self: &Self) -> Self
    {
        let root: T = T::from_f64(2.0);
        Complex { re: (self.re * self.re + self.im * self.im).root(root),
                  im: T::zero() }
    }

    fn is_positive(self: &Self) -> bool
    {
        unimplemented!();
    }

    fn is_negative(&self) -> bool
    {
        unimplemented!();
    }
}


// #[cfg(feature = "native")]
// impl<T> Scalar for Complex<T>
//     where T: Real
// {
// }

// #[cfg(feature = "lapack")]
// impl<T> Scalar for Complex<T>
//     where T: Real,
//         Complex<T>: Lapack + Blas
// {
// }

/// Compares to complex numbers
impl<T> PartialEq for Complex<T>
    where T: Real
{
    fn eq(self: &Self, rhs: &Self) -> bool
    {
        if self.re == rhs.re && self.im == rhs.im
        {
            return true;
        }
        false
    }
}

impl<T> Display for Complex<T> where T: Real /* Display */
{
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result
    {
        write!(f, "{} + {}i", self.re, self.im)
    }
}

/// Returns 0 + i0
impl<T> Zero for Complex<T>
    where T: Real
{
    fn zero() -> Self
    {
        return Complex::new(T::zero(), T::zero());
    }
}

/// Adds two complex numbers
///
/// (a + ib) + (c + id) = (a + c) + i(b + d)
impl<T> Add for Complex<T>
    where T: Real
{
    type Output = Complex<T>;

    fn add(self, rhs: Complex<T>) -> Complex<T>
    {
        Complex { re: self.re + rhs.re,
                  im: self.im + rhs.im }
    }
}

impl<'a, 'b, T> Add<&'b Complex<T>> for &'a Complex<T>
    where T: Real
{
    type Output = Complex<T>;

    fn add(self: Self, rhs: &'b Complex<T>) -> Self::Output
    {
        Complex { re: self.re + rhs.re,
                  im: self.im + rhs.im }
    }
}

impl<T> AddAssign for Complex<T> where T: Real /* AddAssign */
{
    fn add_assign(self: &mut Self, other: Self)
    {
        self.re += other.re;
        self.im += other.im;
    }
}

/// Returns 1 + i0
impl<T> One for Complex<T> where T: Real
{
    fn one() -> Self
    {
        return Complex::new(T::one(), T::zero());
    }
}

/// Multiplies two complex numbers
///
/// (a + ib)(c + id) = (ac - bd) + i(bc + ad)
impl<T> Mul for Complex<T>
    where T: Real
{
    type Output = Complex<T>;

    fn mul(self: Self, other: Self) -> Self::Output
    {
        Complex { re: self.re * other.re - self.im * other.im,
                  im: self.im * other.re + self.re * other.im }
    }
}

impl<'a, 'b, T> Mul<&'b Complex<T>> for &'a Complex<T>
    where T: Real
{
    type Output = Complex<T>;

    fn mul(self: Self, rhs: &'b Complex<T>) -> Self::Output
    {
        Complex { re: self.re * rhs.re - self.im * rhs.im,
                  im: self.im * rhs.re + self.re * rhs.im }
    }
}

impl<T> MulAssign for Complex<T>
    where T: Real
{
    fn mul_assign(self: &mut Self, other: Self)
    {
        let temp: Complex<T> = *self;
        self.re = temp.re * other.re - temp.im * other.im;
        self.im = temp.im * other.re + temp.re * other.im;
    }
}

///Subtracts two complex numbers
///
/// (a + ib) - (c + id) = (a - c) + i(b - d)
impl<T> Sub for Complex<T> where T: Real
{
    type Output = Complex<T>;

    fn sub(self: Self, rhs: Self) -> Self::Output
    {
        &self - &rhs
    }
}

impl<'a, 'b, T> Sub<&'b Complex<T>> for &'a Complex<T>
    where T: Real
{
    type Output = Complex<T>;

    fn sub(self: Self, rhs: &'b Complex<T>) -> Self::Output
    {
        Complex { re: self.re - rhs.re,
                  im: self.im - rhs.im }
    }
}

impl<T> SubAssign for Complex<T> where T: Real
{
    fn sub_assign(self: &mut Self, other: Self)
    {
        self.re -= other.re;
        self.im -= other.im;
    }
}

/// Divides two complex numbers
impl<T> Div for Complex<T>
    where T: Real
{
    type Output = Complex<T>;

    fn div(self: Self, rhs: Self) -> Self::Output
    {
        &self / &rhs
    }
}

impl<'a, 'b, T> Div<&'b Complex<T>> for &'a Complex<T>
    where T: Real
{
    type Output = Complex<T>;

    fn div(self: Self, rhs: &'b Complex<T>) -> Self::Output
    {
        let divisor: T = rhs.re * rhs.re + rhs.im * rhs.im;

        if divisor == T::zero()
        {
            panic!("rhs * rhs.conj() == 0")
        }
        let quot: Complex<T> = self * &Complex::new(rhs.re, -rhs.im);

        return Complex::new(quot.re / divisor, quot.im / divisor);
    }
}

impl<T> DivAssign for Complex<T>
    where T: Real
{
    fn div_assign(self: &mut Self, other: Self)
    {
        *self = *self / other;
    }
}

macro_rules! impl_to_primitive {
    ($ty:ty, $to:ident) => {
        fn $to(&self) -> $ty
        {
            return self.re.$to();
            //            if self.im == T::zero()
            //            {
            //            	self.re.$to()
            //            }
            //            else
            //            {
            //            	None
            //            }
        }
    };
}

/// Returns None if Complex part is non-zero
impl<T: ToPrimitive> ToPrimitive for Complex<T>
    where T: Real
{
    impl_to_primitive!(u8, to_u8);

    impl_to_primitive!(u16, to_u16);

    impl_to_primitive!(u32, to_u32);

    impl_to_primitive!(u64, to_u64);

    impl_to_primitive!(u128, to_u128);

    impl_to_primitive!(i8, to_i8);

    impl_to_primitive!(i16, to_i16);

    impl_to_primitive!(i32, to_i32);

    impl_to_primitive!(i64, to_i64);

    impl_to_primitive!(i128, to_i128);

    impl_to_primitive!(f32, to_f32);

    impl_to_primitive!(f64, to_f64);
}

/// A generic trait for converting a number to a value.
impl<T> FromPrimitive for Complex<T> where T: Real
{
    /// Convert an `i64` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_i64(_n: i64) -> Self
    {
        unimplemented!();
    }

    /// Convert an `i128` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_i128(_n: i128) -> Self
    {
        unimplemented!();
    }

    /// Convert an `u64` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_u64(n: u64) -> Self
    {
        Complex { re: cast::cast(n),
            im: T::zero() }
    }

    /// Convert an `u128` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_u128(n: u128) -> Self
    {
        Complex { re: cast::cast(n),
            im: T::zero() }
    }

    /// Convert a `f64` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_f64(n: f64) -> Self
    {
        Complex { re: cast::cast(n),
            im: T::zero() }
    }
}

/// An interface for casting between machine scalars.
impl<T> NumCast for Complex<T> where T: Real
{
    /// Creates a number from another value that can be converted into
    /// a primitive via the `ToPrimitive` trait.
    fn from<K: ToPrimitive>(n: K) -> Self
    {
        Complex { re: cast::cast(n.to_f64()),
                  im: T::zero() }
    }
}

impl<T> Identity<Addition> for Complex<T> where T: Identity<Addition> + Real
{
    fn id() -> Self
    {
        return Complex::new(Identity::<Addition>::id(), Identity::<Addition>::id());
    }
}

impl<T> Identity<Multiplication> for Complex<T>
    where T: Identity<Multiplication> + Identity<Addition> + Real
{
    fn id() -> Self
    {
        return Complex::new(Identity::<Multiplication>::id(), Identity::<Addition>::id());
    }
}

impl<T> Magma<Addition> for Complex<T> where T: Real
{
    fn operate(self, rhs: Self) -> Self
    {
        return Complex::new(self.re + rhs.re, self.im + rhs.im);
    }
}

impl<T> MagmaAdd for Complex<T> where T: Real
{
}

impl<T> Magma<Multiplication> for Complex<T> where T: Real
{
    fn operate(self, rhs: Self) -> Self
    {
        return Complex::new(self.re * rhs.re - self.im * rhs.im,
                            self.re * rhs.im + self.im * rhs.re);
    }
}

impl<T> MagmaMul for Complex<T> where T: Real
{
}

impl<T> Quasigroup<Addition> for Complex<T> where T: Real
{
}

impl<T> Quasigroup<Multiplication> for Complex<T> where T: Real
{

}

impl<T> AbelianGroup<Addition> for Complex<T> where T: Real
{
}

impl<T> AbelianGroupAdd for Complex<T> where T: Real
{
}

impl<T> AbelianGroup<Multiplication> for Complex<T> where T: Real
{
}
impl<T> AbelianGroupMul for Complex<T> where T: Real
{
}

impl<T> Loop<Addition> for Complex<T> where T: Real
{
}

impl<T> Loop<Multiplication> for Complex<T> where T: Real
{
}

impl<T> CommutativeRing for Complex<T> where T: Real
{
}

impl<T> Ring for Complex<T> where T: Real
{
}

impl<T> Monoid<Addition> for Complex<T> where T: Real
{
}

impl<T> MonoidAdd for Complex<T> where T: Real
{
}

impl<T> Monoid<Multiplication> for Complex<T> where T: Real
{
}
impl<T> MonoidMul for Complex<T> where T: Real
{
}

impl<T> Semigroup<Addition> for Complex<T> where T: Real
{
}

impl<T> SemigroupAdd for Complex<T> where T: Real
{
}

impl<T> Semigroup<Multiplication> for Complex<T> where T: Real
{
}

impl<T> SemigroupMul for Complex<T> where T: Real
{
}

impl<T> Field for Complex<T> where T: Real
{
}

impl<T> Group<Addition> for Complex<T> where T: Real
{
}

impl<T> GroupAdd for Complex<T> where T: Real
{
}

impl<T> Group<Multiplication> for Complex<T> where T: Real
{
}

impl<T> GroupMul for Complex<T>
    where T: Real
{
}

// impl<T> Real for Complex<T>
//     where T: Real
// {
//
// }

impl<T> AbsDiffEq for Complex<T>
    where T: Real
{
    type Epsilon = Self;

    fn default_epsilon() -> Self
    {
        Complex::new(T::default_epsilon(), T::default_epsilon())
    }

    fn abs_diff_eq(&self, other: &Complex<T>, epsilon: Self) -> bool
    {
        self.re.abs_diff_eq(&other.re, epsilon.re) && self.im.abs_diff_eq(&other.im, epsilon.im)
    }
}

impl<T> RelativeEq for Complex<T> where T: Real
{
    fn default_max_relative() -> Self {
        Complex::new(T::default_max_relative(), T::default_max_relative())
    }

    fn relative_eq(&self, other: &Self, epsilon: Self, max_relative: Self) -> bool {
        self.re.relative_eq(&other.re, epsilon.re, max_relative.re) && self.im.relative_eq(&other.im, epsilon.im, max_relative.im)
    }
}


