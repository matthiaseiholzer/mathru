use crate::{
    algebra::abstr::{
        cast,
        cast::{FromPrimitive, NumCast, ToPrimitive},
        AbelianGroup, AbelianGroupAdd, AbelianGroupMul, Addition, CommutativeRing,
        Field, Group, GroupAdd, GroupMul, Identity, Loop, Magma, MagmaAdd,
        MagmaMul, Monoid, MonoidAdd, MonoidMul, Multiplication, One, Quasigroup, Real, Ring,
        Scalar, Semigroup, SemigroupAdd, SemigroupMul, Sign, Zero,
    },
    elementary::{Exponential, Hyperbolic, Power, Trigonometry},
    algebra::abstr::{AbsDiffEq, RelativeEq},
};

use std::{
    cmp::Ordering,
    fmt,
    fmt::Display,
    ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign},
};

#[cfg(feature = "lapack")]
use crate::algebra::abstr::{Blas, Lapack};

/// Complex number in cartesian form
#[derive(Debug, Clone, Copy)]
// #[repr(C)]
pub struct Complex<T>
    where T: Real
{
    /// Real portion of the complex number
    pub re: T,
    /// Imaginary  portion of the complex number
    pub im: T,
}

impl<T> Complex<T> where T: Real
{
    /// Create a new complex
    pub fn new(re: T, im: T) -> Complex<T>
    {
        Complex { re, im }
    }
}

/// Returns -self = -Re(self) - i Im(self)
impl<T> Neg for Complex<T> where T: Real /* Neg<Output = T> */
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
    where T: Real + AbsDiffEq<Epsilon = T>
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
    where T: Real + AbsDiffEq<Epsilon = T>
    ,Complex<T>: Lapack,  Complex<T>: Blas
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


#[cfg(feature = "native")]
impl<T> Scalar for Complex<T>
    where T: Real
{
}

#[cfg(feature = "lapack")]
impl<T> Scalar for Complex<T>
    where T: Real,
    Complex<T>: Lapack + Blas
{
}

/// Compares to complex numbers
impl<T> PartialEq for Complex<T> where T: Real
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
impl<T> Zero for Complex<T> where T: Real /* Zero */
{
    fn zero() -> Self
    {
        return Complex::new(T::zero(), T::zero());
    }
}

/// Adds two complex numbers
///
/// (a + ib) + (c + id) = (a + c) + i(b + d)
impl<T> Add for Complex<T> where T: Real /* Add<T, Output = T> */
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

#[cfg(feature = "native")]
impl<T> Exponential for Complex<T>
    where T: Real + AbsDiffEq<Epsilon = T>
{
    /// Returns the euler number represented as a complex number
    fn e() -> Self
    {
        Complex { re: T::e(),
            im: T::zero() }
    }

    ///Exponential function
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Exponential};
    /// use mathru::algebra::abstr::Complex;
    ///
    /// let z: Complex<f64> = Complex::new(1.0, 2.0);
    /// let a: Complex<f64> = z.exp();
    /// ```
    fn exp(self: Self) -> Self
    {
        let k: T = self.re.exp();
        Complex { re: k * self.im.cos(),
            im: k * self.im.sin() }
    }

    ///Logiarithm function
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Exponential};
    /// use mathru::algebra::abstr::Complex;
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(5.0_f64.powf(0.5_f64).ln(), 2.0_f64.atan());
    ///
    /// assert_eq!(refer, a.ln());
    /// ```
    fn ln(self: Self) -> Self
    {
        Complex { re: self.abs().re.ln(),
            im: self.arg().re }
    }
}

#[cfg(feature = "lapack")]
impl<T> Exponential for Complex<T>
    where T: Real + AbsDiffEq<Epsilon = T>,
    Complex<T>: Lapack + Blas
{
    /// Returns the euler number represented as a complex number
    fn e() -> Self
    {
        Complex { re: T::e(),
                  im: T::zero() }
    }

    ///Exponential function
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Exponential, num::Complex};
    ///
    /// let z: Complex<f64> = Complex::new(1.0, 2.0);
    /// let a: Complex<f64> = z.exp();
    /// ```
    fn exp(self: Self) -> Self
    {
        let k: T = self.re.exp();
        Complex { re: k * self.im.cos(),
                  im: k * self.im.sin() }
    }

    ///Logiarithm function
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Exponential, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(5.0_f64.powf(0.5_f64).ln(), 2.0_f64.atan());
    ///
    /// assert_eq!(refer, a.ln());
    /// ```
    fn ln(self: Self) -> Self
    {
        Complex { re: self.abs().re.ln(),
                  im: self.arg().re }
    }
}

#[cfg(feature = "native")]
impl<T> Trigonometry for Complex<T>
    where T: Real + AbsDiffEq<Epsilon = T>
{
    /// Returns the mathematics constant PI, represented as a complex number
    fn pi() -> Self
    {
        Complex { re: T::pi(),
            im: T::zero() }
    }

    /// Sinus function
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry};
    /// use mathru::algebra::abstr::Complex;
    ///
    /// let a: f64 = 1.0;
    /// let b: f64 = 2.0;
    /// let z: Complex<f64> = Complex::new(a, b);
    /// let re: f64 = (-(-b).exp() * a.sin() - b.exp() * a.sin()) / -2.0;
    /// let im: f64 = ((-b).exp() * a.cos() - b.exp() * a.cos()) / -2.0;
    ///
    /// let uut: Complex<f64> = z.sin();
    /// let refer: Complex<f64> = Complex::new(re, im);
    ///
    /// assert_eq!(refer, uut);
    /// ```
    fn sin(self: Self) -> Self
    {
        let a: Complex<T> = Complex::new(-self.im, self.re);
        let b: Complex<T> = Complex::new(self.im, -self.re);
        let c: Complex<T> = Complex::new(T::zero(), T::one() + T::one());

        (a.exp() - b.exp()) / c
    }

    /// Cosinus function
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry};
    /// use mathru::algebra::abstr::Complex;
    ///
    /// let a: f64 = 1.0;
    /// let b: f64 = 2.0;
    /// let z: Complex<f64> = Complex::new(a, b);
    /// let re: f64 = ((-b).exp() * a.cos() + b.exp() * (-a).cos()) / 2.0;
    /// let im: f64 = ((-b).exp() * a.sin() + b.exp() * (-a).sin()) / 2.0;
    /// let refer: Complex<f64> = Complex::new(re, im);
    ///
    /// let uut: Complex<f64> = z.cos();
    ///
    /// assert_eq!(refer, uut);
    /// ```
    fn cos(self: Self) -> Self
    {
        let a: Complex<T> = Complex::new(-self.im, self.re);
        let b: Complex<T> = Complex::new(self.im, -self.re);
        let c: Complex<T> = Complex::new(T::one() + T::one(), T::zero());

        (a.exp() + b.exp()) / c
    }

    /// Tangens function
    ///
    /// # Arguments
    ///
    /// self \in \mathbb{C} \setminus \{ k\pi + \frac{\pi}{2} | k \in \mathbb{Z}
    /// \}
    ///
    /// # Panics
    ///
    /// if the argument bounds are not fulfilled
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry};
    /// use mathru::algebra::abstr::Complex;
    ///
    /// let a: f64 = 1.0;
    /// let b: f64 = 2.0;
    /// let z: Complex<f64> = Complex::new(a, b);
    /// let refer: Complex<f64> = z.sin() / z.cos();
    ///
    /// let uut: Complex<f64> = z.tan();
    ///
    /// assert_eq!(refer, uut);
    /// ```
    fn tan(self: Self) -> Self
    {
        return self.sin() / self.cos();
    }

    /// Cotangens function
    ///
    /// # Arguments
    ///
    /// self: \mathbb{C} \setminus \{ \frac{k * \pi}{2} | k \in \mathbb{Z} \}
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry};
    /// use mathru::algebra::abstr::Complex;
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(1.0_f64, 0.0_f64) / a.tan();
    ///
    /// assert_eq!(refer, a.cot());
    /// ```
    fn cot(self: Self) -> Self
    {
        Complex::one() / self.tan()
    }

    /// Secant function
    ///
    /// # Arguments
    ///
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry};
    /// use mathru::algebra::abstr::Complex;
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(1.0_f64, 0.0_f64) / a.cos();
    ///
    /// assert_eq!(refer, a.sec());
    /// ```
    fn sec(self: Self) -> Self
    {
        Complex::one() / self.cos()
    }

    /// Cosecant function
    ///
    /// # Arguments
    ///
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry};
    /// use mathru::algebra::abstr::Complex;
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(1.0_f64, 0.0_f64) / a.sin();
    ///
    /// assert_eq!(refer, a.csc());
    /// ```
    fn csc(self: Self) -> Self
    {
        Complex::one() / self.sin()
    }

    /// Inverse sinus function
    ///
    /// # Arguments
    ///
    /// # Panics
    ///
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry};
    /// use mathru::algebra::abstr::Complex;
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(0.4270785863924768, 1.5285709194809995);
    ///
    /// assert_eq!(refer, a.arcsin());
    /// ```
    fn arcsin(self: Self) -> Self
    {
        let mi: Complex<T> = Complex::new(T::zero(), -T::one());
        let iz: Complex<T> = Complex::new(-self.im, self.re);
        let exp: Complex<T> = Complex::new(T::one() / (T::one() + T::one()), T::zero());

        mi * (iz + (Complex::one() - self * self).pow(exp)).ln()
    }

    /// Inverse cosinus function
    ///
    /// # Arguments
    ///
    /// # Panics
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry};
    /// use mathru::algebra::abstr::Complex;
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(std::f64::consts::PI / 2.0_f64, 0.0_f64) - a.arcsin();
    ///
    /// assert_eq!(refer, a.arccos());
    /// ```
    fn arccos(self: Self) -> Self
    {
        Complex::new(T::pi() / (T::one() + T::one()), T::zero()) - self.arcsin()
    }

    /// Inverse tangens function
    ///
    /// # Arguments
    ///
    /// self: Complex numbers without {-i, i}
    ///
    /// # Panics
    ///
    /// iff self = i or self = -i
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry};
    /// use mathru::algebra::abstr::Complex;
    ///
    /// let a: Complex<f64> = Complex::new(0.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(std::f64::consts::PI / 2.0,
    ///                                        (4.0_f64 / 5.0_f64).atanh() / 2.0_f64);
    ///
    /// assert_eq!(refer, a.arctan());
    /// ```
    fn arctan(self: Self) -> Self
    {
        //		let iz: Complex<T> = Complex::new(-self.im, self.re);
        //		let f: Complex<T> = Complex::new(T::zero(), T::one() / (T::one() + T::one()));
        //		((Complex::one() - iz).ln() - (Complex::one() + iz).ln()) * f
        //		let k: Complex<T> = Complex::new(T::zero(), T::one() + T::one());
        //		((Complex::one() + iz).ln() - (Complex::one() - iz).ln()) / k

        let two: T = T::one() + T::one();
        let re: T;

        if self.re == T::zero()
        {
            if self.im.abs() <= T::one()
            {
                re = T::zero()
            }
            else
            {
                if self.im > T::zero()
                {
                    re = T::pi() / two;
                }
                else
                {
                    re = -T::pi() / two;
                }
            }
        }
        else
        {
            if self.re > T::zero()
            {
                re =
                    (((self.re * self.re + self.im * self.im - T::one()) / (two * self.re)).arctan()
                        + T::pi() / two)
                        / two
            }
            else
            {
                re =
                    (((self.re * self.re + self.im * self.im - T::one()) / (two * self.re)).arctan()
                        - T::pi() / two)
                        / two
            }
        }

        let im: T =
            ((two * self.im) / (self.re * self.re + self.im * self.im + T::one())).artanh() / two;

        Complex::new(re, im)
    }

    fn arctan2(self: Self, _other: Self) -> Self
    {
        unimplemented!()
    }

    /// Inverse cotangens function
    ///
    /// # Arguments
    ///
    /// self: Complex numbers without {-i, i}
    ///
    /// # Panics
    ///
    /// iff self = i or self = -i
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{algebra::abstr::{Complex, One}, elementary::Trigonometry };
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = (Complex::one() / a).arctan();
    ///
    /// assert_eq!(refer, a.arccot());
    /// ```
    fn arccot(self: Self) -> Self
    {
        if self.re == T::zero()
        {
            if self.im == T::one() || self.im == -T::one()
            {
                panic!()
            }
        }
        (Complex::one() / self).arctan()
    }

    /// Inverse secant function
    ///
    /// # Arguments
    ///
    /// self: Complex numbers without {-1, 0, 1}
    ///
    /// # Panics
    ///
    /// iff self = -1 or self = 0 or self = 1
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{algebra::abstr::{Complex, One}, elementary::Trigonometry};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = (Complex::one() / a).arccos();
    ///
    /// assert_eq!(refer, a.arcsec());
    /// ```
    fn arcsec(self: Self) -> Self
    {
        if self.im == T::zero()
        {
            if self.re == -T::one() || self.re == T::zero() || self.re == T::one()
            {
                panic!()
            }
        }

        (Complex::one() / self).arccos()
    }

    /// Inverse cosecant function
    ///
    /// # Arguments
    ///
    ///
    /// # Panics
    ///
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{algebra::abstr::{Complex, One}, elementary::Trigonometry};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = (Complex::one() / a).arcsin();
    ///
    /// assert_eq!(refer, a.arccsc());
    /// ```
    fn arccsc(self: Self) -> Self
    {
        (Complex::one() / self).arcsin()
    }
}

#[cfg(feature = "lapack")]
impl<T> Trigonometry for Complex<T>
    where T: Real + AbsDiffEq<Epsilon = T>,
    Complex<T>: Lapack + Blas
{
    /// Returns the mathematics constant PI, represented as a complex number
    fn pi() -> Self
    {
        Complex { re: T::pi(),
                  im: T::zero() }
    }

    /// Sinus function
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry};
    /// use mathru::algebra::abstr::Complex;
    ///
    /// let a: f64 = 1.0;
    /// let b: f64 = 2.0;
    /// let z: Complex<f64> = Complex::new(a, b);
    /// let re: f64 = (-(-b).exp() * a.sin() - b.exp() * a.sin()) / -2.0;
    /// let im: f64 = ((-b).exp() * a.cos() - b.exp() * a.cos()) / -2.0;
    ///
    /// let uut: Complex<f64> = z.sin();
    /// let refer: Complex<f64> = Complex::new(re, im);
    ///
    /// assert_eq!(refer, uut);
    /// ```
    fn sin(self: Self) -> Self
    {
        let a: Complex<T> = Complex::new(-self.im, self.re);
        let b: Complex<T> = Complex::new(self.im, -self.re);
        let c: Complex<T> = Complex::new(T::zero(), T::one() + T::one());

        (a.exp() - b.exp()) / c
    }

    /// Cosinus function
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry, num::Complex};
    ///
    /// let a: f64 = 1.0;
    /// let b: f64 = 2.0;
    /// let z: Complex<f64> = Complex::new(a, b);
    /// let re: f64 = ((-b).exp() * a.cos() + b.exp() * (-a).cos()) / 2.0;
    /// let im: f64 = ((-b).exp() * a.sin() + b.exp() * (-a).sin()) / 2.0;
    /// let refer: Complex<f64> = Complex::new(re, im);
    ///
    /// let uut: Complex<f64> = z.cos();
    ///
    /// assert_eq!(refer, uut);
    /// ```
    fn cos(self: Self) -> Self
    {
        let a: Complex<T> = Complex::new(-self.im, self.re);
        let b: Complex<T> = Complex::new(self.im, -self.re);
        let c: Complex<T> = Complex::new(T::one() + T::one(), T::zero());

        (a.exp() + b.exp()) / c
    }

    /// Tangens function
    ///
    /// # Arguments
    ///
    /// self \in \mathbb{C} \setminus \{ k\pi + \frac{\pi}{2} | k \in \mathbb{Z}
    /// \}
    ///
    /// # Panics
    ///
    /// if the argument bounds are not fulfilled
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry, num::Complex};
    ///
    /// let a: f64 = 1.0;
    /// let b: f64 = 2.0;
    /// let z: Complex<f64> = Complex::new(a, b);
    /// let refer: Complex<f64> = z.sin() / z.cos();
    ///
    /// let uut: Complex<f64> = z.tan();
    ///
    /// assert_eq!(refer, uut);
    /// ```
    fn tan(self: Self) -> Self
    {
        return self.sin() / self.cos();
    }

    /// Cotangens function
    ///
    /// # Arguments
    ///
    /// self: \mathbb{C} \setminus \{ \frac{k * \pi}{2} | k \in \mathbb{Z} \}
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(1.0_f64, 0.0_f64) / a.tan();
    ///
    /// assert_eq!(refer, a.cot());
    /// ```
    fn cot(self: Self) -> Self
    {
        Complex::one() / self.tan()
    }

    /// Secant function
    ///
    /// # Arguments
    ///
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(1.0_f64, 0.0_f64) / a.cos();
    ///
    /// assert_eq!(refer, a.sec());
    /// ```
    fn sec(self: Self) -> Self
    {
        Complex::one() / self.cos()
    }

    /// Cosecant function
    ///
    /// # Arguments
    ///
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(1.0_f64, 0.0_f64) / a.sin();
    ///
    /// assert_eq!(refer, a.csc());
    /// ```
    fn csc(self: Self) -> Self
    {
        Complex::one() / self.sin()
    }

    /// Inverse sinus function
    ///
    /// # Arguments
    ///
    /// # Panics
    ///
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(0.4270785863924768, 1.5285709194809995);
    ///
    /// assert_eq!(refer, a.arcsin());
    /// ```
    fn arcsin(self: Self) -> Self
    {
        let mi: Complex<T> = Complex::new(T::zero(), -T::one());
        let iz: Complex<T> = Complex::new(-self.im, self.re);
        let exp: Complex<T> = Complex::new(T::one() / (T::one() + T::one()), T::zero());

        mi * (iz + (Complex::one() - self * self).pow(exp)).ln()
    }

    /// Inverse cosinus function
    ///
    /// # Arguments
    ///
    /// # Panics
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(std::f64::consts::PI / 2.0_f64, 0.0_f64) - a.arcsin();
    ///
    /// assert_eq!(refer, a.arccos());
    /// ```
    fn arccos(self: Self) -> Self
    {
        Complex::new(T::pi() / (T::one() + T::one()), T::zero()) - self.arcsin()
    }

    /// Inverse tangens function
    ///
    /// # Arguments
    ///
    /// self: Complex numbers without {-i, i}
    ///
    /// # Panics
    ///
    /// iff self = i or self = -i
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Trigonometry, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(0.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = Complex::new(std::f64::consts::PI / 2.0,
    ///                                        (4.0_f64 / 5.0_f64).atanh() / 2.0_f64);
    ///
    /// assert_eq!(refer, a.arctan());
    /// ```
    fn arctan(self: Self) -> Self
    {
        //		let iz: Complex<T> = Complex::new(-self.im, self.re);
        //		let f: Complex<T> = Complex::new(T::zero(), T::one() / (T::one() + T::one()));
        //		((Complex::one() - iz).ln() - (Complex::one() + iz).ln()) * f
        //		let k: Complex<T> = Complex::new(T::zero(), T::one() + T::one());
        //		((Complex::one() + iz).ln() - (Complex::one() - iz).ln()) / k

        let two: T = T::one() + T::one();
        let re: T;

        if self.re == T::zero()
        {
            if self.im.abs() <= T::one()
            {
                re = T::zero()
            }
            else
            {
                if self.im > T::zero()
                {
                    re = T::pi() / two;
                }
                else
                {
                    re = -T::pi() / two;
                }
            }
        }
        else
        {
            if self.re > T::zero()
            {
                re =
                    (((self.re * self.re + self.im * self.im - T::one()) / (two * self.re)).arctan()
                     + T::pi() / two)
                    / two
            }
            else
            {
                re =
                    (((self.re * self.re + self.im * self.im - T::one()) / (two * self.re)).arctan()
                     - T::pi() / two)
                    / two
            }
        }

        let im: T =
            ((two * self.im) / (self.re * self.re + self.im * self.im + T::one())).artanh() / two;

        Complex::new(re, im)
    }

    fn arctan2(self: Self, _other: Self) -> Self
    {
        unimplemented!()
    }

    /// Inverse cotangens function
    ///
    /// # Arguments
    ///
    /// self: Complex numbers without {-i, i}
    ///
    /// # Panics
    ///
    /// iff self = i or self = -i
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{algebra::abstr::One, elementary::Trigonometry, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = (Complex::one() / a).arctan();
    ///
    /// assert_eq!(refer, a.arccot());
    /// ```
    fn arccot(self: Self) -> Self
    {
        if self.re == T::zero()
        {
            if self.im == T::one() || self.im == -T::one()
            {
                panic!()
            }
        }
        (Complex::one() / self).arctan()
    }

    /// Inverse secant function
    ///
    /// # Arguments
    ///
    /// self: Complex numbers without {-1, 0, 1}
    ///
    /// # Panics
    ///
    /// iff self = -1 or self = 0 or self = 1
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{algebra::abstr::One, elementary::Trigonometry, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = (Complex::one() / a).arccos();
    ///
    /// assert_eq!(refer, a.arcsec());
    /// ```
    fn arcsec(self: Self) -> Self
    {
        if self.im == T::zero()
        {
            if self.re == -T::one() || self.re == T::zero() || self.re == T::one()
            {
                panic!()
            }
        }

        (Complex::one() / self).arccos()
    }

    /// Inverse cosecant function
    ///
    /// # Arguments
    ///
    ///
    /// # Panics
    ///
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{algebra::abstr::One, elementary::Trigonometry, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let refer: Complex<f64> = (Complex::one() / a).arcsin();
    ///
    /// assert_eq!(refer, a.arccsc());
    /// ```
    fn arccsc(self: Self) -> Self
    {
        (Complex::one() / self).arcsin()
    }
}

#[cfg(feature = "native")]
impl<T> Power for Complex<T>
    where T: Real + AbsDiffEq<Epsilon = T>
{
    /// Power
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Power};
    /// use mathru::algebra::abstr::Complex;
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let b: Complex<f64> = Complex::new(-2.0_f64, -1.0_f64);
    /// let refer: Complex<f64> = Complex::new(-0.6006033457684014, -0.07399065302898929);
    ///
    /// assert_eq!(refer, a.pow(b));
    /// ```
    fn pow(self: Self, exp: Self) -> Self
    {
        let r: T = self.abs().re;
        let phi: T = self.arg().re;
        let k: T = r.pow(exp.re) * (-exp.im * phi).exp();
        let theta: T = r.ln() * exp.im + exp.re * phi;
        let re: T = k * theta.cos();
        let im: T = k * theta.sin();
        Complex::new(re, im)
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


#[cfg(feature = "lapack")]
impl<T> Power for Complex<T>
    where T: Real + AbsDiffEq<Epsilon = T>,
          Complex<T>: Lapack + Blas

{
    /// Power
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::{elementary::Power, num::Complex};
    ///
    /// let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    /// let b: Complex<f64> = Complex::new(-2.0_f64, -1.0_f64);
    /// let refer: Complex<f64> = Complex::new(-0.6006033457684014, -0.07399065302898929);
    ///
    /// assert_eq!(refer, a.pow(b));
    /// ```
    fn pow(self: Self, exp: Self) -> Self
    {
        let r: T = self.abs().re;
        let phi: T = self.arg().re;
        let k: T = r.pow(exp.re) * (-exp.im * phi).exp();
        let theta: T = r.ln() * exp.im + exp.re * phi;
        let re: T = k * theta.cos();
        let im: T = k * theta.sin();
        Complex::new(re, im)
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

#[cfg(feature = "native")]
impl<T> Hyperbolic for Complex<T>
    where T: Real + AbsDiffEq<Epsilon = T>
{
    /// Hyperbolic sine
    fn sinh(self: Self) -> Self
    {
        Complex::new(T::zero(), -T::one()) * Complex::new(-self.im, self.re).sin()
    }

    /// Hyperbolic cosine
    fn cosh(self: Self) -> Self
    {
        Complex::new(-self.im, self.re).cos()
    }

    /// Hyperbolic tangens
    fn tanh(self: Self) -> Self
    {
        self.sinh() / self.cosh()
    }

    /// Hyperbolic cotangens
    fn coth(self: Self) -> Self
    {
        self.cosh() / self.sinh()
    }

    /// Hyperbolic secant
    fn sech(self: Self) -> Self
    {
        Complex::new(-self.im, self.re).sec()
    }

    /// Hyperbolic cosecant
    fn csch(self: Self) -> Self
    {
        Complex::new(T::zero(), -T::one()) * Complex::new(-self.im, self.re).csc()
    }

    /// Hyperbolic inverse sine
    ///
    /// # Arguments
    ///
    /// # Panics
    fn arsinh(self: Self) -> Self
    {
        let p: Complex<T> = Complex::new(T::one() / (T::one() + T::one()), T::zero());
        (self + (self * self + Complex::one()).pow(p)).ln()
    }

    /// Hyperbolic inverse cosine
    ///
    /// # Argument
    ///
    /// # Panics
    fn arcosh(self: Self) -> Self
    {
        let p: Complex<T> = Complex::new(T::one() / (T::one() + T::one()), T::zero());
        (self + (self * self - Complex::one()).pow(p)).ln()
    }

    /// Inverse hyperbolic tangent
    ///
    /// # Arguments
    ///
    /// # Panics
    fn artanh(self: Self) -> Self
    {
        let f: Complex<T> = Complex::new(T::one() / (T::one() + T::one()), T::zero());
        ((Complex::one() + self) / (Complex::one() - self)).ln() * f
    }

    /// Inverse hyperbolic cosecant
    ///
    /// # Arguments
    ///
    ///
    /// # Panics
    fn arcoth(self: Self) -> Self
    {
        let f: Complex<T> = Complex::new(T::one() / (T::one() + T::one()), T::zero());
        ((self + Complex::one()) / (self - Complex::one())).ln() * f
    }

    /// Hyperbolic inverse secant
    fn arsech(self: Self) -> Self
    {
        (Complex::one() / self).arcosh()
    }

    // Hyperbolic inverse cosecant
    fn arcsch(self: Self) -> Self
    {
        (Complex::one() / self).arsinh()
    }
}


#[cfg(feature = "lapack")]
impl<T> Hyperbolic for Complex<T>
    where T: Real + AbsDiffEq<Epsilon = T>,
    Complex<T>: Lapack + Blas
{
    /// Hyperbolic sine
    fn sinh(self: Self) -> Self
    {
        Complex::new(T::zero(), -T::one()) * Complex::new(-self.im, self.re).sin()
    }

    /// Hyperbolic cosine
    fn cosh(self: Self) -> Self
    {
        Complex::new(-self.im, self.re).cos()
    }

    /// Hyperbolic tangens
    fn tanh(self: Self) -> Self
    {
        self.sinh() / self.cosh()
    }

    /// Hyperbolic cotangens
    fn coth(self: Self) -> Self
    {
        self.cosh() / self.sinh()
    }

    /// Hyperbolic secant
    fn sech(self: Self) -> Self
    {
        Complex::new(-self.im, self.re).sec()
    }

    /// Hyperbolic cosecant
    fn csch(self: Self) -> Self
    {
        Complex::new(T::zero(), -T::one()) * Complex::new(-self.im, self.re).csc()
    }

    /// Hyperbolic inverse sine
    ///
    /// # Arguments
    ///
    /// # Panics
    fn arsinh(self: Self) -> Self
    {
        let p: Complex<T> = Complex::new(T::one() / (T::one() + T::one()), T::zero());
        (self + (self * self + Complex::one()).pow(p)).ln()
    }

    /// Hyperbolic inverse cosine
    ///
    /// # Argument
    ///
    /// # Panics
    fn arcosh(self: Self) -> Self
    {
        let p: Complex<T> = Complex::new(T::one() / (T::one() + T::one()), T::zero());
        (self + (self * self - Complex::one()).pow(p)).ln()
    }

    /// Inverse hyperbolic tangent
    ///
    /// # Arguments
    ///
    /// # Panics
    fn artanh(self: Self) -> Self
    {
        let f: Complex<T> = Complex::new(T::one() / (T::one() + T::one()), T::zero());
        ((Complex::one() + self) / (Complex::one() - self)).ln() * f
    }

    /// Inverse hyperbolic cosecant
    ///
    /// # Arguments
    ///
    ///
    /// # Panics
    fn arcoth(self: Self) -> Self
    {
        let f: Complex<T> = Complex::new(T::one() / (T::one() + T::one()), T::zero());
        ((self + Complex::one()) / (self - Complex::one())).ln() * f
    }

    /// Hyperbolic inverse secant
    fn arsech(self: Self) -> Self
    {
        (Complex::one() / self).arcosh()
    }

    // Hyperbolic inverse cosecant
    fn arcsch(self: Self) -> Self
    {
        (Complex::one() / self).arsinh()
    }
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

impl<T> Magma<Addition> for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
    fn operate(self, rhs: Self) -> Self
    {
        return Complex::new(self.re + rhs.re, self.im + rhs.im);
    }
}

impl<T> MagmaAdd for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> Magma<Multiplication> for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
    fn operate(self, rhs: Self) -> Self
    {
        return Complex::new(self.re * rhs.re - self.im * rhs.im,
                            self.re * rhs.im + self.im * rhs.re);
    }
}

impl<T> MagmaMul for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> Quasigroup<Addition> for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> Quasigroup<Multiplication> for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{

}

impl<T> AbelianGroup<Addition> for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> AbelianGroupAdd for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> AbelianGroup<Multiplication> for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}
impl<T> AbelianGroupMul for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> Loop<Addition> for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> Loop<Multiplication> for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> CommutativeRing for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> Ring for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> Monoid<Addition> for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> MonoidAdd for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> Monoid<Multiplication> for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}
impl<T> MonoidMul for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> Semigroup<Addition> for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> SemigroupAdd for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> Semigroup<Multiplication> for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> SemigroupMul for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> Field for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> Group<Addition> for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> GroupAdd for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> Group<Multiplication> for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> GroupMul for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
{
}

impl<T> AbsDiffEq for Complex<T> where T: Real + AbsDiffEq<Epsilon = T>
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

impl<T> RelativeEq for Complex<T> where T: Real + RelativeEq<Epsilon = T>
{
    fn default_max_relative() -> Self {
        Complex::new(T::default_max_relative(), T::default_max_relative())
    }

    fn relative_eq(&self, other: &Self, epsilon: Self, max_relative: Self) -> bool {
        self.re.relative_eq(&other.re, epsilon.re, max_relative.re) && self.im.relative_eq(&other.im, epsilon.im, max_relative.im)
    }
}


