use crate::algebra::abstr::{Real, Complex};
use crate::elementary::exponential::Exponential;
use crate::algebra::abstr::One;
use crate::elementary::power::Power;


/// Trigonometric functions
///
///<a href="https://en.wikipedia.org/wiki/Trigonometric_functions">https://en.wikipedia
/// .org/wiki/Trigonometric_functions</a>
pub trait Trigonometry
{
    /// Returns the mathematics constant PI
    fn pi() -> Self;

    /// Sinus function
    fn sin(self: Self) -> Self;

    /// Cosinus function
    fn cos(self: Self) -> Self;

    /// Tangens function
    fn tan(self: Self) -> Self;

    /// Cotangens function
    fn cot(self: Self) -> Self;

    /// Secant function
    fn sec(self: Self) -> Self;

    /// Cosecant function
    fn csc(self: Self) -> Self;

    /// Inverse sinus function
    fn arcsin(self: Self) -> Self;

    /// Inverse cosinus function
    fn arccos(self: Self) -> Self;

    /// Inverse tangens function
    fn arctan(self: Self) -> Self;

    fn arctan2(self: Self, other: Self) -> Self;

    /// Inverse cosecant function
    fn arccot(self: Self) -> Self;

    /// Inverse secant function
    fn arcsec(self: Self) -> Self;

    // Inverse cosecant function
    fn arccsc(self: Self) -> Self;
}

macro_rules! trigonometry_impl {
    ($t:ty, $pi: expr) => {
        impl Trigonometry for $t
        {
            /// Returns the mathematic constant PI
            fn pi() -> Self
            {
                $pi
            }

            /// Sinus
            fn sin(self: Self) -> Self
            {
                self.sin()
            }

            /// Cosinus
            fn cos(self: Self) -> Self
            {
                self.cos()
            }

            ///Tangens
            fn tan(self: Self) -> Self
            {
                self.tan()
            }

            //
            fn cot(self: Self) -> Self
            {
                1.0 / self.tan()
            }

            /// Secant
            ///
            /// # Panics
            ///
            /// self = n pi + pi/2 n \in Z
            fn sec(self: Self) -> Self
            {
                1.0 / self.cos()
            }

            fn csc(self: Self) -> Self
            {
                1.0 / self.sin()
            }

            /// Inverse sine function
            ///
            /// # Arguemnts
            ///
            /// -1.0 <= x <= 1.0
            ///
            /// # Panics
            ///
            /// |x| > 1.0
            fn arcsin(self: Self) -> Self
            {
                if self.abs() > 1.0
                {
                    panic!();
                }

                self.asin()
            }

            /// Inverse cosine function
            ///
            /// # Arguemnts
            ///
            /// -1.0 <= x <= 1.0
            ///
            /// # Panics
            ///
            /// |x| > 1.0
            fn arccos(self: Self) -> Self
            {
                if self.abs() > 1.0
                {
                    panic!();
                }

                self.acos()
            }

            /// Computes the arctangent of a number
            fn arctan(self: Self) -> Self
            {
                self.atan()
            }

            /// Computes the arctangent
            fn arctan2(self: Self, other: Self) -> Self
            {
                self.atan2(other)
            }

            fn arccot(self: Self) -> Self
            {
                if self == 0.0
                {
                    return 0.0;
                }

                if self > 0.0
                {
                    return (1.0 / self).atan();
                }
                else
                {
                    return (1.0 / self).atan();
                }
            }

            fn arcsec(self: Self) -> Self
            {
                (1.0 / self).acos()
            }

            fn arccsc(self: Self) -> Self
            {
                (1.0 / self).asin()
            }
        }
    };
}

trigonometry_impl!(f32, std::f32::consts::PI);
trigonometry_impl!(f64, std::f64::consts::PI);

impl<T> Trigonometry for Complex<T>
    where T: Real
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