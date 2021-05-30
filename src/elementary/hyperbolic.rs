/// Hyperbolic functions
///
///<a href="https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions">https://en.wikipedia.org/wiki/Inverse_hyperbolic_functions</a>
pub trait Hyperbolic
{
    /// Hyperbolic sine
    fn sinh(self: Self) -> Self;

    /// Hyperbolic cosine
    fn cosh(self: Self) -> Self;

    /// Hyperbolic tangens
    fn tanh(self: Self) -> Self;

    /// Hyperbolic cotangens
    fn coth(self: Self) -> Self;

    /// Hyperbolic secant
    fn sech(self: Self) -> Self;

    /// Hyperbolic cosecant
    fn csch(self: Self) -> Self;

    /// Inverse hyperbolic  sine
    fn arsinh(self: Self) -> Self;

    /// Inverse hyperbolic cosine
    fn arcosh(self: Self) -> Self;

    /// Inverse hyperbolic tangens
    fn artanh(self: Self) -> Self;

    /// Inverse hyperbolic cosecant
    fn arcoth(self: Self) -> Self;

    /// Inverse hyperbolic secant
    fn arsech(self: Self) -> Self;

    /// Inverse hyperbolic cosecant
    fn arcsch(self: Self) -> Self;
}

macro_rules! hyperbolic_impl {
    ($t:ty) => {
        impl Hyperbolic for $t
        {
            /// Hyperbolic sine
            fn sinh(self: Self) -> Self
            {
                self.sinh()
            }

            /// Hyperbolic cosine
            fn cosh(self: Self) -> Self
            {
                self.cosh()
            }

            /// Hyperbolic tangens
            ///
            /// # Arguments
            ///
            /// * `self` :
            ///
            /// # Example
            ///
            /// ```
            /// use mathru::elementary::Hyperbolic;
            ///
            /// let x: f64 = 0.0_f64;
            ///
            /// let f: f64 = x.tanh();
            /// ```
            fn tanh(self: Self) -> Self
            {
                self.tanh()
            }

            /// Hyperbolic cotangens
            ///
            /// # Arguments
            ///
            /// * `self` : != 0.0
            ///
            /// # Panic
            ///
            /// iff self == 0.0
            ///
            /// # Example
            ///
            /// ```
            /// use mathru::elementary::Hyperbolic;
            ///
            /// let x: f64 = 1.0_f64;
            ///
            /// let f: f64 = x.coth();
            /// ```
            fn coth(self: Self) -> Self
            {
                if self == 0.0
                {
                    panic!();
                }

                self.cosh() / self.sinh()
            }

            /// Hyperbolic secant
            ///
            /// # Arguments
            ///
            /// * `self` :
            ///
            /// # Example
            ///
            /// ```
            /// use mathru::elementary::Hyperbolic;
            ///
            /// let x: f64 = 0.0_f64;
            ///
            /// let f: f64 = x.sech();
            /// ```
            fn sech(self: Self) -> Self
            {
                1.0 / self.cosh()
            }

            /// Hyperbolic cosecant
            ///
            /// # Arguments
            ///
            /// * `self` : != 0.0
            ///
            /// # Panics
            ///
            /// if  self == 0
            ///
            /// # Example
            ///
            ///
            /// ```
            /// use mathru::elementary::Hyperbolic;
            ///
            /// let x: f64 = 1.0_f64;
            ///
            /// let f: f64 = x.csch();
            /// ```
            fn csch(self: Self) -> Self
            {
                if self == 0.0
                {
                    panic!();
                }
                1.0 / self.sinh()
            }

            /// Hyperbolic inverse sine
            fn arsinh(self: Self) -> Self
            {
                self.asinh()
            }

            /// Hyperbolic inverse cosine
            fn arcosh(self: Self) -> Self
            {
                self.acosh()
            }

            /// Hyperbolic inverse tangens
            fn artanh(self: Self) -> Self
            {
                if -1.0 >= self || self >= 1.0
                {
                    panic!();
                }

                self.atanh()
            }

            /// Hyperbolic inverse cotan
            ///
            /// # Arguments
            ///
            /// * `self`  -1.0 > self, self > 1.0
            ///
            /// # Panics
            ///
            /// if  -1.0 <= self && self <= 1.0
            ///
            /// # Example
            ///
            /// ```
            /// use mathru::{
            ///     algebra::abstr::Field,
            ///     elementary::{Exponential, Hyperbolic},
            /// };
            ///
            /// let x: f64 = 2.0_f64;
            /// let f: f64 = x.arcoth();
            /// ```
            fn arcoth(self: Self) -> Self
            {
                if -1.0 <= self && self <= 1.0
                {
                    panic!();
                }

                ((self + 1.0) / (self - 1.0)).ln() / 2.0
            }

            /// Hyperbolic inverse secant
            ///
            /// # Arguments
            ///
            /// * `self`  0.0 < self <= 1.0
            ///
            /// # Panics
            ///
            /// if  0.0 >= self || self > 1.0
            ///
            /// # Example
            ///
            /// ```
            /// use mathru::elementary::{Exponential, Hyperbolic};
            ///
            /// let x: f64 = 0.5_f64;
            /// let f: f64 = x.arsech();
            /// let g: f64 = (1.0 / x).arcosh();
            /// ```
            fn arsech(self: Self) -> Self
            {
                if 0.0 >= self || self > 1.0
                {
                    panic!();
                }

                (1.0 / self).arcosh()
            }

            /// Hyperbolic inverse cosecant
            ///
            /// # Arguments
            ///
            /// * `self`  <> 0.0
            ///
            /// # Panics
            ///
            /// iff self = 0.0
            ///
            /// # Example
            ///
            /// ```
            /// use mathru::{
            ///     algebra::abstr::Field,
            ///     elementary::{Exponential, Hyperbolic},
            /// };
            ///
            /// let x: f64 = 2.0_f64;
            /// let f: f64 = x.arcsch();
            /// ```
            fn arcsch(self: Self) -> Self
            {
                (1.0 / self).arsinh()
            }
        }
    };
}

hyperbolic_impl!(f32);
hyperbolic_impl!(f64);
