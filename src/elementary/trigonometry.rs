/// Trigonometric functions
///
///<a href="https://en.wikipedia.org/wiki/Trigonometric_functions">https://en.wikipedia
/// .org/wiki/Trigonometric_functions</a>
pub trait Trigonometry
{
    /// Returns the mathematic constant PI
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
