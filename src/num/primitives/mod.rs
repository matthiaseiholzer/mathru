use std::{f32, f64, i128, i16, i32, i64, i8, u128, u16, u32, u64, u8};

use crate::{
    algebra::abstr::cast::{AsPrimitive, ToPrimitive},
    elementary::{Exponential, Hyperbolic, Power, Trigonometry},
};

#[cfg(feature = "native")]
pub mod native;

macro_rules! impl_to_primitive_int_to_int {
    ($SrcT:ty, $DstT:ty, $slf:expr) => {{ $slf as $DstT }};
}

macro_rules! impl_to_primitive_int_to_uint {
    ($SrcT:ty, $DstT:ty, $slf:expr) => {{
        return $slf as $DstT;
    }};
}

macro_rules! impl_to_primitive_int {
    ($T:ty) => {
        impl ToPrimitive for $T
        {
            fn to_i8(&self) -> i8
            {
                impl_to_primitive_int_to_int!($T, i8, *self)
            }

            fn to_i16(&self) -> i16
            {
                impl_to_primitive_int_to_int!($T, i16, *self)
            }

            fn to_i32(&self) -> i32
            {
                impl_to_primitive_int_to_int!($T, i32, *self)
            }

            fn to_i64(&self) -> i64
            {
                impl_to_primitive_int_to_int!($T, i64, *self)
            }

            fn to_i128(&self) -> i128
            {
                impl_to_primitive_int_to_int!($T, i128, *self)
            }

            fn to_u8(&self) -> u8
            {
                impl_to_primitive_int_to_uint!($T, u8, *self)
            }

            fn to_u16(&self) -> u16
            {
                impl_to_primitive_int_to_uint!($T, u16, *self)
            }

            fn to_u32(&self) -> u32
            {
                impl_to_primitive_int_to_uint!($T, u32, *self)
            }

            fn to_u64(&self) -> u64
            {
                impl_to_primitive_int_to_uint!($T, u64, *self)
            }

            fn to_u128(&self) -> u128
            {
                impl_to_primitive_int_to_uint!($T, u128, *self)
            }

            fn to_f32(&self) -> f32
            {
                return *self as f32;
            }

            fn to_f64(&self) -> f64
            {
                return *self as f64;
            }
        }
    };
}

impl_to_primitive_int!(i8);
impl_to_primitive_int!(i16);
impl_to_primitive_int!(i32);
impl_to_primitive_int!(i64);
impl_to_primitive_int!(i128);

macro_rules! impl_to_primitive_uint_to_int {
    ($DstT:ty, $slf:expr) => {{
        return $slf as $DstT;
    }};
}

macro_rules! impl_to_primitive_uint_to_uint {
    ($SrcT:ty, $DstT:ty, $slf:expr) => {{
        return $slf as $DstT;
    }};
}

macro_rules! impl_to_primitive_uint {
    ($T:ty) => {
        impl ToPrimitive for $T
        {
            fn to_i8(&self) -> i8
            {
                impl_to_primitive_uint_to_int!(i8, *self)
            }

            fn to_i16(&self) -> i16
            {
                impl_to_primitive_uint_to_int!(i16, *self)
            }

            fn to_i32(&self) -> i32
            {
                impl_to_primitive_uint_to_int!(i32, *self)
            }

            fn to_i64(&self) -> i64
            {
                impl_to_primitive_uint_to_int!(i64, *self)
            }

            fn to_i128(&self) -> i128
            {
                impl_to_primitive_uint_to_int!(i128, *self)
            }

            fn to_u8(&self) -> u8
            {
                impl_to_primitive_uint_to_uint!($T, u8, *self)
            }

            fn to_u16(&self) -> u16
            {
                impl_to_primitive_uint_to_uint!($T, u16, *self)
            }

            fn to_u32(&self) -> u32
            {
                impl_to_primitive_uint_to_uint!($T, u32, *self)
            }

            fn to_u64(&self) -> u64
            {
                impl_to_primitive_uint_to_uint!($T, u64, *self)
            }

            fn to_u128(&self) -> u128
            {
                impl_to_primitive_uint_to_uint!($T, u128, *self)
            }

            fn to_f32(&self) -> f32
            {
                return *self as f32;
            }

            fn to_f64(&self) -> f64
            {
                return *self as f64;
            }
        }
    };
}

impl_to_primitive_uint!(u8);
impl_to_primitive_uint!(u16);
impl_to_primitive_uint!(u32);
impl_to_primitive_uint!(u64);
impl_to_primitive_uint!(u128);
//impl_to_primitive_uint!(f32);
//impl_to_primitive_uint!(f64);

macro_rules! impl_to_primitive_float_to_float {
    ($SrcT:ident, $DstT:ident, $slf:expr) => {
        return $slf as $DstT;
    };
}

macro_rules! impl_to_primitive_float {
    ($T:ident) => {
        impl ToPrimitive for $T
        {
            fn to_i8(&self) -> i8
            {
                return *self as i8;
            }

            fn to_i16(&self) -> i16
            {
                return *self as i16;
            }

            fn to_i32(&self) -> i32
            {
                return *self as i32;
            }

            fn to_i64(&self) -> i64
            {
                return *self as i64;
            }

            fn to_i128(&self) -> i128
            {
                return *self as i128;
            }

            fn to_u8(&self) -> u8
            {
                return *self as u8;
            }

            fn to_u16(&self) -> u16
            {
                return *self as u16;
            }

            fn to_u32(&self) -> u32
            {
                return *self as u32;
            }

            fn to_u64(&self) -> u64
            {
                return *self as u64;
            }

            fn to_u128(&self) -> u128
            {
                return *self as u128;
            }

            fn to_f32(&self) -> f32
            {
                impl_to_primitive_float_to_float!($T, f32, *self)
            }

            fn to_f64(&self) -> f64
            {
                impl_to_primitive_float_to_float!($T, f64, *self)
            }
        }
    };
}

impl_to_primitive_float!(f32);
impl_to_primitive_float!(f64);

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

macro_rules! exponential_impl {
    ($t:ty, $e: expr) => {
        impl Exponential for $t
        {
            ///
            fn e() -> Self
            {
                $e
            }

            ///Exponential function
            fn exp(self: Self) -> Self
            {
                self.exp()
            }

            ///Logarithm function
            fn ln(self: Self) -> Self
            {
                self.ln()
            }
        }
    };
}

exponential_impl!(f32, f32::consts::E);
exponential_impl!(f64, f64::consts::E);

macro_rules! power_impl {
    ($t:ty) => {
        impl Power for $t
        {
            fn pow(self: Self, exp: Self) -> Self
            {
                return self.powf(exp);
            }

            fn root(self: Self, root: Self) -> Self
            {
                return self.powf(1.0 / root);
            }

            fn sqrt(self: Self) -> Self
            {
                return self.powf(0.5);
            }
        }
    };
}

power_impl!(f32);
power_impl!(f64);

macro_rules! power_impl_integer {
    ($t:ty) => {
        impl Power for $t
        {
            fn pow(self: Self, _exp: Self) -> Self
            {
                unimplemented!();
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
    };
}

power_impl_integer!(u8);
power_impl_integer!(u16);
power_impl_integer!(u32);
power_impl_integer!(u64);
power_impl_integer!(u128);
power_impl_integer!(usize);
power_impl_integer!(i8);
power_impl_integer!(i16);
power_impl_integer!(i32);
power_impl_integer!(i64);
power_impl_integer!(i128);
power_impl_integer!(isize);

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

macro_rules! impl_as_primitive
{
    (@ $T: ty => $(#[$cfg:meta])* impl $U: ty ) =>
    {
        $(#[$cfg])*
        impl AsPrimitive<$U> for $T
        {
            #[inline]
            fn as_(self) -> $U
            {
            	self as $U
            }
        }
    };

    (@ $T: ty => { $( $U: ty ),* } ) =>
    {
    	$(
        	impl_as_primitive!(@ $T => impl $U);
    	)*
    };
    ($T: ty => { $( $U: ty ),* } ) =>
    {
        impl_as_primitive!(@ $T => { $( $U ),* });
        impl_as_primitive!(@ $T => { u8, u16, u32, u64, usize });
        impl_as_primitive!(@ $T => #[cfg(has_i128)] impl u128);
        impl_as_primitive!(@ $T => { i8, i16, i32, i64, isize });
        impl_as_primitive!(@ $T => #[cfg(has_i128)] impl i128);
    };
}

impl_as_primitive!(u8 => { char, f32, f64 });
impl_as_primitive!(i8 => { f32, f64 });
impl_as_primitive!(u16 => { f32, f64 });
impl_as_primitive!(i16 => { f32, f64 });
impl_as_primitive!(u32 => { f32, f64 });
impl_as_primitive!(i32 => { f32, f64 });
impl_as_primitive!(u64 => { f32, f64 });
impl_as_primitive!(i64 => { f32, f64 });
#[cfg(has_i128)]
impl_as_primitive!(u128 => { f32, f64 });
#[cfg(has_i128)]
impl_as_primitive!(i128 => { f32, f64 });
impl_as_primitive!(usize => { f32, f64 });
impl_as_primitive!(isize => { f32, f64 });
impl_as_primitive!(f32 => { f32, f64 });
impl_as_primitive!(f64 => { f32, f64 });
impl_as_primitive!(char => { char });
impl_as_primitive!(bool => {});

