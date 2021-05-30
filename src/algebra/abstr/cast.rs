//! Cast module

//Copied from https://github.com/rust-num/num-traits/blob/master/src/cast.rs

/// A generic trait for converting a value to a number.
pub trait ToPrimitive
{
    /// Converts the value of `self` to an `i8`.
    fn to_i8(&self) -> i8
    {
        self.to_i64().to_i8()
    }

    /// Converts the value of `self` to an `i16`.
    fn to_i16(&self) -> i16
    {
        self.to_i64().to_i16()
    }

    /// Converts the value of `self` to an `i32`.
    fn to_i32(&self) -> i32
    {
        self.to_i64().to_i32()
    }

    /// Converts the value of `self` to an `i64`.
    fn to_i64(&self) -> i64;

    /// Converts the value of `self` to an `i128`.
    fn to_i128(&self) -> i128;

    /// Converts the value of `self` to an `u8`.
    fn to_u8(&self) -> u8
    {
        self.to_u64().to_u8()
    }

    /// Converts the value of `self` to an `u16`.
    fn to_u16(&self) -> u16
    {
        self.to_u64().to_u16()
    }

    /// Converts the value of `self` to an `u32`.
    fn to_u32(&self) -> u32
    {
        self.to_u64().to_u32()
    }

    /// Converts the value of `self` to an `u64`.
    fn to_u64(&self) -> u64;

    /// Converts the value of `self` to an `u128`.
    fn to_u128(&self) -> u128;

    /// Converts the value of `self` to an `f32`.
    fn to_f32(&self) -> f32
    {
        self.to_f64().to_f32()
    }

    /// Converts the value of `self` to an `f64`.
    fn to_f64(&self) -> f64;
    //	{
    //		self.to_f64().and_then(|x: f64| x.to_f64())
    //	}
}

/// A generic trait for converting a number to a value.
pub trait FromPrimitive: Sized
{
    /// Convert an `i8` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_i8(n: i8) -> Self
    {
        FromPrimitive::from_i64(n as i64)
    }

    /// Convert an `i16` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_i16(n: i16) -> Self
    {
        FromPrimitive::from_i64(n as i64)
    }

    /// Convert an `i32` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_i32(n: i32) -> Self
    {
        FromPrimitive::from_i64(n as i64)
    }

    /// Convert an `i64` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_i64(n: i64) -> Self;

    /// Convert an `i128` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_i128(n: i128) -> Self;

    /// Convert an `u8` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_u8(n: u8) -> Self
    {
        FromPrimitive::from_u64(n as u64)
    }

    /// Convert an `u16` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.

    fn from_u16(n: u16) -> Self
    {
        FromPrimitive::from_u64(n as u64)
    }

    /// Convert an `u32` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_u32(n: u32) -> Self
    {
        FromPrimitive::from_u64(n as u64)
    }

    /// Convert an `u64` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_u64(n: u64) -> Self;

    /// Convert an `u128` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_u128(n: u128) -> Self;

    /// Convert a `f32` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_f32(n: f32) -> Self
    {
        FromPrimitive::from_f64(n as f64)
    }

    /// Convert a `f64` to return an optional value of this type. If the
    /// type cannot be represented by this value, the `None` is returned.
    fn from_f64(n: f64) -> Self;
    //	{
    //		FromPrimitive::from_f64(n as f64)
    //	}
}

/// Cast from one machine scalar to another.
pub fn cast<T: NumCast, U: NumCast>(n: T) -> U
{
    NumCast::from(n)
}

/// An interface for casting between machine scalars.
pub trait NumCast: Sized + ToPrimitive
{
    /// Creates a number from another value that can be converted into
    /// a primitive via the `ToPrimitive` trait.
    fn from<T: ToPrimitive>(n: T) -> Self;
}

pub trait AsPrimitive<T>: 'static + Copy
    where T: 'static + Copy
{
    /// Convert a value to another, using the `as` operator.
    fn as_(self) -> T;
}

macro_rules! impl_from_primitive {
    ($T:ty, $to_ty:ident) => {
        impl FromPrimitive for $T
        {
            fn from_i8(n: i8) -> $T
            {
                n.$to_ty()
            }

            fn from_i16(n: i16) -> $T
            {
                n.$to_ty()
            }

            fn from_i32(n: i32) -> $T
            {
                n.$to_ty()
            }

            fn from_i64(n: i64) -> $T
            {
                n.$to_ty()
            }

            fn from_i128(n: i128) -> $T
            {
                n.$to_ty()
            }

            fn from_u8(n: u8) -> $T
            {
                n.$to_ty()
            }

            fn from_u16(n: u16) -> $T
            {
                n.$to_ty()
            }

            fn from_u32(n: u32) -> $T
            {
                n.$to_ty()
            }

            fn from_u64(n: u64) -> $T
            {
                n.$to_ty()
            }

            fn from_u128(n: u128) -> $T
            {
                n.$to_ty()
            }

            fn from_f32(n: f32) -> $T
            {
                n.$to_ty()
            }

            fn from_f64(n: f64) -> $T
            {
                n.$to_ty()
            }
        }
    };
}

impl_from_primitive!(i8, to_i8);
impl_from_primitive!(i16, to_i16);
impl_from_primitive!(i32, to_i32);
impl_from_primitive!(i64, to_i64);
impl_from_primitive!(i128, to_i128);
impl_from_primitive!(u8, to_u8);
impl_from_primitive!(u16, to_u16);
impl_from_primitive!(u32, to_u32);
impl_from_primitive!(u64, to_u64);
impl_from_primitive!(u128, to_u128);
impl_from_primitive!(f32, to_f32);
impl_from_primitive!(f64, to_f64);

macro_rules! impl_num_cast {
    ($T:ty, $conv:ident) => {
        impl NumCast for $T
        {
            fn from<N: ToPrimitive>(n: N) -> $T
            {
                n.$conv()
            }
        }
    };
}

impl_num_cast!(u8, to_u8);
impl_num_cast!(u16, to_u16);
impl_num_cast!(u32, to_u32);
impl_num_cast!(u64, to_u64);
impl_num_cast!(u128, to_u128);
impl_num_cast!(i8, to_i8);
impl_num_cast!(i16, to_i16);
impl_num_cast!(i32, to_i32);
impl_num_cast!(i64, to_i64);
impl_num_cast!(i128, to_i128);
impl_num_cast!(f32, to_f32);
impl_num_cast!(f64, to_f64);

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
