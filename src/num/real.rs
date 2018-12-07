use std::ops::{Add, AddAssign, Mul, MulAssign, Sub, SubAssign, Div, DivAssign, Neg};
use std::fmt;
use std::fmt::{Display};
use std::cmp::{Ordering};
use algebra::abstr::{Number, Field};
use algebra::abstr::Real as RealT;
use algebra::abstr::{Semiring, Zero, One};
use elementary::{Exponential, Trigonometry, Power, Hyperbolic};

#[macro_export]
macro_rules! Real
{
    ($v:expr) =>
    {
		(Real::new($v));
    };
}

#[macro_export]
macro_rules! Real32
{
    ($v:expr) =>
    {
		(Real::new($v));
    };
}

#[macro_export]
macro_rules! Real64
{
    ($v:expr) =>
    {
		(Real::new($v));
    };
}


pub type Real32 = Real<f32>;
pub type Real64 = Real<f64>;


#[derive(Debug, Copy, Clone)]
pub struct Real<T>
{
	num: T
}

impl<T> Real<T>
	where T: Neg<Output = T> + PartialOrd + Copy + Zero
{
	pub fn new(num: T) -> Self
	{
		Real
		{
			num: num,
		}
	}

    pub fn abs(self: Self) -> Real<T>
    {
        if self < Real::zero()
        {
           	return self.neg();
        }
        self
    }

}

impl<T> PartialEq for Real<T>
    where T: PartialEq
{
    fn eq<'a, 'b>(self: &'a Self, other: &Self) -> bool
    {
        if self.num == other.num
        {
            return true;
        }
        false
    }
}

impl<T> PartialOrd for Real<T>
	where T: PartialOrd
{
    fn partial_cmp(self: &Self, other: &Self) -> Option<Ordering>
    {
        self.num.partial_cmp(&other.num)
    }
}

impl<T> Display for Real<T>
	where T: Display
{
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result
    {
        // Write strictly the first element into the supplied output
        // stream: `f`. Returns `fmt::Result` which indicates whether the
        // operation succeeded or failed. Note that `write!` uses syntax which
        // is very similar to `println!`.
        write!(f, "{}", self.num)
    }
}

impl<T> Zero for Real<T>
	where T: Zero + Copy
{
    // This cannot be an associated constant, because of bignums.
    fn zero() -> Self
	{
		Real
		{
			num: T::zero()
		}
	}
}

//impl<T> Neg for Real<T>
//	where T: Neg
//{
//	type Output = Real<T>;
//
//    fn neg(self) -> Real<T>
//    {
//    	Real
//		{
//			num: -self.num
//		}
//    }
//}

impl<T> Add for Real<T>
    where T: Add<T, Output = T> + Copy
{
    type Output = Real<T>;

    fn add(self: Self, other: Self) -> Self::Output
    {
        &self + &other
    }
}


impl< 'a, 'b, T> Add<&'b Real<T>> for &'a Real<T>
    where T: Add<T, Output = T> + Copy
{
    type Output = Real<T>;

    fn add(self: Self, rhs: &'b Real<T>) -> Self::Output
    {
        Real{
            num: self.num + rhs.num
        }
    }
}


impl<T> AddAssign for Real<T>
    where T: AddAssign
{
    fn add_assign<'a>(self: &'a mut Real<T>, other: Real<T>)
    {
        self.num += other.num;
    }
}

impl<T> Mul for Real<T>
    where T: Mul<T, Output = T> + Copy
{
    type Output = Real<T>;

    fn mul(self: Self, other: Self) -> Self::Output
    {
        (&self).mul(&other)
    }
}

impl<'a, 'b, T> Mul<&'b Real<T>> for &'a Real<T>
    where T: Mul<T, Output = T> + Clone + Copy
{
    type Output = Real<T>;

    fn mul(self: Self, rhs: &'b Real<T>) -> Self::Output
    {
        Real{
            num: self.num.mul(rhs.num)
        }
    }
}

impl<T>  MulAssign for Real<T>
    where T: MulAssign
{
    fn mul_assign<'a>(self: &'a mut Self, rhs: Real<T>)
    {
        self.num *= rhs.num;
    }
}

impl<T> Sub for Real<T>
    where T: Sub<T, Output = T> + Copy
{
    type Output = Real<T>;

    fn sub(self: Self, rhs: Self) -> Self::Output
    {
        (&self).sub(&rhs)
    }
}

impl<'a, 'b, T> Sub<&'b Real<T>> for &'a Real<T>
    where T: Sub<T, Output = T> + Copy
{
    type Output = Real<T>;

    fn sub(self: Self, rhs: &'b Real<T>) -> Self::Output
    {
        Real{
            num: self.num.sub(rhs.num)
        }
    }
}


impl<T> SubAssign for Real<T>
    where T: SubAssign
{
    fn sub_assign<'a>(self: &'a mut Self, other: Real<T>)
    {
        self.num -= other.num;
    }
}

impl<T> Div for Real<T>
    where T: Div<T, Output = T> + Copy
{
    type Output = Real<T>;

    fn div(self: Self, rhs: Self) ->  Self::Output
    {
        (&self).div(&rhs)
    }
}

impl<'a, 'b, T> Div<&'b Real<T>> for &'a Real<T>
    where T: Div<T, Output = T> + Copy
{
    type Output = Real<T>;

    fn div(self: Self, rhs: &'b Real<T>) -> Self::Output
    {
        Real{
            num: self.num.div(rhs.num)
        }
    }
}


impl<T> DivAssign for Real<T>
    where T: DivAssign
{
    fn div_assign<'a>(self: &'a mut Self, other: Self)
    {
        self.num.div_assign(other.num);
    }
}


impl<T> Neg for Real<T>
    where T: Neg<Output = T>
{
    type Output = Real<T>;

    fn neg(self: Self) -> Self::Output
    {
        Real
        {
            num: -self.num
        }
    }
}





impl<T> One for Real<T>
    where T: One + Copy
{
    fn one() -> Self
    {
        Real
        {
            num: T::one()
        }
    }
}


//impl<T> Abs for Real<T>
//    where T: Abs
//{
//    fn abs<'a>(self: &'a Self) -> Self
//    {
//        Real
//        {
//            num: self.num.abs()
//        }
//    }
//}

impl<T> Trigonometry for Real<T>
	where T: Trigonometry
{
	/// Returns the mathematic constant PI
	fn pi() -> Self
	{
		Real
		{
			num: T::pi()
		}
	}

	/// Sinus
	fn sin(self: Self) -> Self
	{
		Real
		{
			num: self.num.sin()
		}
	}

	/// Cosinus
	fn cos(self: Self) -> Self
	{
		Real
		{
			num: self.num.cos()
		}
	}

	///Tangens
	fn tan(self: Self) -> Self
	{
		Real
		{
			num: self.num.tan()
		}
	}

	fn cot(self: Self) -> Self
	{
		Real
		{
			num: self.num.cot()
		}
	}

	fn sec(self: Self) -> Self
	{
		Real
		{
			num: self.num.sec()
		}
	}

	fn csc(self: Self) -> Self
	{
		Real
		{
			num: self.num.csc()
		}
	}

	fn arcsin(self: Self) -> Self
	{
		Real
		{
			num: self.num.arcsin()
		}
	}

	fn arccos(self: Self) -> Self
	{
		Real
		{
			num: self.num.arccos()
		}
	}
	fn arctan(self: Self) -> Self
	{
		Real
		{
			num: self.num.arctan()
		}
	}

	fn arccot(self: Self) -> Self
	{
		unimplemented!();
	}

	fn arcsec(self: Self) -> Self
	{
		Real
		{
			num: self.num.arccos()
		}
	}

	fn arccsc(self: Self) -> Self
	{
		Real
		{
			num: self.num.arccsc()
		}
	}
}

impl<T> Exponential for Real<T>
	where T: Exponential
{
	fn e() -> Self
	{
		Real
		{
			num: T::e()
		}
	}

	///Exponential function
	fn exp(self: Self) -> Self
	{
		Real
		{
			num: self.num.exp()
		}
	}

	///Logiarithm function
	fn ln(self: Self) -> Self
	{
		Real
		{
			num: self.num.ln()
		}
	}
}

impl<T> Power for Real<T>
	where T: Power
{
	fn pow(self: Self, exp: Self) -> Self
	{
		Real
		{
			num: self.num.pow(exp.num)
		}
	}

	fn root(self: Self, root: Self) -> Self
	{
		Real
		{
			num: self.num.root(root.num)
		}
	}
}

impl<T> Hyperbolic for Real<T>
	where T: Hyperbolic + Field
{
	/// Hyperbolic sine
	fn sinh(self: Self) -> Self
	{
		Real
		{
			num: self.num.sinh()
		}
	}

	/// Hyperbolic cosine
	fn cosh(self: Self) -> Self
	{
		Real
		{
			num: self.num.cosh()
		}
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
    /// extern crate mathru;
    /// use mathru::num::{Real};
   	/// use mathru::elementary::Hyperbolic;
    ///
    /// let x: Real<f64> = Real::new(0.0_f64);
	///
	/// let f: Real<f64> = x.tanh();
	/// let g: Real<f64> = Real::new(0.0_f64);
	/// let abs_difference: Real<f64> = (f - g).abs();
	///
	/// assert!(abs_difference < Real::new(1.0e-10));
    /// ```
	fn tanh(self: Self) -> Self
	{
		Real
		{
			num: self.num.tanh()
		}
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
    /// extern crate mathru;
    /// use mathru::num::{Real};
   	/// use mathru::elementary::Hyperbolic;
    ///
    /// let x: Real<f64> = Real::new(1.0_f64);
	///
	/// let f: Real<f64> = x.coth();
	/// let g: Real<f64> = x.cosh() / x.sinh();
	/// let abs_difference: Real<f64> = (f - g).abs();
	///
	/// assert!(abs_difference < Real::new(1.0e-10));
    /// ```
	fn coth(self: Self) -> Self
	{
		if self == Real::zero()
		{
			panic!();
		}

		Real
		{
			num: self.num.coth()
		}
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
    /// extern crate mathru;
    /// use mathru::num::{Real};
    /// use mathru::elementary::Hyperbolic;
    ///
    /// let x: Real<f64> = Real::new(0.0_f64);
	///
	/// let f: Real<f64> = x.sech();
	/// let g: Real<f64> = Real::new(1.0);
	/// let abs_difference: Real<f64> = (f - g).abs();
	///
	/// assert!(abs_difference < Real::new(1.0e-10));
    /// ```
	fn sech(self: Self) -> Self
	{
		Real
		{
			num: self.num.sech()
		}
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
    /// ```
    /// extern crate mathru;
    /// use mathru::num::{Real};
    /// use mathru::elementary::Hyperbolic;
    ///
    /// let x: Real<f64> = Real::new(1.0_f64);
	///
	/// let f: Real<f64> = x.csch();
	/// let g: Real<f64> = Real::new(1.0) / x.sinh();
	/// let abs_difference: Real<f64> = (f - g).abs();
	///
	/// assert!(abs_difference < Real::new(1.0e-10));
    /// ```
	fn csch(self: Self) -> Self
	{
		Real
		{
			num: self.num.csch()
		}
	}

	/// Hyperbolic inverse sine
	fn arsinh(self: Self) -> Self
	{
		Real
		{
			num: self.num.arsinh()
		}
	}

	/// Hyperbolic inverse cosine
	fn arcosh(self: Self) -> Self
	{
		Real
		{
			num: self.num.arcosh()
		}
	}

	/// Hyperbolic inverse tangens
	fn artanh(self: Self) -> Self
	{
		Real
		{
			num: self.num.artanh()
		}
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
    /// extern crate mathru;
    /// use mathru::num::Real;
   	/// use mathru::algebra::abstr::{Field};
    /// use mathru::elementary::{Exponential, Hyperbolic};
    ///
    /// let x: Real<f64> = Real::new(2.0_f64);
	/// let f: Real<f64> = x.arcoth();
	/// let g: Real<f64> = ((x + Real::new(1.0)) / ( x - Real::new(1.0))).ln() / Real::new(2.0);
	/// let abs_difference: Real<f64> = (f - g).abs();
	///
	/// assert!(abs_difference < Real::new(1.0e-10));
    /// ```
	fn arcoth(self: Self) -> Self
	{
		if -Real::one() <= self && self <= Real::one()
		{
			panic!();
		}

		Real
		{
			num: self.num.arcoth()
		}
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
    /// extern crate mathru;
    /// use mathru::num::Real;
   	/// use mathru::algebra::abstr::{Field};
    /// use mathru::elementary::{Exponential, Hyperbolic};
    ///
    /// let x: Real<f64> = Real::new(0.5_f64);
	/// let f: Real<f64> = x.arsech();
	/// let g: Real<f64> = (Real::new(1.0) / x).arcosh();
	/// let abs_difference: Real<f64> = (f - g).abs();
	///
	/// assert!(abs_difference < Real::new(1.0e-10));
    /// ```
	fn arsech(self: Self) -> Self
	{
		if Real::zero() >= self || self > Real::one()
		{
			panic!();
		}

		Real
		{
			num: self.num.arsech()
		}
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
    /// extern crate mathru;
    /// use mathru::num::Real;
   	/// use mathru::algebra::abstr::{Field};
    /// use mathru::elementary::{Exponential, Hyperbolic};
    ///
    /// let x: Real<f64> = Real::new(2.0_f64);
	/// let f: Real<f64> = x.arcsch();
	/// let g: Real<f64> = (Real::new(1.0) / x).arsinh();
	/// let abs_difference: Real<f64> = (f - g).abs();
	///
	/// assert!(abs_difference < Real::new(1.0e-10));
    /// ```
	fn arcsch(self: Self) -> Self
	{
		if self == Real::zero()
		{
			panic!()
		}

		Real
		{
			num: self.num.arcsch()
		}
	}
}