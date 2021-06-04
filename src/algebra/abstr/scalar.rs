use crate::algebra::abstr::{
    cast::{FromPrimitive, NumCast, ToPrimitive},
    Sign,
};
use std::fmt::{Debug, Display};
use crate::algebra::abstr::Complex;

#[cfg(feature = "lapack")]
use crate::algebra::abstr::Zero;

/// comparisons, basic numeric operations, and string conversion.
#[cfg(feature = "native")]
pub trait Scalar<Rhs = Self, Output = Self>:
    Sized + Display + ToPrimitive + FromPrimitive + NumCast + Debug + Copy + PartialOrd + Sign
{
}

/// comparisons, basic numeric operations, and string conversion.
#[cfg(feature = "lapack")]
pub trait Scalar<Rhs = Self, Output = Self>:
    Sized
    + Display
    + ToPrimitive
    + FromPrimitive
    + NumCast
    + Debug
    + Copy
    + PartialOrd
    + Sign
    + Blas
    + Lapack
{
}

macro_rules! impl_scalar {
    ($t:ty) => {
        impl Scalar for $t
        {
        }
    };
}

impl_scalar!(f32);
impl_scalar!(f64);
impl_scalar!(Complex<f32>);
impl_scalar!(Complex<f64>);


#[cfg(feature = "lapack")]
pub trait Lapack: Sized + Zero
{
    fn xgehrd(n: i32,
              ilo: i32,
              ihi: i32,
              a: &mut [Self],
              lda: i32,
              tau: &mut [Self],
              work: &mut [Self],
              lwork: i32,
              info: &mut i32);

    fn xgehrd_work_size(n: i32,
                        ilo: i32,
                        ihi: i32,
                        a: &mut [Self],
                        lda: i32,
                        tau: &mut [Self],
                        info: &mut i32)
                        -> i32;

    fn xorghr(n: i32,
              ilo: i32,
              ihi: i32,
              a: &mut [Self],
              lda: i32,
              tau: &[Self],
              work: &mut [Self],
              lwork: i32,
              info: &mut i32);

    fn xorghr_work_size(n: i32,
                        ilo: i32,
                        ihi: i32,
                        a: &mut [Self],
                        lda: i32,
                        tau: &[Self],
                        info: &mut i32)
                        -> i32;

    fn xgeev(jobvl: u8,
             jobvr: u8,
             n: i32,
             a: &mut [Self],
             lda: i32,
             w: &mut [Self],
             vl: &mut [Self],
             ldvl: i32,
             vr: &mut [Self],
             ldvr: i32,
             work: &mut [Self],
             lwork: i32,
             info: &mut i32);

    fn xgeev_work_size(jobvl: u8,
                       jobvr: u8,
                       n: i32,
                       a: &mut [Self],
                       lda: i32,
                       w: &mut [Self],
                       vl: &mut [Self],
                       ldvl: i32,
                       vr: &mut [Self],
                       ldvr: i32,
                       info: &mut i32)
                       -> i32;

    fn xgetrf(m: i32, n: i32, a: &mut [Self], lda: i32, ipiv: &mut [i32], info: &mut i32);

    fn xgeqrf(m: i32,
              n: i32,
              a: &mut [Self],
              lda: i32,
              tau: &mut [Self],
              work: &mut [Self],
              lwork: i32,
              info: &mut i32);

    fn xgeqrf_work_size(m: i32,
                        n: i32,
                        a: &mut [Self],
                        lda: i32,
                        tau: &mut [Self],
                        info: &mut i32)
                        -> i32;

    fn xorgqr(m: i32,
              n: i32,
              k: i32,
              a: &mut [Self],
              lda: i32,
              tau: &mut [Self],
              work: &mut [Self],
              lwork: i32,
              info: &mut i32);

    fn xorgqr_work_size(m: i32,
                        n: i32,
                        k: i32,
                        a: &mut [Self],
                        lda: i32,
                        tau: &mut [Self],
                        info: &mut i32)
                        -> i32;

    fn xgetri(n: i32,
              a: &mut [Self],
              lda: i32,
              ipiv: &[i32],
              work: &mut [Self],
              lwork: i32,
              info: &mut i32);

    fn xgetri_work_size(n: i32, a: &mut [Self], lda: i32, ipiv: &[i32], info: &mut i32) -> i32;

    fn xpotrf(uplo: char, n: i32, a: &mut [Self], lda: i32, info: &mut i32);

    fn xgetrs(n: i32,
              nrhs: i32,
              a: &mut [Self],
              lda: i32,
              ipiv: &mut [i32],
              b: &mut [Self],
              ldb: i32,
              info: &mut i32);
}

#[cfg(feature = "lapack")]
pub trait Blas: Sized + Zero
{
    fn xgemm(transa: u8,
             transb: u8,
             m: i32,
             n: i32,
             k: i32,
             alpha: Self,
             a: &[Self],
             lda: i32,
             b: &[Self],
             ldb: i32,
             beta: Self,
             c: &mut [Self],
             ldc: i32);

    fn xtrsm(side: char,
             uplo: char,
             transa: char,
             diag: char,
             m: i32,
             n: i32,
             alpha: Self,
             a: &[Self],
             lda: i32,
             b: &mut [Self],
             ldb: i32);

    fn xscal(n: i32, a: Self, x: &mut [Self], inc: i32 );

    fn xaxpy(n: i32, a: Self, x: &[Self], incx: i32, y: &mut [Self], incy: i32);
}
