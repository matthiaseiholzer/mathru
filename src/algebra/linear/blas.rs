use crate::algebra::abstr::{Blas, Complex};
use blas_sys as ffi;

macro_rules! blas_impl (
    ($T: ty, $xgemm: path, $xtrsm: path, $xscal: path, $xaxpy: path)
    => (
        impl Blas for $T
       	{
       		fn xgemm(transa: u8, transb: u8, m: i32, n: i32, k: i32, alpha: Self,
    		a: &[Self],
    		lda: i32,
    		b: &[Self],
    		ldb: i32,
    		beta: Self,
    		c: &mut [Self],
    		ldc: i32 )
			{
				unsafe
				{
					$xgemm(
						&(transa as i8),
						&(transb as i8),
						&m,
						&n ,
						&k,
						&alpha as *const _ as *const _,
						a.as_ptr() as *const _,
						&lda,
						b.as_ptr() as *const _,
						&ldb,
						&beta as *const _ as *const _,
						c.as_mut_ptr() as *mut _,
						&ldc
						)
				}
			}

			fn xtrsm(side: char, uplo: char, transa: char, diag: char, m: i32, n: i32, alpha: Self, a: &[Self], lda: i32, b: &mut
			[Self], ldb: i32)
			{
				unsafe
				{
					$xtrsm(
						&(side as i8),
						&(uplo as i8),
						&(transa as i8),
						&(diag as i8),
						&m,
						&n,
						&alpha as *const _ as *const _,
						a.as_ptr() as *const _,
						&lda,
						b.as_mut_ptr() as *mut _,
						&ldb
					);
				}
			}

			fn xscal(n: i32, a: Self, x: &mut [Self], inc: i32)
			{
			    unsafe
			    {
			        $xscal(
			        	&n,
			        	&a as *const _ as *const _,
			        	x.as_mut_ptr() as *mut _,
			        	&inc
					);
			    }
			}

			fn xaxpy(n: i32, a: Self, x: &[Self], incx: i32, y: &mut [Self], incy: i32)
			{
			    unsafe
			    {
                    $xaxpy(
                        &n,
                        &a as *const _ as *const _,
                        x.as_ptr() as *const _,
                        &incx,
                        y.as_mut_ptr() as *mut _,
                        &incy,
                    )
                }
			}
		}
	)
);

blas_impl!(f32, ffi::sgemm_, ffi::strsm_, ffi::sscal_, ffi::saxpy_);
blas_impl!(f64, ffi::dgemm_, ffi::dtrsm_, ffi::dscal_, ffi::daxpy_);
blas_impl!(Complex<f32>, ffi::cgemm_, ffi::ctrsm_, ffi::cscal_, ffi::caxpy_);
blas_impl!(Complex<f64>, ffi::zgemm_, ffi::ztrsm_, ffi::zscal_, ffi::zaxpy_);
