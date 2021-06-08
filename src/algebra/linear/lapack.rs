use crate::algebra::abstr::{Lapack, Zero, Complex};
use lapack;
use lapack_sys as ffi;
use std::os::raw::c_char;


macro_rules! lapack_real (
    ($T: ty, $xgehrd: path, $xorghr: path, $xgeev: path, $xgetrf: path, $xgeqrf: path, $xorgqr: path, $xgetri: path, $xpotrf: path,
    $xgetrs: path)
    => (
        impl Lapack for $T
       	{

       		//Hessenberg
       		fn xgehrd(n: i32, ilo: i32, ihi: i32, a: &mut [Self], lda: i32,
                      tau: &mut [Self], work: &mut [Self], lwork: i32, info: &mut i32)
           	{
                unsafe { $xgehrd(n, ilo, ihi, a, lda, tau, work, lwork, info) }
			}

            fn xgehrd_work_size(n: i32, ilo: i32, ihi: i32, a: &mut [Self], lda: i32,
                                tau: &mut [Self], info: &mut i32) -> i32
            {
                let mut work = [<$T>::zero()];
                let lwork = -1 as i32;

                unsafe { $xgehrd(n, ilo, ihi, a, lda, tau, &mut work, lwork, info) };

                work[0] as i32
            }

         	fn xorghr(n: i32, ilo: i32, ihi: i32, a: &mut [Self], lda: i32, tau: &[Self],
                      work: &mut [Self], lwork: i32, info: &mut i32)
          	{
                unsafe { $xorghr(n, ilo, ihi, a, lda, tau, work, lwork, info) }
            }

            fn xorghr_work_size(n: i32, ilo: i32, ihi: i32, a: &mut [Self], lda: i32, tau: &[Self], info: &mut i32) -> i32 {
                let mut work = [<$T>::zero() ];
                let lwork = -1 as i32;

                unsafe { $xorghr(n, ilo, ihi, a, lda, tau, &mut work, lwork, info) };

                work[0] as i32
            }

            fn xgeev(jobvl: u8, jobvr: u8, n: i32, a: &mut [Self], lda: i32,
                     w: &mut [Self],
                     vl: &mut [Self], ldvl: i32, vr: &mut [Self], ldvr: i32,
                     work: &mut [Self], lwork: i32, info: &mut i32)
          	{
          		let mut wi = Vec::with_capacity(n as usize);
                unsafe { $xgeev(jobvl, jobvr, n, a, lda, w, &mut wi, vl, ldvl, vr, ldvr, work, lwork, info) }
            }


            fn xgeev_work_size(jobvl: u8, jobvr: u8, n: i32, a: &mut [Self], lda: i32,
                               w: &mut [Self], vl: &mut [Self], ldvl: i32,
                               vr: &mut [Self], ldvr: i32, info: &mut i32) -> i32
          	{
                let mut work = [<$T>::zero()];
                let mut wi = Vec::with_capacity(n as usize);
                let lwork = -1 as i32;

                unsafe { $xgeev(jobvl, jobvr, n, a, lda, w, &mut wi, vl, ldvl, vr, ldvr, &mut work, lwork, info) };
                work[0] as i32
			}

			//LU decomposition
			fn xgetrf(m: i32, n: i32, a: &mut [Self], lda: i32, ipiv: &mut [i32], info: &mut i32)
			{
                unsafe { $xgetrf(m, n, a, lda, ipiv, info)}
			}

			fn xgeqrf(m: i32, n: i32, a: &mut [Self], lda: i32, tau: &mut [Self], work: &mut [Self], lwork: i32,
			info: &mut i32)
			{
				unsafe { $xgeqrf(m, n, a, lda, tau, work, lwork, info) };
			}

  			fn xgeqrf_work_size(m: i32, n: i32, a: &mut [Self], lda: i32, tau: &mut [Self], info: &mut i32) -> i32
			{
				let mut work = [<$T>::zero() ];
                let lwork = -1 as i32;

                unsafe { $xgeqrf(m, n, a, lda, tau, &mut work, lwork, info) };
                work[0] as i32
			}

			fn xorgqr(m: i32, n: i32, k: i32, a: &mut [Self], lda: i32, tau: &mut [Self], work: &mut [Self], lwork:
			i32,
			info: &mut i32)
			{
				unsafe { $xorgqr(m, n, k, a, lda, tau, work, lwork, info) };
			}

  			fn xorgqr_work_size(m: i32, n: i32, k: i32, a: &mut [Self], lda: i32, tau: &mut [Self], info: &mut i32) -> i32
			{
				let mut work = [<$T>::zero() ];
                let lwork = -1 as i32;

                unsafe { $xorgqr(m, n, k, a, lda, tau, &mut work, lwork, info) };
                work[0] as i32
			}

			//
			fn xgetri(n: i32, a: &mut [Self], lda: i32, ipiv: &[i32], work: &mut [Self], lwork: i32, info: &mut i32)
			{
                unsafe { $xgetri(n, a, lda, ipiv, work, lwork, info)}
			}

			fn xgetri_work_size(n: i32, a: &mut [Self], lda: i32, ipiv: &[i32], info: &mut i32) -> i32
			{
				let mut work = [ <$T>::zero() ];
                let lwork = -1 as i32;
				unsafe { $xgetri(n, a, lda, ipiv, &mut work, lwork, info) };

				work[0] as i32
			}

			//cholesky
			fn xpotrf(uplo: char, n: i32, a: &mut [Self], lda: i32, info: &mut i32)
			{
				unsafe
				{
					$xpotrf(uplo as u8, n, a, lda, info);
				}
			}

			// solve
			fn xgetrs(n: i32, nrhs: i32, a: &mut [Self], lda: i32, ipiv: &mut [i32], b: &mut [Self], ldb: i32, info: &mut i32)
			{
				unsafe
				{
					$xgetrs('N' as u8, n, nrhs, a, lda, ipiv, b, ldb, info);
				}
			}
      	}
    )
);

lapack_real!(f32,
             lapack::sgehrd,
             lapack::sorghr,
             lapack::sgeev,
             lapack::sgetrf,
             lapack::sgeqrf,
             lapack::sorgqr,
             lapack::sgetri,
             lapack::spotrf,
             lapack::sgetrs);

lapack_real!(f64,
             lapack::dgehrd,
             lapack::dorghr,
             lapack::dgeev,
             lapack::dgetrf,
             lapack::dgeqrf,
             lapack::dorgqr,
             lapack::dgetri,
             lapack::dpotrf,
             lapack::dgetrs);


macro_rules! lapack_complex (
    ($T: ty, $xgehrd: path, $xorghr: path, $xgeev: path, $xgetrf: path, $xgeqrf: path, $xorgqr: path, $xgetri: path, $xpotrf: path, $xgetrs: path)
    => (
		impl Lapack for Complex<$T>
		{
			//Hessenberg
       		fn xgehrd(n: i32, ilo: i32, ihi: i32, a: &mut [Self], lda: i32, tau: &mut [Self], work: &mut [Self], lwork: i32, info: &mut i32)
           	{
                unsafe { $xgehrd(&n, &ilo, &ihi, a.as_mut_ptr() as *mut _, &lda, tau.as_mut_ptr() as *mut _, work.as_mut_ptr() as *mut _, &lwork, info as *mut _) }
			}

            fn xgehrd_work_size(n: i32, ilo: i32, ihi: i32, a: &mut [Self], lda: i32, tau: &mut [Self], info: &mut i32) -> i32
            {
                let mut work = [Self::zero()];
                let lwork = -1 as i32;

                unsafe { $xgehrd(&n, &ilo, &ihi, a.as_mut_ptr() as *mut _, &lda, tau.as_mut_ptr() as *mut _, work.as_mut_ptr() as *mut _, &lwork, info as *mut _) };

                work[0].re as i32
            }


			fn xorghr(n: i32, ilo: i32, ihi: i32, a: &mut [Self], lda: i32, tau: &[Self], work: &mut [Self], lwork: i32, info: &mut i32)
          	{
                unsafe { $xorghr(&n, &ilo, &ihi, a.as_mut_ptr() as *mut _, &lda, tau.as_ptr() as *const _, work.as_mut_ptr() as *mut _, &lwork, info as *mut _) }
            }

            fn xorghr_work_size(n: i32, ilo: i32, ihi: i32, a: &mut [Self], lda: i32, tau: &[Self], info: &mut i32) -> i32 {
                let mut work = [Self::zero() ];
                let lwork = -1 as i32;

                unsafe { $xorghr(&n, &ilo, &ihi, a.as_mut_ptr() as *mut _, &lda, tau.as_ptr() as *const _, work.as_mut_ptr() as *mut _, &lwork, info as *mut _) };

                work[0].re as i32
            }

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
					 info: &mut i32)
			{
				let mut rwork: Vec<$T> = Vec::with_capacity(2*n as usize);
				unsafe
				{
					$xgeev(&(jobvl as c_char), &(jobvr as c_char), &n, a.as_mut_ptr() as *mut _, &lda, w.as_mut_ptr() as *mut _, vl.as_mut_ptr() as *mut _, &ldvl, vr.as_mut_ptr() as *mut _, &ldvr, work.as_mut_ptr() as *mut _, &lwork, rwork.as_mut_ptr() as *mut _, info as *mut _)
				};
			}

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
							   -> i32
			{
				let mut work = [Self::zero()];
				let lwork = -1 as i32;
				let mut rwork: Vec<$T> = Vec::with_capacity(2*n as usize);

				unsafe
				{
					$xgeev(&(jobvl as c_char), &(jobvr as c_char), &n, a.as_mut_ptr() as *mut _, &lda, w.as_mut_ptr() as *mut _, vl.as_mut_ptr() as *mut _, &ldvl, vr.as_mut_ptr() as *mut _, &ldvr, work.as_mut_ptr() as *mut _, &lwork, rwork.as_mut_ptr() as *mut _, info as *mut _)
				};
				work[0].re as i32
			}

			fn xgetrf(m: i32, n: i32, a: &mut [Self], lda: i32, ipiv: &mut [i32], info: &mut i32)
			{
				unsafe
				{
					$xgetrf(&m, &n, a.as_mut_ptr() as *mut _, &lda, ipiv.as_mut_ptr(), info);
				}
			}

			fn xgeqrf(m: i32, n: i32, a: &mut [Self], lda: i32, tau: &mut [Self], work: &mut [Self], lwork: i32, info: &mut i32)
			{
				unsafe
				{
					$xgeqrf(&m, &n, a.as_mut_ptr() as *mut _, &lda, tau.as_mut_ptr() as *mut _, work.as_mut_ptr() as *mut _, &lwork, info as *mut _)
				};
			}

  			fn xgeqrf_work_size(m: i32, n: i32, a: &mut [Self], lda: i32, tau: &mut [Self], info: &mut i32) -> i32
			{
				let mut work = [Self::zero() ];
                let lwork = -1 as i32;

                unsafe {
                	$xgeqrf(&m, &n, a.as_mut_ptr() as *mut _, &lda, tau.as_mut_ptr() as *mut _, work.as_mut_ptr() as *mut _, &lwork, info as *mut _)
				};
                work[0].re as i32
			}

			fn xorgqr(m: i32, n: i32, k: i32, a: &mut [Self], lda: i32, tau: &mut [Self], work: &mut [Self], lwork: i32, info: &mut i32)
			{
				unsafe
				{
					$xorgqr(&m, &n, &k, a.as_mut_ptr() as *mut _, &lda, tau.as_mut_ptr() as *mut _, work.as_mut_ptr() as *mut _, &lwork, info as *mut _)
				};
			}

  			fn xorgqr_work_size(m: i32, n: i32, k: i32, a: &mut [Self], lda: i32, tau: &mut [Self], info: &mut i32) -> i32
			{
				let mut work = [Self::zero()];
                let lwork = -1 as i32;

                unsafe
                {
                	$xorgqr(&m, &n, &k, a.as_mut_ptr() as *mut _, &lda, tau.as_mut_ptr() as *mut _, work.as_mut_ptr() as *mut _, &lwork, info as *mut _)
				};
                work[0].re as i32
			}

			fn xgetri(n: i32,
					  a: &mut [Self],
					  lda: i32,
					  ipiv: &[i32],
					  work: &mut [Self],
					  lwork: i32,
					  info: &mut i32)
			{
				unsafe
				{
					$xgetri(&n, a.as_mut_ptr() as *mut _, &lda, ipiv.as_ptr(), work.as_mut_ptr() as *mut _, &lwork, info)
				};
			}

			fn xgetri_work_size(n: i32,
								a: &mut [Self],
								lda: i32,
								ipiv: &[i32],
								info: &mut i32)
								-> i32
			{
				let mut work = [Self::zero()];
				let lwork = -1 as i32;
				unsafe
				{
					$xgetri(&n, a.as_mut_ptr() as *mut _, &lda, ipiv.as_ptr(), work.as_mut_ptr() as *mut _, &lwork, info)
				};

				work[0].re as i32
			}

			fn xpotrf(uplo: char, n: i32, a: &mut [Self], lda: i32, info: &mut i32)
			{
				unsafe
				{
					$xpotrf(&(uplo as i8), &n, a.as_mut_ptr() as *mut _, &lda, info as *mut _);
				}
			}

			fn xgetrs(n: i32,
					  nrhs: i32,
					  a: &mut [Self],
					  lda: i32,
					  ipiv: &mut [i32],
					  b: &mut [Self],
					  ldb: i32,
					  info: &mut i32)
			{
				unsafe
				{
					$xgetrs(&('N' as c_char), &n, &nrhs, a.as_ptr() as *const _, &lda, ipiv.as_ptr() as *const _, b.as_mut_ptr() as *mut _, &ldb, info as *mut _);
				}
			}
		}
	)
);

lapack_complex!(f32,
             ffi::cgehrd_,
             ffi::cunghr_,
             ffi::cgeev_,
             ffi::cgetrf_,
             ffi::cgeqrf_,
             ffi::cungqr_,
             ffi::cgetri_,
             ffi::cpotrf_,
             ffi::cgetrs_);

lapack_complex!(f64,
             ffi::zgehrd_,
             ffi::zunghr_,
             ffi::zgeev_,
             ffi::zgetrf_,
             ffi::zgeqrf_,
             ffi::zungqr_,
             ffi::zgetri_,
             ffi::zpotrf_,
             ffi::zgetrs_);