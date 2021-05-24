use crate::algebra::abstr::{Lapack, Zero, Complex};
use lapack;
use lapack_sys as ffi;


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

            fn xorghr_work_size(n: i32, ilo: i32, ihi: i32, a: &mut [Self], lda: i32,
                                tau: &[Self], info: &mut i32) -> i32 {
                let mut work = [<$T>::zero() ];
                let lwork = -1 as i32;

                unsafe { $xorghr(n, ilo, ihi, a, lda, tau, &mut work, lwork, info) };

                work[0] as i32
            }

            fn xgeev(jobvl: u8, jobvr: u8, n: i32, a: &mut [Self], lda: i32,
                     wr: &mut [Self], wi: &mut [Self],
                     vl: &mut [Self], ldvl: i32, vr: &mut [Self], ldvr: i32,
                     work: &mut [Self], lwork: i32, info: &mut i32)
          	{
                unsafe { $xgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, work, lwork, info) }
            }


            fn xgeev_work_size(jobvl: u8, jobvr: u8, n: i32, a: &mut [Self], lda: i32,
                               wr: &mut [Self], wi: &mut [Self], vl: &mut [Self], ldvl: i32,
                               vr: &mut [Self], ldvr: i32, info: &mut i32) -> i32
          	{
                let mut work = [ <$T>::zero() ];
                let lwork = -1 as i32;

                unsafe { $xgeev(jobvl, jobvr, n, a, lda, wr, wi, vl, ldvl, vr, ldvr, &mut work, lwork, info) };
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

  			fn xorgqr_work_size(m: i32, n: i32, k: i32, a: &mut [Self], lda: i32, tau: &mut [Self], info: &mut i32) ->
  			 i32
			{
				let mut work = [<$T>::zero() ];
                let lwork = -1 as i32;

                unsafe { $xorgqr(m, n, k, a, lda, tau, &mut work, lwork, info) };
                work[0] as i32
			}

			//
			fn xgetri(n: i32, a: &mut [Self], lda: i32, ipiv: &mut [i32], work: &mut [Self], lwork: i32, info: &mut i32)
			{
                unsafe { $xgetri(n, a, lda, ipiv, work, lwork, info)}
			}

			fn xgetri_work_size(n: i32, a: &mut [Self], lda: i32, ipiv: &mut [i32], info: &mut i32) -> i32
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


impl Lapack for Complex<f64>
{
	fn xgehrd(_n: i32,
			  _ilo: i32,
			  _ihi: i32,
			  _a: &mut [Self],
			  _lda: i32,
			  _tau: &mut [Self],
			  _work: &mut [Self],
			  _lwork: i32,
			  _info: &mut i32)
	{
		unimplemented!();
	}

	fn xgehrd_work_size(_n: i32,
						_ilo: i32,
						_ihi: i32,
						_a: &mut [Self],
						_lda: i32,
						_tau: &mut [Self],
						_info: &mut i32)
						-> i32
	{
		unimplemented!();
	}

	fn xorghr(_n: i32,
			  _ilo: i32,
			  _ihi: i32,
			  _a: &mut [Self],
			  _lda: i32,
			  _tau: &[Self],
			  _work: &mut [Self],
			  _lwork: i32,
			  _info: &mut i32)
	{
		unimplemented!();
	}

	fn xorghr_work_size(_n: i32,
						_ilo: i32,
						_ihi: i32,
						_a: &mut [Self],
						_lda: i32,
						_tau: &[Self],
						_info: &mut i32)
						-> i32
	{
		unimplemented!();
	}

	fn xgeev(_jobvl: u8,
			 _jobvr: u8,
			 _n: i32,
			 _a: &mut [Self],
			 _lda: i32,
			 _wr: &mut [Self],
			 _wi: &mut [Self],
			 _vl: &mut [Self],
			 _ldvl: i32,
			 _vr: &mut [Self],
			 _ldvr: i32,
			 _work: &mut [Self],
			 _lwork: i32,
			 _info: &mut i32)
	{
		unimplemented!();
	}

	fn xgeev_work_size(_jobvl: u8,
					   _jobvr: u8,
					   _n: i32,
					   _a: &mut [Self],
					   _lda: i32,
					   _wr: &mut [Self],
					   _wi: &mut [Self],
					   _vl: &mut [Self],
					   _ldvl: i32,
					   _vr: &mut [Self],
					   _ldvr: i32,
					   _info: &mut i32)
					   -> i32
	{
		unimplemented!();
	}

	fn xgetrf(m: i32, n: i32, a: &mut [Self], lda: i32, ipiv: &mut [i32], info: &mut i32)
	{
		unsafe {
			ffi::zgetrf_(&m, &n, a.as_mut_ptr() as *mut _, &lda, ipiv.as_mut_ptr(), info);
		}
	}

	fn xgeqrf(_m: i32,
			  _n: i32,
			  _a: &mut [Self],
			  _lda: i32,
			  _tau: &mut [Self],
			  _work: &mut [Self],
			  _lwork: i32,
			  _info: &mut i32)
	{
		unimplemented!();
	}

	fn xgeqrf_work_size(_m: i32,
						_n: i32,
						_a: &mut [Self],
						_lda: i32,
						_tau: &mut [Self],
						_info: &mut i32)
						-> i32
	{
		unimplemented!();
	}

	fn xorgqr(_m: i32,
			  _n: i32,
			  _k: i32,
			  _a: &mut [Self],
			  _lda: i32,
			  _tau: &mut [Self],
			  _work: &mut [Self],
			  _lwork: i32,
			  _info: &mut i32)
	{
		unimplemented!();
	}

	fn xorgqr_work_size(_m: i32,
						_n: i32,
						_k: i32,
						_a: &mut [Self],
						_lda: i32,
						_tau: &mut [Self],
						_info: &mut i32)
						-> i32
	{
		unimplemented!();
	}

	fn xgetri(_n: i32,
			  _a: &mut [Self],
			  _lda: i32,
			  _ipiv: &mut [i32],
			  _work: &mut [Self],
			  _lwork: i32,
			  _info: &mut i32)
	{
		unimplemented!();
	}

	fn xgetri_work_size(_n: i32,
						_a: &mut [Self],
						_lda: i32,
						_ipiv: &mut [i32],
						_info: &mut i32)
						-> i32
	{
		unimplemented!();
	}

	fn xpotrf(_uplo: char, _n: i32, _a: &mut [Self], _lda: i32, _info: &mut i32)
	{
		unimplemented!();
	}

	fn xgetrs(_n: i32,
			  _nrhs: i32,
			  _a: &mut [Self],
			  _lda: i32,
			  _ipiv: &mut [i32],
			  _b: &mut [Self],
			  _ldb: i32,
			  _info: &mut i32)
	{
		unimplemented!();
	}
}

impl Lapack for Complex<f32>
{
	fn xgehrd(_n: i32,
			  _ilo: i32,
			  _ihi: i32,
			  _a: &mut [Self],
			  _lda: i32,
			  _tau: &mut [Self],
			  _work: &mut [Self],
			  _lwork: i32,
			  _info: &mut i32)
	{
		unimplemented!()
	}

	fn xgehrd_work_size(_n: i32,
						_ilo: i32,
						_ihi: i32,
						_a: &mut [Self],
						_lda: i32,
						_tau: &mut [Self],
						_info: &mut i32)
						-> i32
	{
		unimplemented!();
	}

	fn xorghr(_n: i32,
			  _ilo: i32,
			  _ihi: i32,
			  _a: &mut [Self],
			  _lda: i32,
			  _tau: &[Self],
			  _work: &mut [Self],
			  _lwork: i32,
			  _info: &mut i32)
	{
		unimplemented!();
	}

	fn xorghr_work_size(_n: i32,
						_ilo: i32,
						_ihi: i32,
						_a: &mut [Self],
						_lda: i32,
						_tau: &[Self],
						_info: &mut i32)
						-> i32
	{
		unimplemented!();
	}

	fn xgeev(_jobvl: u8,
			 _jobvr: u8,
			 _n: i32,
			 _a: &mut [Self],
			 _lda: i32,
			 _wr: &mut [Self],
			 _wi: &mut [Self],
			 _vl: &mut [Self],
			 _ldvl: i32,
			 _vr: &mut [Self],
			 _ldvr: i32,
			 _work: &mut [Self],
			 _lwork: i32,
			 _info: &mut i32)
	{
		unimplemented!();
	}

	fn xgeev_work_size(_jobvl: u8,
					   _jobvr: u8,
					   _n: i32,
					   _a: &mut [Self],
					   _lda: i32,
					   _wr: &mut [Self],
					   _wi: &mut [Self],
					   _vl: &mut [Self],
					   _ldvl: i32,
					   _vr: &mut [Self],
					   _ldvr: i32,
					   _info: &mut i32)
					   -> i32
	{
		unimplemented!();
	}

	fn xgetrf(m: i32, n: i32, a: &mut [Self], lda: i32, ipiv: &mut [i32], info: &mut i32)
	{
		unsafe {
			ffi::cgetrf_(&m, &n, a.as_mut_ptr() as *mut _, &lda, ipiv.as_mut_ptr(), info);
		}
	}

	fn xgeqrf(_m: i32,
			  _n: i32,
			  _a: &mut [Self],
			  _lda: i32,
			  _tau: &mut [Self],
			  _work: &mut [Self],
			  _lwork: i32,
			  _info: &mut i32)
	{
		unimplemented!();
	}

	fn xgeqrf_work_size(_m: i32,
						_n: i32,
						_a: &mut [Self],
						_lda: i32,
						_tau: &mut [Self],
						_info: &mut i32)
						-> i32
	{
		unimplemented!();
	}

	fn xorgqr(_m: i32,
			  _n: i32,
			  _k: i32,
			  _a: &mut [Self],
			  _lda: i32,
			  _tau: &mut [Self],
			  _work: &mut [Self],
			  _lwork: i32,
			  _info: &mut i32)
	{
		unimplemented!();
	}

	fn xorgqr_work_size(_m: i32,
						_n: i32,
						_k: i32,
						_a: &mut [Self],
						_lda: i32,
						_tau: &mut [Self],
						_info: &mut i32)
						-> i32
	{
		unimplemented!();
	}

	fn xgetri(_n: i32,
			  _a: &mut [Self],
			  _lda: i32,
			  _ipiv: &mut [i32],
			  _work: &mut [Self],
			  _lwork: i32,
			  _info: &mut i32)
	{
		unimplemented!();
	}

	fn xgetri_work_size(_n: i32,
						_a: &mut [Self],
						_lda: i32,
						_ipiv: &mut [i32],
						_info: &mut i32)
						-> i32
	{
		unimplemented!();
	}

	fn xpotrf(_uplo: char, _n: i32, _a: &mut [Self], _lda: i32, _info: &mut i32)
	{
		unimplemented!();
	}

	fn xgetrs(_n: i32,
			  _nrhs: i32,
			  _a: &mut [Self],
			  _lda: i32,
			  _ipiv: &mut [i32],
			  _b: &mut [Self],
			  _ldb: i32,
			  _info: &mut i32)
	{
		unimplemented!();
	}
}