use crate::{
    algebra::{
        linear::{matrix::CholeskyDec, Matrix},
    },
    elementary::Power,
};
use crate::algebra::abstr::{Complex, Real, Scalar};
use crate::algebra::abstr::Zero;

impl<T> Matrix<T>
    where T: Real
{
    /// Decomposes the symmetric, positive definite quadratic matrix A into a
    /// lower triangular matrix L A = L L^T
    ///
    /// # Arguments
    ///
    /// A has to be symmetric and positive definite
    ///
    /// # Panics
    ///
    /// If A is not a quadratic matrix
    ///
    /// # Example
    ///
    /// ```
    /// # #[macro_use]
    /// # extern crate mathru;
    /// # fn main()
    /// # {
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = matrix![   2.0, -1.0, 0.0;
    ///                                -1.0, 2.0, -1.0;
    ///                                 0.0, -1.0,  2.0];
    ///
    /// let l: (Matrix<f64>) = a.dec_cholesky().unwrap().l();
    /// # }
    /// ```
    pub fn dec_cholesky(self: &Self) -> Result<CholeskyDec<T>, ()>
    {
        let (m, n): (usize, usize) = self.dim();
        assert_eq!(m, n);

        let (m, n) = self.dim();
        let mut l: Matrix<T> = Matrix::zero(m, n);

        for i in 0..n
        {
            for j in 0..i + 1
            {
                let mut sum = T::zero();
                for k in 0..j
                {
                    sum += *l.get(i, k) * *l.get(j, k);
                }

                if i == j
                {
                    *l.get_mut(i, j) = (*self.get(i, i) - sum).sqrt();
                }
                else
                {
                    *l.get_mut(i, j) = (*self.get(i, j) - sum) / *l.get(j, j);
                }
            }
        }
        return Ok(CholeskyDec::new(l));
    }
}

impl<T> Matrix<Complex<T>>
    where T: Real, Complex<T>: Scalar
{
    /// Decomposes the symmetric, positive definite quadratic matrix A into a
    /// lower triangular matrix L A = L L^T
    ///
    /// # Arguments
    ///
    /// A has to be symmetric and positive definite
    ///
    /// # Panics
    ///
    /// If A is not a quadratic matrix
    ///
    /// # Example
    ///
    /// ```
    /// # #[macro_use]
    /// # extern crate mathru;
    /// # fn main()
    /// # {
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = matrix![   2.0, -1.0, 0.0;
    ///                                -1.0, 2.0, -1.0;
    ///                                 0.0, -1.0,  2.0];
    ///
    /// let l: (Matrix<f64>) = a.dec_cholesky().unwrap().l();
    /// # }
    /// ```
    pub fn dec_cholesky(self: &Self) -> Result<CholeskyDec<Complex<T>>, ()>
    {
        let (m, n): (usize, usize) = self.dim();
        assert_eq!(m, n);

        let (m, n) = self.dim();
        let mut l: Matrix<Complex<T>> = Matrix::zero(m, n);

        for i in 0..n
        {
            for j in 0..i + 1
            {
                let mut sum = Complex::<T>::zero();
                for k in 0..j
                {
                    sum += *l.get(i, k) * l.get(j, k).conj();
                }

                if i == j
                {
                    *l.get_mut(i, j) = (*self.get(i, i) - sum).sqrt();
                }
                else
                {
                    *l.get_mut(i, j) = (*self.get(i, j) - sum) / *l.get(j, j);
                }
            }
        }
        return Ok(CholeskyDec::new(l));
    }
}
