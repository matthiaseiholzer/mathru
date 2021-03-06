use crate::algebra::{
    abstr::{Field, Scalar},
    linear::Matrix,
};
use std::ops::Add;

impl<T> Add<Self> for Matrix<T> where T: Field + Scalar
{
    type Output = Matrix<T>;

    /// Adds two matrices
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    ///
    /// let c: Matrix<f64> = a + b;
    /// ```
    fn add(self: Self, rhs: Self) -> Self::Output
    {
        (&self).add(&rhs)
    }
}

///Adds two matrices
impl<'a, 'b, T> Add<&'b Matrix<T>> for &'a Matrix<T> where T: Field + Scalar
{
    type Output = Matrix<T>;

    /// Adds two matrices
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    ///
    /// let c: Matrix<f64> = &b + &a;
    /// ```
    fn add(self: Self, rhs: &'b Matrix<T>) -> Self::Output
    {
        assert_eq!(self.dim(), rhs.dim());
        let (m, n) = self.dim();
        Matrix { m,
                 n,
                 data: self.data
                           .iter()
                           .zip(rhs.data.iter())
                           .map(|(x, y)| *x + *y)
                           .collect::<Vec<T>>() }
    }
}

///
/// Add scalar to matrix
impl<'a, 'b, T> Add<&'b T> for &'a Matrix<T> where T: Field + Scalar
{
    type Output = Matrix<T>;

    /// Add a scalar to the matrix
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b: Matrix<f64> = &a + &-4.0;
    /// ```
    fn add(self: Self, rhs: &T) -> Self::Output
    {
        return self.apply(&|x: &T| -> T { x.clone() + rhs.clone() });
    }
}

///
/// Add scalar to matrix
impl<T> Add<T> for Matrix<T> where T: Field + Scalar
{
    type Output = Matrix<T>;

    /// Add a scalar to the matrix
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b: Matrix<f64> = a + -4.0;
    /// ```
    fn add(self: Self, rhs: T) -> Self::Output
    {
        return (&self).add(&rhs);
    }
}
