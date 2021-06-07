use crate::algebra::{
    abstr::{Field, Scalar},
    linear::{
        matrix::{Inverse, Solve, Substitute},
        Matrix, Vector,
    },
};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::clone::Clone;

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, Clone)]
pub struct LUDec<T>
{
    l: Matrix<T>,
    u: Matrix<T>,
    p: Matrix<T>,
}

impl<T> LUDec<T>
{
    pub(super) fn new(l: Matrix<T>, u: Matrix<T>, p: Matrix<T>) -> LUDec<T>
    {
        return LUDec { l, u, p };
    }

    /// Return l Matrix of LU decomposition
    pub fn l(self: Self) -> Matrix<T>
    {
        return self.l;
    }

    pub fn u(self: Self) -> Matrix<T>
    {
        return self.u;
    }

    pub fn p(self: Self) -> Matrix<T>
    {
        return self.p;
    }

    /// Return l, u, and p matrix of the LU decomposition
    pub fn lup(self: Self) -> (Matrix<T>, Matrix<T>, Matrix<T>)
    {
        return (self.l, self.u, self.p);
    }
}

impl<T> Solve<Vector<T>> for LUDec<T> where T: Field + Scalar
{
    /// Solves Ax = y
    ///  where A \in R^{m * n}, x \in R^n, y \in R^m
    fn solve(self: &Self, rhs: &Vector<T>) -> Result<Vector<T>, ()>
    {
        let b_hat: Vector<T> = &self.p * rhs;

        let y: Vector<T> = self.l.substitute_forward(b_hat);

        let x: Vector<T> = self.u.substitute_backward(y);

        return Ok(x);
    }
}

// TODO
impl<T> Inverse<T> for LUDec<T>
    where T: Field + Scalar
{
    /// Inverse Matrix
    ///
    /// PAX = LUX = I
    /// X = (PA)^-1
    /// X = A^-1P^-1
    /// XP = A^-1
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::{matrix::Inverse, Matrix};
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b_inv: Matrix<f64> = a.inv().unwrap();
    /// ```
    fn inv(self: &Self) -> Result<Matrix<T>, ()>
    {
        let b = Matrix::one(self.p.nrows());
        let x: Matrix<T> = self.solve(&b)?;
        return Ok(x);
    }
}

// TODO
impl<T> Solve<Matrix<T>> for LUDec<T>
    where T: Field + Scalar
{
    fn solve(self: &Self, rhs: &Matrix<T>) -> Result<Matrix<T>, ()>
    {
        let b_hat: Matrix<T> = &self.p * rhs;

        let y: Matrix<T> = self.l.substitute_forward(b_hat);
        let x: Matrix<T> = self.u.substitute_backward(y);

        return Ok(x);
    }
}
