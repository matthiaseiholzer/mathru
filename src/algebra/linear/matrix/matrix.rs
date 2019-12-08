/// Matrix

use std::clone::Clone;
use crate::algebra::abstr::{Scalar, Semiring, Zero, One, Sign};
use crate::algebra::linear::Vector;
use crate::algebra::abstr::Real;
use std::fmt::Display;
use std::fmt;
use rand;
use rand::Rng;
use serde::{Deserialize, Serialize};
use super::{MatrixIntoIterator, MatrixIterator, MatrixIteratorMut, MatrixRowIterator, MatrixRowIteratorMut, MatrixColumnIterator,
MatrixColumnIteratorMut};
use std::cmp::min;
use crate::algebra::linear::matrix::{Substitute};

/// Macro to construct matrices
///
/// ```
/// #[macro_use]
/// extern crate mathru;
/// fn main()
/// {
///     use mathru::algebra::linear::Matrix;
///
///     // Construct a 2x3 matrix of f32
///     let mat: Matrix<f32> = matrix![1.0, 2.0, 3.0; 4.0, 5.0, 6.0];
///     let (m, n) = mat.dim();
///
///     assert_eq!(m, 2);
///     assert_eq!(n, 3);
/// }
/// ```
#[macro_export]
macro_rules! matrix
{
    ($( $( $x: expr ),*);*) =>
    {
        {
            let data_nested_array = [ $( [ $($x),* ] ),* ];
            let rows = data_nested_array.len();
            let cols = data_nested_array[0].len();
            let mut data_array: Vec<_> = Vec::with_capacity(rows * cols);
            for j in 0..cols
            {
                for i in 0..rows
                {
                    data_array.push(data_nested_array[i][j]);
                }
            }
            Matrix::new(rows, cols, data_array)
        }
    }
}


#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Matrix<T>
{
    /// Num of rows which the matrix has
    pub(crate) m: usize,
    /// Num of columns which the matrix ha
    pub(crate) n: usize,
    /// Matrix entries
    pub(crate) data:  Vec<T>
}

impl<T> IntoIterator for Matrix<T>
    where T: Real
{
    type Item = T;
    type IntoIter = MatrixIntoIterator<T>;

    fn into_iter(self) -> Self::IntoIter
    {
        MatrixIntoIterator
        {
            iter: self.data.into_iter(),
            //_phantom: PhantomData::default()
        }
    }
}


impl<T> Matrix<T>
{
    pub fn iter(self: &Self) -> MatrixIterator<T>
    {
        MatrixIterator
        {
            iter: self.data.iter()
        }
    }

    pub fn iter_mut(self: &mut Self) -> MatrixIteratorMut<T>
    {
        MatrixIteratorMut
        {
            iter: self.data.iter_mut()
        }
    }

    pub fn row_iter(self: &Self) -> MatrixRowIterator<T>
    {
        MatrixRowIterator
        {
            iter: self.data.iter()
        }
    }

    pub fn row_iter_mut(self: &mut Self) -> MatrixRowIteratorMut<T>
    {
        MatrixRowIteratorMut
        {
            iter: self.data.iter_mut()
        }
    }

    pub fn column_iter(self: &Self) -> MatrixColumnIterator<T>
    {
        MatrixColumnIterator
        {
            iter: self.data.iter()
        }
    }

    pub fn column_iter_mut(self: &mut Self) -> MatrixColumnIteratorMut<T>
    {
        MatrixColumnIteratorMut
        {
            iter: self.data.iter_mut()
        }
    }
}


impl<T> Matrix<T>
    where T: Clone
{
    /// Applies the function f on every element in the matrix
    ///
    pub fn apply_mut(mut self: Matrix<T>, f: &dyn Fn(&T) -> T) -> Matrix<T>
    {
        self.data = self.data.iter().map(f).collect::<Vec<T>>();
        self
    }

    pub fn apply(self: &Matrix<T>, f: &dyn Fn(&T) -> T) -> Matrix<T>
    {
        return (self.clone()).apply_mut(f);
    }
}


//impl<T> Matrix<T>
 //   where T: Real
//{
    /// Returns the transposed matrix
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::{Matrix};
    ///
    /// let mut a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, -2.0, 3.0, -7.0]);
    /// a = a.transpose();
    ///
    /// let a_trans: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 3.0, -2.0, -7.0]);
    ///
    /// assert_eq!(a_trans, a);
    /// ```
//    pub fn transpose(self: &Self) -> Self
//    {
//        let mut trans : Matrix<T> = Matrix::zero(self.n, self.m);
//
//        for j in 0..self.n
//        {
//            for i in 0..self.m
//            {
//                *trans.get_mut(j, i) = self.get(i, j).clone()
//            }
//        }
//
//        trans
//    }
//}

impl<T> Matrix<T>
    where T: Semiring + Sign
{

    pub fn gcd(mut m: usize, mut n: usize) -> usize
    {
        while m != 0
        {
            let old_m: usize = m;
            m = n % m;
            n = old_m;
        }
        n
    }

    /// Function to transpose a matrix without allocating memory for the transposed matrix
    ///
    /// catanzaro.name/papers/PPoPP-2014.pdf
    /// TODO
    /// # Example
	///
	/// ```
	/// use mathru::algebra::linear::{Matrix};
	///
	/// let mut uut: Matrix<f64> = Matrix::new(4, 2, vec![1.0, 0.0, 3.0, 0.0, 1.0, -7.0, 0.5, 0.25]);
    /// uut = uut.transpose();
    ///
    /// let refer: Matrix<f64> = Matrix::new(2, 4, vec![1.0, 1.0, 0.0, -7.0, 3.0, 0.5, 0.0, 0.25]);
    /// assert_eq!(refer, uut);
    /// ```
    pub fn transpose(self: Self) -> Matrix<T>
    {
        let (m, n) = self.dim();
        let mut matrix_t: Matrix<T> = Matrix::zero(n, m);

        for i in 0..m
        {
            for j in 0..n
            {
                *matrix_t.get_mut(j, i)  = *self.get(i, j);
            }
        }

        return matrix_t;
    }
//    pub fn transpose_r(mut self: Self) -> Matrix<T>
//    {
//        let gcdiv: usize = Matrix::<T>::gcd(self.m, self.n);
//
//        let a: usize = self.m / gcdiv;
//        let b: usize = self.n / gcdiv;
//
//        let bigger: usize = self.m.max(self.n);
//
//        let mut temp : Vec<T> = Vec::with_capacity(bigger);
//        //Bad
//        unsafe {temp.set_len(bigger)};
//
//        if gcdiv > 1
//        {
//            for j in 0..self.n
//            {
//                for i in 0..self.m
//                {
//                    temp[i]  = *self.get(self.gather(i, j, b), j)
//                    //temp[i] = self.data[self.gather(i, j, b) * self.n + j];
//                }
//
//                for i in 0..self.m
//                {
//                    *self.get_mut(i, j) = temp[i];
//                    //self.data[i * self.n + j] = temp[i];
//                }
//            }
//        }
//
//        for i in 0..self.m
//        {
//            for j in 0..self.n
//            {
//                temp[self.scatter_da(i, j, b)] = *self.get(i, j);
//                //temp[self.scatter_da(i, j, b)] = self.data[i * self.n + j];
//            }
//
//            for j in 0..self.n
//            {
//                *self.get_mut(i, j) = temp[j];
//                //self.data[i * self.n + j] = temp[j];
//            }
//        }
//
//        for j in 0..self.n
//        {
//            for i in 0..self.m
//            {
//                temp[i] = *self.get(self.gather_sa(i, j, a), j);
//                //temp[i] =  self.data[self.gather_sa(i, j, a) * self.n + j];
//            }
//
//            for i in 0..self.m
//            {
//                *self.get_mut(i, j) = temp[i];
//                //self.data[i * self.n + j] = temp[i];
//            }
//        }
//
//        let temp_m : usize = self.m;
//        self.m = self.n;
//        self.n = temp_m;
//
//        return self;
//    }
//
//    fn gather<'a>(self: &'a Self, i: usize, j: usize, b: usize) -> usize
//    {
//        return (i + (j / b)) % self.m;
//    }
//
//    fn gather_sa<'a>(self: &'a Self, i: usize, j: usize, a: usize) -> usize
//    {
//        return (j + i * self.n - i / a) % self.m;
//    }
//
//    fn scatter_da<'a>(self: &'a Self, i: usize, j: usize, b: usize) -> usize
//    {
//        return (((i + j / b) % self.m) + j * self.m) % self.n;
//    }

}

impl<T> Substitute<Vector<T>> for Matrix<T>
    where T: Real
{

    fn substitute_forward(self: &Self, b: Vector<T>) -> Vector<T>
    {
        return self.substitute_forward_vector_r(b);
    }

    fn substitute_backward(self: &Self, b: Vector<T>) -> Vector<T>
    {
        return self.substitute_backward_vector_r(b);
    }
}

impl<T> Substitute<Matrix<T>> for Matrix<T>
    where T: Real
{

    fn substitute_forward(self: &Self, b: Matrix<T>) -> Matrix<T>
    {
        return self.substitute_forward_matrix_r(b);
    }

    fn substitute_backward(self: &Self, b: Matrix<T>) -> Matrix<T>
    {
        return self.substitute_backward_matrix_r(b);
    }
}

impl<T> Matrix<T>
    where T: Real
{

    #[cfg(feature = "native")]
    pub fn substitute_forward_vector_r(self: &Self, mut b: Vector<T>) -> Vector<T>
    {
        for k in 0..self.n
        {
            for l in 0..k
            {
                *b.get_mut(k) = *b.get(k) - *self.get(k, l) * *b.get(l);
            }
            *b.get_mut(k) = *b.get(k) / *self.get(k, k);
        }

        return b;
    }


    #[cfg(feature = "native")]
    pub fn substitute_forward_matrix_r(self: &Self, mut b: Matrix<T>) -> Matrix<T>
    {
        let min: usize = min(self.m, self.n);
        for k in 0..min
        {
            for l in 0..k
            {
                b.set_row(&(b.get_row(k) - (b.get_row(l) * *self.get(k, l))), k);
            }
            b.set_row(&(b.get_row(k) / *self.get(k, k)), k);
        }

        return b;
    }

    #[cfg(feature = "native")]
    pub fn substitute_backward_vector_r(self: &Self, mut b: Vector<T>) -> Vector<T>
    {
        for k in (0..self.n).rev()
        {
            for l in (k+1..self.n).rev()
            {
                *b.get_mut(k) = *b.get(k) - *self.get(k, l) * *b.get(l);
            }
            *b.get_mut(k) = *b.get(k) / *self.get(k, k);
        }

        return b;
    }

    #[cfg(feature = "native")]
    pub fn substitute_backward_matrix_r(self: &Self, mut b: Matrix<T>) -> Matrix<T>
    {
        let min = min(self.m, self.n);

        for k in (0..min).rev()
        {
            for l in (k+1..min).rev()
            {
                b.set_row(&(b.get_row(k) - (b.get_row(l) * *self.get(k, l))), k);
            }
            b.set_row( &(b.get_row(k) / *self.get(k, k)), k);
        }

        return b;
    }


    #[cfg(feature = "blaslapack")]
    pub fn substitute_forward_matrix_r(self: &Self, mut b: Matrix<T>) -> Matrix<T>
    {
        T::xtrsm('L', 'L', 'N', 'N', b.m as i32, b.n as i32, T::one(), self.data.as_slice(), self.m as i32, b.data.as_mut_slice(), b
        .m as i32);
        return b;
    }

    #[cfg(feature = "blaslapack")]
    pub fn substitute_forward_vector_r(self: &Self, b: Vector<T>) -> Vector<T>
    {
        let (b_m, b_n): (usize, usize) = b.dim();
        let mut b_data = b.convert_to_vec();
        T::xtrsm('L', 'L', 'N', 'N', b_m as i32, b_n as i32, T::one(), self.data.as_slice(), self.m as i32, b_data
        .as_mut_slice(), b_m as i32);
        return Vector::new_column(b_m, b_data);
    }

    #[cfg(feature = "blaslapack")]
    pub fn substitute_backward_matrix_r(self: &Self, mut b: Matrix<T>) -> Matrix<T>
    {
        T::xtrsm('L', 'U', 'N', 'N', b.m as i32, b.n as i32, T::one(), self.data.as_slice(), self.m as i32, b.data.as_mut_slice(), b
        .m as i32);
        return b;
    }

    #[cfg(feature = "blaslapack")]
    pub fn substitute_backward_vector_r(self: &Self, b: Vector<T>) -> Vector<T>
    {
        let (b_m, b_n): (usize, usize) = b.dim();
        let mut b_data = b.convert_to_vec();
        T::xtrsm('L', 'U', 'N', 'N', b_m as i32, b_n as i32, T::one(), self.data.as_slice(), self.m as i32, b_data
        .as_mut_slice(), b_m as i32);
        return Vector::new_column(b_m, b_data);
    }
}


impl<T> Matrix<T>
    where T: Real
{

    /// Calculates the determinant
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::{Matrix};
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, -2.0, 3.0, -7.0]);
    /// let determinant_a: f64 = a.det();
    ///
    /// assert_eq!(-1.0, determinant_a);
    /// ```
    pub fn det<'a>(self: &'a Self) -> T
    {
        assert_eq!(self.m, self.n);

        if self.m == 1
        {
            return self.get(0, 0).clone();
        }

        if self.m == 2
        {
            let a_11: T = self.get(0, 0).clone();
            let a_12: T = self.get(0, 1).clone();
            let a_21: T = self.get(1, 0).clone();
            let a_22: T = self.get(1, 1).clone();
            return a_11 * a_22 - a_12 * a_21;
        }

        let (_l, u, p): (Matrix<T>, Matrix<T>, Matrix<T>) = self.dec_lu().lup();
        let mut det: T = T::one();

        for i in 0..self.m
        {
			det *= u.get(i, i).clone();
        }

        let mut counter: usize = 0;
        for i in 0..self.m
        {
            if *p.get(i, i) != T::one()
            {
                counter += 1;
            }
        }

        let mut perm: T = T::one();
        if counter != 0
        {
            perm = (-T::one()).pow( &T::from_usize(counter - 1).unwrap());
        }

        perm * det
    }


    /// Calculates the trace of a matrix
    ///
    /// # Arguments
    ///
    /// self: square matrix
    ///
    /// # Panics
    ///
    /// if it is not a square matrix
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::{Matrix};
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, -2.0, 3.0, -7.0]);
    /// let tr: f64 = a.trace();
    ///
    /// assert_eq!(-6.0, tr);
    /// ```
    pub fn trace(self: &Self) -> T
        where T: Real
    {
        let (m, n): (usize, usize) = self.dim();
        if m != n
        {
            panic!("matrix is not square");
        }

        let mut sum: T = T::zero();
        for i in 0..m
        {
            sum += self.get(i, i).clone();
        }

        return sum;
    }
}




#[cfg(feature = "native")]
impl<T> Matrix<T>
    where T: Clone
{
    pub(super) fn swap_rows<'a>(self: &'a mut Self, i: usize, j: usize)
    {
        for k in 0..self.n
        {
            let temp: T = self.get(i, k).clone();
            *(self.get_mut(i, k)) = self.get(j, k).clone();
            *(self.get_mut(j, k)) = temp;
        }
    }
}




impl<T> Matrix<T>
    where T: Real
{
    //
    // returns column vector
    pub fn get_column<'a>(self: &'a Self, i: usize) -> Vector<T>
    {
        let mut v: Vector<T> = Vector::zero(self.m);

        for k in 0..self.m
        {
            *(v.get_mut(k)) = *(self.get(k, i));
        }

        v
    }

    ///
    /// return row vector
    ///
    /// i: row
    pub fn get_row<'a>(self: &'a Self, i: usize) -> Vector<T>
    {
        let mut v: Vector<T> = Vector::zero(self.n);
        v = v.transpose();

        for k in 0..self.n
        {
            *(v.get_mut(k)) = *(self.get(i, k));
        }

        v
    }


    /// set column
    pub fn set_column(mut self: Self, column: &Vector<T>, i: usize) -> Matrix<T>
    {
        let (m, _n) = column.dim();
        if m != self.m
        {
            panic!("Dimensions do not match");
        }

        for k in 0..self.m
        {
            *self.get_mut(k, i) = *column.get(k);
        }
        self
    }

    /// set row
    ///
    /// # Arguments
    /// * 'row'
    /// * 'i'
    ///
    /// # Panics
    ///
    ///
    pub fn set_row(self: &mut Self, row: &Vector<T>, i: usize)
    {
        let (_m, n): (usize, usize) = row.dim();
        if n != self.n
        {
            panic!("Dimensions do not match");
        }

        for k in 0..self.n
        {
            *(self.get_mut(i, k)) = *row.get(k);
        }
    }


}

impl<T> Matrix<T>
    where T: Real
{

    ///
    ///
    pub fn givens<'d, 'e>(m: usize, i: usize, j: usize, c: T, s: T) -> Self
    {
        if i >= m || j >= m
        {
            panic!("Index out of bounds");
        }

        let mut givens: Matrix<T> = Matrix::one(m);
        *(givens.get_mut(i,i)) = c.clone();
        *(givens.get_mut(j,j)) = c;
        *(givens.get_mut(i,j)) = s.clone();
        *(givens.get_mut(j,i)) = -s;
        givens
    }

    #[cfg(feature = "native")]
    /// function [c,s] = Givens(a,b)
    /// Givens rotation computation
    /// Determines cosine-sine pair (c,s) so that [c s;-s c]'*[a;b] = [r;0]
    /// GVL4: Algorithm 5.1.3
    pub(super) fn givens_cosine_sine_pair(a: T,b: T) -> (T, T)
    {
        let exponent: T = T::from_f64(2.0).unwrap();
        let exponent_sqrt: T = T::from_f64(0.5).unwrap();

        let c: T;
        let s: T;

        if b == T::zero()
        {
            c = T::one();
            s = T::zero();
        }
        else
        {
            if b.abs() > a.abs()
            {
                let tau: T = -a / b;
                s = T::one() / (T::one() + tau.pow(&exponent)).pow(&exponent_sqrt);
                c = s * tau;
            }
            else
            {
                let tau: T = -b / a;
			    c = T::one() / (T::one() + tau.pow(&exponent)).pow(&exponent_sqrt);
                s = c * tau;
            }
        }

        return (c, s);
    }
}



impl<T> Matrix<T>
    where T: Real
{
    /// Returns the householder matrix
    ///
    /// # Arguments
    ///
    /// v: Column vector
    /// k: index 0 <= k < m
    ///
    /// # Panics
    ///
    /// if index out of bounds
    pub fn householder(v: & Vector<T>, k: usize) -> Self
    {
        let (v_m, v_n): (usize, usize) = v.dim();
        if k >= v_m
        {
            panic!("Index k out of bounds");
        }

        if v_m == 0
        {
            panic!();
        }

        if v_m == 1
        {
            return Matrix::one(v_m);
        }

        let d: Vector<T> = v.get_slice(k, v_m -1);

        let norm: T = T::from_f64(2.0).unwrap();

        let alpha: T;

        let d_0: T = *d.get(0);

        if  d_0 >= T::zero()
        {
            alpha = -d.p_norm(&norm);
        }
        else
        {
            alpha = d.p_norm(&norm);
        }

        if alpha == T::zero()
        {
            let h: Matrix<T> = Matrix::one(v_n);
            return h;
        }

        let (d_m, _d_n) = d.dim();

        let mut v: Vector<T> = Vector::zero(d_m);

        * v.get_mut(0) = (T::from_f64(0.5).unwrap() * (T::one() - d_0 / alpha)).pow(&T::from_f64(0.5).unwrap());
        let p: T = -alpha * *v.get(0);


        if d_m - 1 >= 1
        {
            let temp: Vector<T> = d.get_slice(1, d_m - 1).apply(&|e: &T| -> T {*e / (T::from_f64(2.0).unwrap() * p)});
            v.set_slice(&temp, 1);
        }

        let mut w: Vector<T> = Vector::zero(v_m);

        w.set_slice(&v, k);

        let ident : Matrix<T> = Matrix::one(v_m);

        let two: T = T::from_f64(2.0).unwrap();
        let w_dyadp: Matrix<T> = w.dyadp(&w);
        let h : Matrix<T> = &ident - &(&w_dyadp * &two);

        h
    }
}

impl<T> Matrix<T>
    where T: Real
{
    /// Computes the singular value decomposition
    ///
    /// M = U * S * V*
    ///
    /// # Return
    ///
    /// (u, s, v)
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::{Matrix};
    ///
    /// let a: Matrix<f64> = Matrix::new(4, 4, vec![4.0, 1.0, -2.0, 2.0, 1.0, 2.0, 0.0, -2.0, 0.0, 3.0, -2.0, 2.0,
    /// 2.0, 1.0, -2.0, -1.0]);
    ///
    /// let (u, s, v): (Matrix<f64>, Matrix<f64>, Matrix<f64>) = a.dec_sv();
    /// ```
    pub fn dec_sv<'a>(self: &'a Self) -> (Self, Self, Self)
    {
        let (mut u, mut b, mut v): (Matrix<T>, Matrix<T>, Matrix<T>) = self.householder_bidiag();

        let (_m, n): (usize, usize) = b.dim();
        let max_iterations: usize = 500 * n * n;


        let mut u_k: Matrix<T> = Matrix::one(n);

        for _k in 0.. max_iterations
        {
            let (u_ks, b_k, v_k): (Matrix<T>, Matrix<T>, Matrix<T>) = Matrix::msweep(u_k, b, v);
            u_k = u_ks;
            b = b_k;
            v = v_k;
        }

        u = u * u_k.transpose();

        let (b_m, _b_n): (usize, usize) = b.dim();

        // check that all singular values are positive
        for l in 0..b_m
        {
            if b.get(l, l) < &T::zero()
            {
                *b.get_mut(l, l) = - *b.get(l, l);
                let mut column_l: Vector<T> = u.get_column(l);
                column_l = &column_l * &-T::one();
                u = u.set_column(&column_l, l);
            }
        }

        // null all values beneath the diagonal
        for l in 0..b_m
        {
            for k in 0..b_m
            {
                if k != l
                {
                    *b.get_mut(k, l) = T::zero();
                }
            }
        }

        // sort singular values in descending order
        return (u, b, v)
    }

    fn msweep(mut u: Matrix<T>, mut b: Matrix<T>, mut v: Matrix<T>) -> (Matrix<T>, Matrix<T>, Matrix<T>)
    {
        let (_m, n): (usize, usize) = b.dim();

        for k in 0..n-1
        {
            let mut q : Matrix<T> = Matrix::one(n);

            // Construct matrix Q and multiply on the right by Q'.
            // Q annihilates both B(k-1,k+1) and B(k,k+1)
            // but makes B(k+1,k) non-zero.
            let (c_r, s_r, _r_r): (T, T, T) = Matrix::rot(*b.get(k, k), *b.get(k, k + 1));

            *q.get_mut(k, k) = c_r;
            *q.get_mut(k, k +1) = s_r;
            *q.get_mut(k + 1, k) = -s_r;
            *q.get_mut(k + 1, k + 1) = c_r;

            let q_t: Matrix<T>= q.clone().transpose();
            b = &b * &q_t;
            v = &v * &q_t;

            //B(find(abs(B)<1.e-13))=0;

            // Construct matrix Q and multiply on the left by Q.
            // Q annihilates B(k+1,k) but makes B(k,k+1) and
            // B(k,k+2) non-zero.
            let (c_l, s_l, _r_l): (T, T, T) = Matrix::rot(*b.get(k, k), *b.get(k + 1, k));

            *q.get_mut(k, k) = c_l;
            *q.get_mut(k, k + 1) = s_l;
            *q.get_mut(k + 1, k) = -s_l;
            *q.get_mut(k + 1, k + 1) = c_l;

            b = &q * &b;
		    u = &q * &u;
            //B(find(abs(B)<1.e-13))=0;
        }

        return (u, b, v);
    }


    pub fn rot(f: T, g: T) -> (T, T, T)
    {
        if f == T::zero()
        {
            return (T::zero(), T::one(), g);
        }
        else
        {
            let expo: T = T::from_f64(2.0).unwrap();
            let sqrt: T = T::from_f64(0.5).unwrap();
            if f.abs() > g.abs()
            {
                let t: T = g / f;
                let t1: T = (T::one() + t.pow(& expo)).pow(& sqrt);

                return (T::one() / t1, t / t1, f * t1);
            }
            else
            {
                let t: T = f / g;
                let t1: T = (T::one() + t.pow(& expo)).pow(& sqrt);

                return (t / t1, T::one() / t1, g * t1);
            }
        }
    }



    ///
    /// self is an m times n matrix with m >= n
    /// A = UBV^{T}
    /// U \in T^{m \times n}
    /// B \in T^{n \times n}
    /// V \in T^{n \times n}
    ///
    pub fn householder_bidiag<'a>(self: &'a Self) -> (Self, Self, Self)
    {
        let (m, n) : (usize, usize) = self.dim();
        if m < n
        {
            panic!("Read the API");
        }

        let mut u: Matrix<T> = Matrix::one(m);
        let mut v: Matrix<T> = Matrix::one(n);
        //let mut b_i : Matrix<T> = Matrix::zero(&n, &n);
        let mut a_i : Matrix<T> = self.clone();


        for i in 0..n-1
        {
            // eliminate non-zeros below the diagonal
            // Keep the product U*B unchanged
            let u_x: Vector<T> = a_i.clone().get_column(i);
            let u_slice: Vector<T> = u_x.get_slice(i, m - 1);

            let u_i: Matrix<T> = Matrix::householder(&u_slice, 0);

            let a_i_slice = &u_i * &a_i.clone().get_slice(i, m - 1, i, n - 1);
            a_i = a_i.set_slice(&a_i_slice, i, i);
            let mut u_mi: Matrix<T> = Matrix::one(m);
            u_mi = u_mi.set_slice(&u_i, i, i);


            u = &u * &u_mi;

            //eliminate non-zeros to the right of the
            //superdiagonal by working with the transpose
            // Keep the product B*V' unchanged
            //B_T = B';
            if i < (n - 1)
            {
                let v_x: Vector<T> = a_i.get_row(i);
                let v_x_trans: Vector<T> = v_x.transpose();
                let v_x_trans_slice: Vector<T> = v_x_trans.get_slice(i+1, n - 1);

                let v_i: Matrix<T> = Matrix::householder(&v_x_trans_slice, 0);

                let mut v_ni: Matrix<T> = Matrix::one(n);
                v_ni = v_ni.set_slice(&v_i, i+1, i+1);
                //let a_i_slice = &a_i.clone().get_slice(i+1, m - 1, i+1, n - 1) * &v_i;
                //a_i = a_i.set_slice(&a_i_slice, i+1, i+1);
                a_i = &a_i * &v_ni;

                v = &v * &v_ni;
            }

        }

        //Null all elements beneath the diagonal, and superdiagonal
        for i in 0..m
        {
            for k in 0..n
            {
                if k != i && k != (i + 1)
                {
                    *a_i.get_mut(i, k)  = T::zero();
                }
            }
        }
        (u, a_i, v)

    }

}

impl<T> Matrix<T>
    where T: Real
{
    /// Returns a slice of the matrix
    ///
    /// # Arugments
    ///
    /// 0 <= row_s < m \
    /// 0 <= row_e < m \
    /// 0 <= column_s < n \
    /// 0 <= column_e <= n \
    ///
    /// row_s: start row \
    /// row_e: end row \
    /// column_s: start column \
    /// column_e: end column \
    ///
    /// # Example
    ///
    /// ```
    /// # #[macro_use]
    /// # extern crate mathru;
    /// # fn main()
    /// # {
    /// use mathru::algebra::linear::{Matrix};
    ///
    /// let mut a: Matrix<f64> = matrix![1.0, -2.0; 3.0, -7.0];
    /// a = a.get_slice(0, 0, 1, 1);
    ///
    /// let a_ref: Matrix<f64> = Matrix::new(1, 1, vec![-2.0]);
    ///
    /// assert_eq!(a_ref, a);
    /// # }
    /// ```
    pub fn get_slice(self: &Self, row_s: usize, row_e: usize, column_s: usize, column_e: usize) -> Matrix<T>
    {
        assert!(row_s < self.m);
        assert!(row_e < self.m);
        assert!(column_s < self.n);
        assert!(column_e < self.n);

        let mut slice: Matrix<T> = Matrix::zero(row_e - row_s + 1, column_e - column_s + 1);

        for r in row_s..(row_e + 1)
        {
            for c in column_s..(column_e + 1)
            {
                  *slice.get_mut(r - row_s, c - column_s) = *self.get(r, c)
            }
        }
        return slice;
    }

    /// Replaces parts of the matrix with the given values
    ///
    /// # Arugments
    ///
    /// 0 <= row < m \
    /// 0 <= column < n
    ///
    /// # Example
    ///
    /// ```
    /// # #[macro_use]
    /// # extern crate mathru;
    /// # fn main()
    /// # {
    /// use mathru::algebra::linear::{Matrix};
    ///
    /// let mut a: Matrix<f64> = matrix![   1.0, 0.0;
    ///                                     3.0, -7.0];
    ///
    /// let b: Matrix<f64> = matrix![2.0, -1.0];
    /// a = a.set_slice(&b, 0, 0);
    ///
    /// let a_updated: Matrix<f64> = matrix![   2.0, -1.0;
    ///                                         3.0, -7.0];
    ///
    /// assert_eq!(a_updated, a);
    /// # }
    /// ```
    pub fn set_slice(mut self: Self, slice: &Self, row: usize, column: usize) -> Matrix<T>
    {
        let (s_m, s_n): (usize, usize) = slice.dim();
        let (m, n): (usize, usize) = self.dim();
        assert!(row + s_m <= m);
        assert!(column + s_n <= n);

        for r in row..(row + s_m)
        {
            for c in column..(column + s_n)
            {
                *self.get_mut(r, c) = *slice.get(r - row, c - column);
            }
        }
        self
    }
}


impl<T> Display for Matrix<T>
    where T: Display
{
    fn fmt(self: &Self, f: &mut fmt::Formatter) -> fmt::Result
    {
        write!(f, "\n").unwrap();
        for i in 0..self.m
        {
            for j in 0..self.n
            {
                write!(f, "{} ", self.get(i, j)).unwrap();
            }
            write!(f, "\n").unwrap();
        }
        write!(f, "\n")
    }
}


impl<T> Matrix<T>
{
    /// Returns the mutual element a_ij from the matrix
    ///
    /// # Example
    ///
    /// ```
    /// # #[macro_use]
    /// # extern crate mathru;
    /// # fn main()
    /// # {
    /// use mathru::algebra::linear::{Matrix};
    ///
    /// let mut a: Matrix<f64> = matrix![1.0, 0.0; 3.0, -7.0];
    /// *a.get_mut(1, 0) = -8.0;
    ///
    /// let a_updated: Matrix<f64> = matrix![1.0, 0.0; -8.0, -7.0];
    /// assert_eq!(a_updated, a);
    /// # }
    /// ```
    pub fn get_mut<'a>(self: &'a mut Self, i: usize, j: usize) -> &'a mut T
    {
        assert!(i < self.m);
        assert!(j < self.n);
        & mut(self.data[j * self.m + i])
    }


    /// Returns the element a_ij from the matrix
    ///
    /// # Example
    ///
    /// ```
    /// # #[macro_use]
    /// # extern crate mathru;
    /// # fn main()
    /// # {
    /// use mathru::algebra::linear::{Matrix};
    ///
    /// let a: Matrix<f64> = matrix![   1.0, 0.0;
    ///                                 3.0, -7.0];
    ///
    /// let a_ref: f64 = 3.0;
    /// let element: f64 = *a.get(1, 0);
    ///
    /// assert_eq!(a_ref, element);
    /// # }
    /// ```
    pub fn get(self: &Self, i: usize, j: usize) -> & T
    {
        assert!(i < self.m);
        assert!(j < self.n);

        return & self.data[j * self.m + i];
    }
}


impl<T> PartialEq for Matrix<T>
    where T: PartialEq
{
    /// Checks if two matrices are equal
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::{Matrix};
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    ///
    /// assert_eq!(true, a==b);
    /// ```
    fn eq<'a, 'b>(self: &'a Self, other: &'b Self) -> bool
    {

        if self.dim() != other.dim()
        {
            return false;
        }

        if self.data == other.data
        {
            return true;
        }

        false
    }
}

impl<T> Matrix<T>
    where T: Clone + Copy + Zero + One
{
    /// Creates a new Matrix object
    ///
    /// Fortran like, column wise
    ///
    /// [
    ///   0, 1, 2]
    ///   3, 4, 5,
    ///   6, 7, 8
    /// ] => vec![ 0, 3, 6, 1, 4, 7, 2, 5, 8]
    ///
    pub fn new<'a, 'b>(m: usize, n: usize, data: Vec<T>) -> Self
    {
        assert_eq!(m * n, data.len());
        Matrix{m: m, n: n, data: data}
    }

}

impl<T> Matrix<T>
{
    pub fn convert_to_vec(self) -> Vec<T>
    {
        return self.data
    }
}

impl<T> Matrix<T>
    where T: Scalar + Clone + Copy + Zero + One
{
    pub fn new_random(m: usize, n: usize) -> Matrix<T>
    {
        let mut rng = rand::thread_rng();
        let data: Vec<T> = vec![T::from_f64(rng.gen()).unwrap() ; m * n];
        Matrix::new(m, n, data)
    }
}



impl<T> Matrix<T>
    where T: Zero + Clone
{
    /// Returns the zero matrix(additive neutral element)
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::{Matrix};
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b: Matrix<f64> = &a + &Matrix::zero(2, 2);
    ///
    /// assert_eq!(a, b);
    /// ```
    pub fn zero(m: usize, n: usize) -> Self
    {
        Matrix {
            m: m,
            n: n,
            data: vec![T::zero(); m * n],
        }
    }
}


impl<T> Matrix<T>
    where T: One + Zero + Clone
{
    /// Returns the eye matrix(multiplicative neutral element)
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::{Matrix};
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b: Matrix<f64> = &a * &Matrix::one(2);
    ///
    /// assert_eq!(a, b);
    /// ```
    pub fn one(size: usize) -> Self
    {
        let mut data = vec![T::zero(); size * size];

        for i in 0..size
        {
            data[i * size + i] = T::one();
        }

        Matrix
        {
            m: size,
            n: size,
            data: data,
        }
    }

    pub fn ones(m: usize, n: usize) -> Self
    {
        Matrix
        {
            m: m,
            n: n,
            data: vec![T::one(); m * n]
        }
    }
}


impl<T> Matrix<T>
    where T: Real
{
    /// Calculates the pseudo inverse matrix
    ///
    /// A^+ = (A^TA)^-1A^T
    pub fn pinv(self: &Self) -> Matrix<T>
    {
        let r: Matrix<T> = self.dec_qr().r();
        let x: Matrix<T> = r.clone().transpose().substitute_forward(self.clone().transpose());
        let a_pinv: Matrix<T> = r.substitute_backward(x);
        return a_pinv;
    }
}

impl<T> Matrix<T>
{
    /// Returns the matrix dimension
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::{Matrix};
    ///
    /// let a: Matrix<f64> = Matrix::new(4, 2, vec![1.0, 0.0, 3.0, 0.0, 1.0, -7.0, 0.5, 0.25]);
    /// let (m, n) = a.dim();
    ///
    /// assert_eq!(4, m);
    /// assert_eq!(2, n);
    /// ```
    pub fn dim(&self) -> (usize, usize)
    {
        return (self.m, self.n)
    }

    /// Returns the number of rows
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::{Matrix};
    ///
    /// let a: Matrix<f64> = Matrix::new(4, 2, vec![1.0, 0.0, 3.0, 0.0, 1.0, -7.0, 0.5, 0.25]);
    /// let m: usize = a.nrow();
    ///
    /// assert_eq!(4, m);
    pub fn nrow(&self) -> usize
    {
        return self.m
    }

    /// Returns the number of columns
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::{Matrix};
    ///
    /// let a: Matrix<f64> = Matrix::new(4, 2, vec![1.0, 0.0, 3.0, 0.0, 1.0, -7.0, 0.5, 0.25]);
    /// let n: usize = a.ncol();
    ///
    /// assert_eq!(2, n);
    /// ```
    pub fn ncol(&self) -> usize
    {
        return self.n
    }
}

impl<T> Matrix<T>
    where T: Real
{
    /// Compares to matrices
    ///
    /// Checks if all elememnts in the self matrix are in a epsilon neighbourhood of exp
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::{Matrix};
    ///
    /// let a: Matrix<f64> = Matrix::new(4, 2, vec![1.0, 0.0, 3.0, 0.0, 1.0, -7.0, 0.5, 0.25]);
    ///
    pub fn compare_neighbourhood(self: &Self, exp: &Self, epsilon: T) -> bool
    {
        let (exp_m, exp_n): (usize, usize) = exp.dim();
        let (self_m, self_n): (usize, usize) = self.dim();

        assert!(exp_m == self_m);
        assert!(exp_n == self_n);

        for i in 0..exp_m
        {
            for k in 0..exp_n
            {
                if (*exp.get(i, k) - *self.get(i, k)).abs() > epsilon
                {
                    println!("exp: {}, act: {} exp - act: {}", exp, self, (exp - self));
                    return false;
                }
            }
        }

        return true;
    }
}
