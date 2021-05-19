//! Matrix
use super::{
    //MatrixColumnIterator,
    //MatrixColumnIteratorMut,
    MatrixColumnIntoIterator,
    MatrixIntoIterator,
    MatrixIterator,
    MatrixIteratorMut,
    //MatrixRowIterator,
    //MatrixRowIteratorMut,
    MatrixRowIntoIterator
};
use crate::{
    algebra::{
        abstr::{Addition, Field, Identity, Multiplication, Scalar},
        linear::{matrix::Substitute, Vector, matrix::Transpose},
        abstr::{AbsDiffEq, RelativeEq},
    },
    elementary::Power,
};
use rand::{self, Rng};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::clone::Clone;
use std::{fmt, fmt::Display};
use std::convert::From;

/// Macro to construct matrices
///
/// ```
/// # #[macro_use]
/// # extern crate mathru;
/// # fn main()
/// # {
/// use mathru::algebra::linear::Matrix;
///
/// // Construct a 2x3 matrix of f32
/// let mat: Matrix<f32> = matrix![1.0, 2.0, 3.0; 4.0, 5.0, 6.0];
/// let (m, n) = mat.dim();
///
/// assert_eq!(m, 2);
/// assert_eq!(n, 3);
/// # }
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

#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Debug, Clone)]
pub struct Matrix<T>
{
    /// Num of rows which the matrix has
    pub(super) m: usize,
    /// Num of columns which the matrix ha
    pub(super) n: usize,
    /// Matrix entries
    pub(super) data: Vec<T>,
}

impl<T> From<Vector<T>> for Matrix<T> where T: Field + Scalar
{
    fn from(v: Vector<T>) -> Self
    {
        let (v_m, v_n): (usize, usize) = v.dim();

        return Matrix::new(v_m, v_n, v.convert_to_vec());
    }
}

impl<T> IntoIterator for Matrix<T> where T: Field + Scalar
{
    type Item = T;
    type IntoIter = MatrixIntoIterator<T>;

    fn into_iter(self) -> Self::IntoIter
    {
        MatrixIntoIterator { iter: self.data.into_iter() }
    }
}

impl<T> Matrix<T>
{
    pub fn iter(self: &Self) -> MatrixIterator<T>
    {
        MatrixIterator::new( self.data.iter() )
    }

    pub fn iter_mut(self: &mut Self) -> MatrixIteratorMut<T>
    {
        MatrixIteratorMut::new(self.data.iter_mut())
    }

    pub fn row_into_iter(self: &Self) -> MatrixRowIntoIterator<T>
    {
        MatrixRowIntoIterator::new(self)
    }
    // pub fn row_iter(self: &Self) -> MatrixRowIterator<T>
    //     where T: Zero
    // {
    //     MatrixRowIterator::new(self.data.iter())
    // }

    // pub fn row_iter_mut(self: &mut Self) -> MatrixRowIteratorMut<T>
    //     where T: Zero
    // {
    //     MatrixRowIteratorMut::new(self.data.iter_mut())
    // }

    pub fn column_into_iter(self: &Self) -> MatrixColumnIntoIterator<T>
    {
         MatrixColumnIntoIterator::new(self)
    }

    // pub fn column_iter(self: &Self) -> MatrixColumnIterator<T>
    // {
    //     MatrixColumnIterator::new(self.data.iter())
    // }

    // pub fn column_iter_mut(self: &mut Self) -> MatrixColumnIteratorMut<T>
    // {
    //     MatrixColumnIteratorMut::new(self.data.iter_mut())
    // }
}

impl<T> Matrix<T> where T: Clone
{
    /// Applies the function f on every element in the matrix
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

impl<T> Matrix<T> where T: Field + Scalar + Power
{
    /// Calculates the determinant
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, -2.0, 3.0, -7.0]);
    /// let det: f64 = a.det();
    /// ```
    pub fn det(self: &Self) -> T
    {
        assert_eq!(self.m, self.n);

        if self.m == 1
        {
            return *self.get(0, 0);
        }

        if self.m == 2
        {
            let a_11: T = *self.get(0, 0);
            let a_12: T = *self.get(0, 1);
            let a_21: T = *self.get(1, 0);
            let a_22: T = *self.get(1, 1);
            return a_11 * a_22 - a_12 * a_21;
        }

        let (_l, u, p) = match self.dec_lu()
        {
            Err(_e) => return T::zero(),
            Ok(dec) => dec.lup(),
        };

        let mut det: T = T::one();

        for i in 0..self.m
        {
            det *= *u.get(i, i);
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
            perm = (-T::one()).pow(T::from_u128(counter as u128 - 1));
        }

        return perm * det;
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
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, -2.0, 3.0, -7.0]);
    /// let tr: f64 = a.trace();
    ///
    /// assert_eq!(-6.0, tr);
    /// ```
    pub fn trace(self: &Self) -> T
        where T: Field + Scalar
    {
        let (m, n): (usize, usize) = self.dim();
        if m != n
        {
            panic!("matrix is not square");
        }

        let mut sum: T = T::zero();
        for i in 0..m
        {
            sum += *self.get(i, i);
        }

        return sum;
    }
}

#[cfg(feature = "native")]
impl<T> Matrix<T> where T: Scalar
{
    pub(super) fn swap_rows<'a>(self: &mut Self, i: usize, j: usize)
    {
        for k in 0..self.n
        {
            let temp: T = *self.get(i, k);
            *(self.get_mut(i, k)) = *self.get(j, k);
            *(self.get_mut(j, k)) = temp;
        }
    }
}

impl<T> Matrix<T> where T: Field + Scalar
{
    // returns column vector
    pub fn get_column(self: &Self, i: usize) -> Vector<T>
    {
        assert!(i < self.n);

        let mut v: Vector<T> = Vector::zero(self.m);

        for k in 0..self.m
        {
            *(v.get_mut(k)) = *(self.get(k, i));
        }

        v
    }

    /// return row vector
    ///
    /// i: row
    pub fn get_row(self: &Self, i: usize) -> Vector<T>
    {
        assert!(i < self.m);

        let mut v: Vector<T> = Vector::zero(self.n);
        v = v.transpose();

        for k in 0..self.n
        {
            *(v.get_mut(k)) = *(self.get(i, k));
        }

        v
    }

    /// set column
    pub fn set_column(self: &mut Self, column: &Vector<T>, i: usize)
    {
        let (m, _n) = column.dim();
        if m != self.m
        {
            panic!("Dimensions do not match");
        }

        for k in 0..self.m
        {
            *(self.get_mut(k, i)) = *column.get(k);
        }
    }

    /// set row
    ///
    /// # Arguments
    /// * 'row'
    /// * 'i'
    ///
    /// # Panics
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

impl<T> Matrix<T> where T: Field + Scalar + Power
{
    ///
    pub fn givens<'d, 'e>(m: usize, i: usize, j: usize, c: T, s: T) -> Self
    {
        if i >= m || j >= m
        {
            panic!("Index out of bounds");
        }

        let mut givens: Matrix<T> = Matrix::one(m);
        *(givens.get_mut(i, i)) = c;
        *(givens.get_mut(j, j)) = c;
        *(givens.get_mut(i, j)) = s;
        *(givens.get_mut(j, i)) = -s;
        givens
    }

    #[cfg(feature = "native")]
    /// function [c,s] = Givens(a,b)
    /// Givens rotation computation
    /// Determines cosine-sine pair (c,s) so that [c s;-s c]'*[a;b] = [r;0]
    /// GVL4: Algorithm 5.1.3
    pub fn givens_cosine_sine_pair(a: T, b: T) -> (T, T)
    {
        let exponent: T = T::from_f64(2.0);
        let exponent_sqrt: T = T::from_f64(0.5);

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
                s = T::one() / (T::one() + tau.pow(exponent)).pow(exponent_sqrt);
                c = s * tau;
            }
            else
            {
                let tau: T = -b / a;
                c = T::one() / (T::one() + tau.pow(exponent)).pow(exponent_sqrt);
                s = c * tau;
            }
        }

        return (c, s);
    }
}

impl<T> Matrix<T> where T: Field + Scalar + Power
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
    pub fn householder(v: &Vector<T>, k: usize) -> Self
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

        let d: Vector<T> = v.get_slice(k, v_m - 1);

        let norm: T = T::from_f64(2.0);

        let alpha: T;

        let d_0: T = *d.get(0);

        if d_0 >= T::zero()
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

        *v.get_mut(0) = (T::from_f64(0.5) * (T::one() - d_0 / alpha)).pow(T::from_f64(0.5));
        let p: T = -alpha * *v.get(0);

        if d_m - 1 >= 1
        {
            let temp: Vector<T> = d.get_slice(1, d_m - 1)
                                   .apply(&|e: &T| -> T { *e / (T::from_f64(2.0) * p) });
            v.set_slice(&temp, 1);
        }

        let mut w: Vector<T> = Vector::zero(v_m);

        w.set_slice(&v, k);

        let ident: Matrix<T> = Matrix::one(v_m);

        let two: T = T::from_f64(2.0);
        let w_dyadp: Matrix<T> = w.dyadp(&w);
        let h: Matrix<T> = &ident - &(&w_dyadp * &two);

        h
    }
}

impl<T> Matrix<T> where T: Field + Scalar + Power
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
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(4,
    ///                                  4,
    ///                                  vec![4.0, 1.0, -2.0, 2.0, 1.0, 2.0, 0.0, -2.0, 0.0, 3.0,
    ///                                       -2.0, 2.0, 2.0, 1.0, -2.0, -1.0]);
    ///
    /// let (u, s, v): (Matrix<f64>, Matrix<f64>, Matrix<f64>) = a.dec_sv();
    /// ```
    pub fn dec_sv(self: &Self) -> (Self, Self, Self)
    {
        let (mut u, mut b, mut v): (Matrix<T>, Matrix<T>, Matrix<T>) = self.householder_bidiag();

        let (_m, n): (usize, usize) = b.dim();
        let max_iterations: usize = 500 * n * n;

        let mut u_k: Matrix<T> = Matrix::one(n);

        for _k in 0..max_iterations
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
                *b.get_mut(l, l) = -*b.get(l, l);
                let mut column_l: Vector<T> = u.get_column(l);
                column_l = &column_l * &-T::one();
                u.set_column(&column_l, l);
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
        return (u, b, v);
    }

    fn msweep(mut u: Matrix<T>,
              mut b: Matrix<T>,
              mut v: Matrix<T>)
              -> (Matrix<T>, Matrix<T>, Matrix<T>)
    {
        let (_m, n): (usize, usize) = b.dim();

        for k in 0..n - 1
        {
            let mut q: Matrix<T> = Matrix::one(n);

            // Construct matrix Q and multiply on the right by Q'.
            // Q annihilates both B(k-1,k+1) and B(k,k+1)
            // but makes B(k+1,k) non-zero.
            let (c_r, s_r, _r_r): (T, T, T) = Matrix::rot(*b.get(k, k), *b.get(k, k + 1));

            *q.get_mut(k, k) = c_r;
            *q.get_mut(k, k + 1) = s_r;
            *q.get_mut(k + 1, k) = -s_r;
            *q.get_mut(k + 1, k + 1) = c_r;

            let q_t: Matrix<T> = q.clone().transpose();
            b = &b * &q_t;
            v = &v * &q_t;


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
        }

        return (u, b, v);
    }

    pub fn rot(f: T, g: T) -> (T, T, T)
    {
        return if f == T::zero()
        {
            (T::zero(), T::one(), g)
        }
        else
        {
            let expo: T = T::from_f64(2.0);
            let sqrt: T = T::from_f64(0.5);
            if f.abs() > g.abs()
            {
                let t: T = g / f;
                let t1: T = (T::one() + t.pow(expo)).pow(sqrt);

                (T::one() / t1, t / t1, f * t1)
            }
            else
            {
                let t: T = f / g;
                let t1: T = (T::one() + t.pow(expo)).pow(sqrt);

                (t / t1, T::one() / t1, g * t1)
            }
        }
    }

    ///
    /// self is an m times n matrix with m >= n
    /// A = UBV^{T}
    /// U \in T^{m \times n}
    /// B \in T^{n \times n}
    /// V \in T^{n \times n}
    pub fn householder_bidiag(self: &Self) -> (Self, Self, Self)
    {
        let (m, n): (usize, usize) = self.dim();
        if m < n
        {
            panic!("Read the API");
        }

        let mut u: Matrix<T> = Matrix::one(m);
        let mut v: Matrix<T> = Matrix::one(n);
        let mut a_i: Matrix<T> = self.clone();

        for i in 0..n - 1
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
                let v_x_trans_slice: Vector<T> = v_x_trans.get_slice(i + 1, n - 1);

                let v_i: Matrix<T> = Matrix::householder(&v_x_trans_slice, 0);

                let mut v_ni: Matrix<T> = Matrix::one(n);
                v_ni = v_ni.set_slice(&v_i, i + 1, i + 1);
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
                    *a_i.get_mut(i, k) = T::zero();
                }
            }
        }
        (u, a_i, v)
    }
}

impl<T> Matrix<T> where T: Field + Scalar
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
    pub fn get_slice(self: &Self,
                     row_s: usize,
                     row_e: usize,
                     column_s: usize,
                     column_e: usize)
                     -> Matrix<T>
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

impl<T> Display for Matrix<T> where T: Display
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
    pub fn get_mut(self: &mut Self, i: usize, j: usize) -> &mut T
    {
        assert!(i < self.m);
        assert!(j < self.n);
        &mut (self.data[j * self.m + i])
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
    pub fn get(self: &Self, i: usize, j: usize) -> &T
    {
        assert!(i < self.m);
        assert!(j < self.n);

        return &self.data[j * self.m + i];
    }
}

impl<T> PartialEq for Matrix<T> where T: PartialEq
{
    /// Checks if two matrices are equal
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    ///
    /// assert_eq!(true, a == b);
    /// ```
    fn eq(self: &Self, other: &Self) -> bool
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

impl<T> Matrix<T> where T: Clone + Copy
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
    pub fn new(m: usize, n: usize, data: Vec<T>) -> Self
    {
        assert_eq!(m * n, data.len());
        Matrix { m, n, data }
    }
}

impl<T> Matrix<T>
{
    pub fn convert_to_vec(self) -> Vec<T>
    {
        return self.data;
    }
}

impl<T> Matrix<T> where T: Scalar + Clone + Copy
{
    pub fn new_random(m: usize, n: usize) -> Matrix<T>
    {
        let mut rng = rand::thread_rng();
        let data: Vec<T> = vec![T::from_f64(rng.gen()); m * n];
        Matrix::new(m, n, data)
    }
}

impl<T> Matrix<T> where T: Field + Scalar
{
    /// Returns the zero matrix(additive neutral element)
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b: Matrix<f64> = &a + &Matrix::zero(2, 2);
    /// ```
    pub fn zero(m: usize, n: usize) -> Self
    {
        return Matrix { m,
                        n,
                        data: vec![T::zero(); m * n] };
    }
}

impl<T> Identity<Addition> for Matrix<T> where T: Identity<Addition>
{
    /// Returns the additive neutral element)
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// ```
    fn id() -> Self
    {
        //        Matrix {
        //            m: m,
        //            n: n,
        //            data: vec![Identity::<Addition>::id(); m * n],
        //        }

        unimplemented!();
    }
}

impl<T> Matrix<T> where T: Field + Scalar
{
    /// Returns the eye matrix(multiplicative neutral element)
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 0.0, 3.0, -7.0]);
    /// let b: Matrix<f64> = &a * &Matrix::one(2);
    /// ```
    pub fn one(size: usize) -> Self
    {
        let mut data: Vec<T> = vec![Identity::<Addition>::id(); size * size];

        for i in 0..size
        {
            data[i * size + i] = Identity::<Multiplication>::id();
        }

        Matrix { m: size,
                 n: size,
                 data }
    }

    pub fn ones(m: usize, n: usize) -> Self
    {
        Matrix { m,
                 n,
                 data: vec![Identity::<Multiplication>::id(); m * n] }
    }
}

impl<T> Matrix<T> where T: Field + Scalar + Power
{
    /// Calculates the pseudo inverse matrix
    ///
    /// A^+ = (A^TA)^-1A^T
    pub fn pinv(self: &Self) -> Matrix<T>
    {
        let r: Matrix<T> = self.dec_qr().r();
        let x: Matrix<T> = r.clone()
                            .transpose()
                            .substitute_forward(self.clone().transpose());
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
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(4, 2, vec![1.0, 0.0, 3.0, 0.0, 1.0, -7.0, 0.5, 0.25]);
    /// let (m, n) = a.dim();
    ///
    /// assert_eq!(4, m);
    /// assert_eq!(2, n);
    /// ```
    pub fn dim(&self) -> (usize, usize)
    {
        return (self.m, self.n);
    }

    /// Returns the number of rows
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::{Matrix};
    ///
    /// let a: Matrix<f64> = Matrix::new(4, 2, vec![1.0, 0.0, 3.0, 0.0, 1.0, -7.0, 0.5, 0.25]);
    /// let m: usize = a.nrows();
    ///
    /// assert_eq!(4, m);
    pub fn nrows(self: &Self) -> usize
    {
        return self.m;
    }

    /// Returns the number of columns
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::Matrix;
    ///
    /// let a: Matrix<f64> = Matrix::new(4, 2, vec![1.0, 0.0, 3.0, 0.0, 1.0, -7.0, 0.5, 0.25]);
    /// let n: usize = a.ncols();
    ///
    /// assert_eq!(2, n);
    /// ```
    pub fn ncols(self: &Self) -> usize
    {
        return self.n;
    }
}

impl<T> AbsDiffEq for Matrix<T>
    where T: Field + Scalar + AbsDiffEq<Epsilon = T>
{
    type Epsilon = T;

    fn default_epsilon() -> T
    {
        T::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Matrix<T>, epsilon: T) -> bool
    {
        for (a, b) in self.iter().zip(other.iter())
        {
            if a.abs_diff_ne(b, epsilon)
            {
                return false;
            }
        }

        return true;
    }
}

impl<T> RelativeEq for Matrix<T>
    where T: Field + Scalar + AbsDiffEq<Epsilon = T> + RelativeEq
{
    fn default_max_relative() -> T
    {
        T::default_max_relative()
    }

    /// A test for equality that uses a relative comparison if the values are far apart.
    fn relative_eq(&self, other: &Matrix<T>, epsilon: Self::Epsilon, max_relative: Self::Epsilon) -> bool
    {
        for (a, b) in self.iter().zip(other.iter())
        {
            if a.relative_ne(b, epsilon, max_relative)
            {
                return false;
            }
        }
        return true;
    }
}