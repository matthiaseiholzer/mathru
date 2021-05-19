use crate::algebra::{
    linear::{
        matrix::{Transpose, EigenDec, Solve},
        Matrix, Vector,
    },
};
use crate::algebra::abstr::{Field, Scalar, AbsDiffEq};
use crate::elementary::Power;

impl<T> Matrix<T> where T: Field + Scalar + Power + AbsDiffEq<Epsilon = T>
{
    /// Computes the eigenvalues of a real matrix
    ///
    ///
    /// # Arguments
    ///
    ///
    ///
    ///
    /// # Example
    ///
    /// ```
    /// use mathru::algebra::linear::{matrix::EigenDec, Matrix, Vector};
    ///
    /// let a: Matrix<f64> = Matrix::new(3, 3, vec![1.0, -3.0, 3.0, 3.0, -5.0, 3.0, 6.0, -6.0, 4.0]);
    /// let eigen: EigenDec<f64> = a.dec_eigen();
    /// ```
    pub fn dec_eigen(self: Self) -> EigenDec<T>
    {
        let (m, n): (usize, usize) = self.dim();
        assert_eq!(m, n, "Unable to compute the eigen value of a non-square matrix");
        assert_ne!(m, 0, "Unable to compute the eigen value of an empty matrix.");

        let value: Vector<T> = self.eigenvalue_r();
        let vector: Matrix<T> = self.eigenvector_r(&value);

        return EigenDec::new(value, vector);
    }

    pub fn eigenvalue_r(self: &Self) -> Vector<T>
    {
        let (m, _n): (usize, usize) = self.dim();

        let h: Matrix<T> = self.dec_hessenberg().h();

        let (_u, t): (Matrix<T>, Matrix<T>) = h.francis();

        let mut eig: Vector<T> = Vector::zero(m);

        for i in 0..m
        {
            *eig.get_mut(i) = *t.get(i, i);
        }

        return eig;
    }

    fn francis(mut self: Self) -> (Matrix<T>, Matrix<T>)
    {
        let epsilon: T = T::default_epsilon();

        let (m, n): (usize, usize) = self.dim();

        let mut u: Matrix<T> = Matrix::one(m);

        let mut p: usize = n;
        let mut q: usize;

        while p > 2
        {
            q = p - 1;

            // Bulge generating
            let s: T = *self.get(q - 1, q - 1) + *self.get(p - 1, p - 1);
            let t: T = *self.get(q - 1, q - 1) * *self.get(p - 1, p - 1)
                       - *self.get(q - 1, p - 1) * *self.get(p - 1, q - 1);

            // compute first 3 elements of first column of M
            let mut x: T = self.get(0, 0).pow(T::from_f64(2.0))
                           + *self.get(0, 1) * *self.get(1, 0)
                           - s * *self.get(0, 0)
                           + t;
            let mut y: T = *self.get(1, 0) * (*self.get(0, 0) + *self.get(1, 1) - s);
            let mut z: T = *self.get(1, 0) * *self.get(2, 1);

            for k in 0..(p - 2)
            {
                let b: Vector<T> = Vector::new_column(3, vec![x, y, z]);
                let h: Matrix<T> = Matrix::householder(&b, 0);

                //Determine the Householder reflector P with P [x; y; z] = αe1 ;
                {
                    let r: usize = k.max(1);

                    let temp = &h * &self.get_slice(k, k + 2, r - 1, n - 1);
                    self = self.set_slice(&temp, k, r - 1);
                }

                {
                    let h_trans: Matrix<T> = h.transpose();
                    let r: usize = p.min(k + 4);
                    let temp: Matrix<T> = &self.get_slice(0, r - 1, k, k + 2) * &h_trans;
                    self = self.set_slice(&temp, 0, k);

                    let temp1: Matrix<T> = &u.get_slice(0, n - 1, k, k + 2) * &h_trans;

                    u = u.set_slice(&temp1, 0, k);
                }

                x = *self.get(k + 1, k);
                y = *self.get(k + 2, k);
                if k < (p - 3)
                {
                    z = *self.get(k + 3, k);
                }
            }

            // Determine the Givens rotation P with P [x; y]T = αe1 ;
            let (c, s): (T, T) = Matrix::givens_cosine_sine_pair(x, y);
            let g: Matrix<T> = Matrix::givens(2, 0, 1, c, s);

            {
                let temp: Matrix<T> = &g * &self.get_slice(q - 1, p - 1, p - 3, n - 1);
                self = self.set_slice(&temp, q - 1, p - 3);
            }

            {
                let g_trans: Matrix<T> = g.transpose();
                let temp: Matrix<T> = &self.get_slice(0, p - 1, p - 2, p - 1) * &g_trans;
                self = self.set_slice(&temp, 0, p - 2);

                let u_slice = &self.get_slice(0, n - 1, p - 2, p - 1) * &g_trans;
                u = u.set_slice(&u_slice, 0, p - 2);
            }

            // check for convergence
            let m: T = self.get(q - 1, q - 1).abs();
            let n: T = self.get(p - 1, p - 1).abs();
            if self.get(p - 1, q - 1).abs() < epsilon * (m + n)
            {
                *self.get_mut(p - 1, q - 1) = T::zero();
                p = p - 1;
            }
            else
            {
                let k: T = self.get(q - 2, q - 2).abs();
                let l: T = self.get(q - 1, q - 1).abs();
                if self.get(p - 2, q - 2).abs() < epsilon * (k + l)
                {
                    *self.get_mut(p - 2, q - 2) = T::zero();
                    p = p - 2;
                }
            }
            p = p - 1;
        }

        return (u, self);
    }

    pub fn eigenvector_r(self: &Self, value: &Vector<T>) -> Matrix<T>
    {
        let eye: Matrix<T> = Matrix::one(self.m);
        let zero_vector: Vector<T> = Vector::zero(self.m);
        let mut vectors: Matrix<T> = Matrix::zero(self.m, self.m);

        for (c, val) in value.iter().enumerate()
        {
            let diff: Matrix<T> = self - &(&eye * val);
            let vec: Vector<T> = diff.solve(&zero_vector).unwrap();
            vectors.set_column(&vec, c);
        }

        return vectors;
    }
}
