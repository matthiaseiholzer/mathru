use crate::{
    algebra::{
        abstr::{Field, Scalar},
        linear::{matrix::EigenDec, Matrix, Vector, matrix::Transpose},
    },
    elementary::Power,
};

impl<T> Matrix<T> where T: Field + Scalar + Power
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
    pub fn dec_eigen(self: &Self) -> EigenDec<T>
    {
        let (m, n): (usize, usize) = self.dim();
        assert_eq!(m, n, "Unable to compute the eigen value of a non-square matrix");
        assert_ne!(m, 0, "Unable to compute the eigen value of an empty matrix.");

        let (m, n): (usize, usize) = self.dim();

        let mut self_data: Vec<T> = self.clone().transpose().data;
        let n_i32: i32 = n as i32;

        let mut info: i32 = 0;

        let mut w: Vec<T> = vec![T::zero(); n];

        let mut temp1 = [T::zero()];
        let mut temp2 = vec![T::zero(); n * n];

        let lwork = T::xgeev_work_size('N' as u8,
                                       'V' as u8,
                                       n_i32,
                                       &mut self_data[..],
                                       n_i32,
                                       w.as_mut_slice(),
                                       &mut temp1,
                                       n_i32,
                                       &mut temp2,
                                       n_i32,
                                       &mut info);

        let mut work: Vec<T> = vec![T::zero(); lwork as usize];

        T::xgeev('N' as u8,
                 'V' as u8,
                 n_i32,
                 &mut self_data[..],
                 n_i32,
                 w.as_mut_slice(),
                 &mut temp1,
                 1 as i32,
                 &mut temp2,
                 n_i32,
                 &mut work,
                 lwork,
                 &mut info);

        assert_eq!(0, info);

        return EigenDec::new(Vector::new_column(m, w), Matrix::new(n, n, temp2));
    }
}
