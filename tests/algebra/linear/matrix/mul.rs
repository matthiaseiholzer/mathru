#[cfg(test)]
mod mul
{
    use mathru::algebra::linear::{Matrix, Vector};

    #[test]
    fn mul_1()
    {
        let size: usize = 23;

        let zero: Matrix<f32> = Matrix::zero(size, size);
        let one: Matrix<f32> = Matrix::one(size);

        let res: Matrix<f32> = zero * one;

        assert_eq!(res, Matrix::zero(size, size));
    }

    #[test]
    fn mul_2()
    {
        let size: usize = 23;

        let i1: Matrix<f32> = Matrix::one(size);
        let i2: Matrix<f32> = Matrix::one(size);

        let res: Matrix<f32> = i1 * i2;

        assert_eq!(res, Matrix::one(size));
    }

    #[test]
    fn mul_3()
    {
        let a: Matrix<f64> = matrix![   1.0, 2.0, 5.0;
                                        3.0, 4.0, 6.0];

        let b: Matrix<f64> = matrix![   5.0, 8.0;
                                        6.0, 9.0;
                                        7.0, 10.0];

        let reference: Matrix<f64> = matrix![   52.0, 76.0;
                                                81.0, 120.0];

        let res: Matrix<f64> = &a * &b;

        assert_eq!(reference, res);
    }

    #[test]
    fn mul_4()
    {
        let a: Matrix<f64> = matrix![1.0, 2.0, 5.0];
        let b: Matrix<f64> = matrix![   5.0, 8.0;
                                        6.0, 9.0;
                                        7.0, 10.0];

        let reference: Matrix<f64> = matrix![52.0, 76.0];

        let res: Matrix<f64> = &a * &b;

        assert_eq!(reference, res);
    }

    #[test]
    fn mul_5()
    {
        let a: Matrix<f64> = matrix![1.0, 2.0, 5.0; 3.0, 4.0, 6.0];
        let b: Matrix<f64> = matrix![5.0; 6.0; 7.0];
        let reference: Matrix<f64> = matrix![52.0; 81.0];

        let res: Matrix<f64> = &a * &b;

        assert_eq!(reference, res);
    }

    #[test]
    fn mul_6()
    {
        let a: Matrix<f64> = matrix![1.0, 2.0, 5.0; 3.0, 4.0, 6.0; 7.0, 8.0, 9.0];
        let b: Matrix<f64> = matrix![5.0, 8.0; 6.0, 9.0; 7.0, 10.0];
        let reference: Matrix<f64> = matrix![52.0, 76.0; 81.0, 120.0; 146.0, 218.0];

        let res: Matrix<f64> = &a * &b;

        assert!(res.compare_neighbourhood(&reference, 10e-10));
    }

    #[test]
    fn mul_scalar_0()
    {
        let a: Matrix<f64> = matrix![1.0, -2.0, 5.0; 3.0, 4.0, 6.0; 7.0, -8.0, 9.0];
        let reference = matrix![-2.0, 4.0, -10.0; -6.0, -8.0, -12.0; -14.0, 16.0, -18.0];

        let res: Matrix<f64> = &a * &-2.0;

        assert_eq!(reference, res);
    }

    #[test]
    fn scalar_mul()
    {
        let m = matrix![1.0, 2.0; 3.0, 4.0];
        let prod_ref = matrix![-0.5, -1.0; -1.5, -2.0];

        let res = m * -0.5;

        assert_eq!(prod_ref, res);
    }

    #[test]
    fn vector_mul()
    {
        let m = matrix![1.0, 2.0; 3.0, 4.0];

        let v = vector![2.0; 4.0];
        let prod_ref = vector![10.0; 22.0];

        let res = m * v;

        assert_eq!(prod_ref, res);
    }
}
