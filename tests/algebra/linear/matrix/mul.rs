
use mathru::algebra::linear::{Matrix, Vector};
use mathru::algebra::abstr::Complex;
use mathru::algebra::abstr::cast::FromPrimitive;

#[test]
fn mul_0()
{
    let size: usize = 23;

    let zero: Matrix<f32> = Matrix::zero(size, size);
    let one: Matrix<f32> = Matrix::one(size);

    let res: Matrix<f32> = zero * one;

    assert_relative_eq!(res, Matrix::zero(size, size));
}


#[test]
fn mul_1()
{
    let a: Matrix<f64> = matrix![   1.0, 2.0, 5.0;
                                    3.0, 4.0, 6.0];

    let b: Matrix<f64> = matrix![   5.0, 8.0;
                                    6.0, 9.0;
                                    7.0, 10.0];

    let reference: Matrix<f64> = matrix![   52.0, 76.0;
                                            81.0, 120.0];

    let res: Matrix<f64> = &a * &b;

    assert_relative_eq!(reference, res);
}

#[test]
fn mul_2()
{
    let a: Matrix<f64> = matrix![1.0, 2.0, 5.0];
    let b: Matrix<f64> = matrix![   5.0, 8.0;
                                    6.0, 9.0;
                                    7.0, 10.0];

    let reference: Matrix<f64> = matrix![52.0, 76.0];

    let res: Matrix<f64> = &a * &b;

    assert_relative_eq!(reference, res);
}

#[test]
fn mul_3()
{
    let a: Matrix<f64> = matrix![   1.0, 2.0, 5.0;
                                    3.0, 4.0, 6.0];
    let b: Matrix<f64> = matrix![   5.0;
                                    6.0;
                                    7.0];
    let reference: Matrix<f64> = matrix![   52.0;
                                            81.0];

    let res: Matrix<f64> = &a * &b;

    assert_relative_eq!(reference, res, epsilon=0.00001, max_relative=1.0e-10);

}

#[test]
fn mul_4()
{
    let a: Matrix<f64> = matrix![   1.0, 2.0, 5.0;
                                    3.0, 4.0, 6.0;
                                    7.0, 8.0, 9.0];

    let b: Matrix<f64> = matrix![   5.0, 8.0;
                                    6.0, 9.0;
                                    7.0, 10.0];

    let reference: Matrix<f64> = matrix![   52.0, 76.0;
                                            81.0, 120.0;
                                            146.0, 218.0];

    let res: Matrix<f64> = &a * &b;

    assert_relative_eq!(reference, res, epsilon=0.00001, max_relative=1.0e-10);
}

#[test]
fn mul_complex_f32()
{
    let a: Matrix<Complex<f32>> = matrix![  Complex::new(1.0, 1.0), Complex::new(-2.0, 2.0) ;
                                            Complex::new(-4.0, 3.0), Complex::new(1.0, -5.0)];

    let b: Matrix<Complex<f32>> = matrix![  Complex::new(-2.0, 3.0), Complex::new(-7.0, 5.0);
                                            Complex::new(-5.0, 2.0), Complex::new(-1.0, -4.0)];

    let mul_ref: Matrix<Complex<f32>> = matrix![Complex::new(1.0, -13.0), Complex::new(-2.0, 4.0);
                                                Complex::new(4.0, 9.0), Complex::new(-8.0, -40.0)];

    assert_relative_eq!(mul_ref, &a * &b);
}

#[test]
fn mul_complex_f64()
{
    let a: Matrix<Complex<f64>> = matrix![  Complex::new(1.0, 1.0), Complex::new(-2.0, 2.0) ;
                                            Complex::new(-4.0, 3.0), Complex::new(1.0, -5.0)];

    let b: Matrix<Complex<f64>> = matrix![  Complex::new(-2.0, 3.0), Complex::new(-7.0, 5.0);
                                            Complex::new(-5.0, 2.0), Complex::new(-1.0, -4.0)];

    let mul_ref: Matrix<Complex<f64>> = matrix![Complex::new(1.0, -13.0), Complex::new(-2.0, 4.0);
                                                Complex::new(4.0, 9.0), Complex::new(-8.0, -40.0)];

    assert_relative_eq!(mul_ref, &a * &b);
}


#[test]
fn mul_scalar_f32()
{
    let a: Matrix<f32> = matrix![   1.0, -2.0, 5.0;
                                    3.0, 4.0, 6.0;
                                    7.0, -8.0, 9.0];
    let reference: Matrix<f32> = matrix![   -2.0, 4.0, -10.0;
                                            -6.0, -8.0, -12.0;
                                            -14.0, 16.0, -18.0];

    let res: Matrix<f32> = &a * &-2.0;

    assert_relative_eq!(reference, res);
}

#[test]
fn mul_scalar_f64()
{
    let a: Matrix<f64> = matrix![   1.0, -2.0, 5.0;
                                    3.0, 4.0, 6.0;
                                    7.0, -8.0, 9.0];
    let reference: Matrix<f64> = matrix![   -2.0, 4.0, -10.0;
                                            -6.0, -8.0, -12.0;
                                            -14.0, 16.0, -18.0];

    let res: Matrix<f64> = &a * &-2.0;

    assert_relative_eq!(reference, res);
}

#[test]
fn mul_scalar_complex_f32()
{
    let a: Matrix<Complex<f32>> = matrix![  Complex::new(1.0, 1.0), Complex::new(-2.0, 2.0);
                                            Complex::new(-4.0, 3.0), Complex::new(1.0, -5.0)];

    let reference: Matrix<Complex<f32>> = matrix![  Complex::new(-4.0, 0.0), Complex::new(0.0, -8.0);
                                                    Complex::new(2.0, -14.0), Complex::new(8.0, 12.0)];

    let res: Matrix<Complex<f32>> = &a * &Complex::new(-2.0, 2.0);

    assert_relative_eq!(reference, res);
}

#[test]
fn matrix_mul_scalar_f32()
{
    let m: Matrix<f32> = matrix![   1.0, 2.0;
                                    3.0, 4.0];
    let prod_ref: Matrix<f32> = matrix![    -0.5, -1.0;
                                            -1.5, -2.0];

    assert_relative_eq!(prod_ref, m * -0.5);
}

#[test]
fn matrix_mul_scalar_f64()
{
    let m: Matrix<f64> = matrix![   1.0, 2.0;
                                    3.0, 4.0];
    let prod_ref: Matrix<f64> = matrix![    -0.5, -1.0;
                                            -1.5, -2.0];

    assert_relative_eq!(prod_ref, m * -0.5);
}

#[test]
fn matrix_mul_vector_f32()
{
    let m: Matrix<f32> = matrix![   1.0, 2.0;
                                    3.0, 4.0];

    let v: Vector<f32> = vector![   2.0;
                                    4.0];
    let prod_ref: Vector<f32> = vector![    10.0;
                                            22.0];

    assert_relative_eq!(prod_ref, m * v);
}

#[test]
fn matrix_mul_vector_f64()
{
    let m: Matrix<f64> = matrix![   1.0, 2.0;
                                    3.0, 4.0];

    let v: Vector<f64> = vector![   2.0;
                                    4.0];
    let prod_ref: Vector<f64> = vector![    10.0;
                                            22.0];

    assert_relative_eq!(prod_ref, m * v);
}

#[test]
fn matrix_mul_vector_complex_f32()
{
    let m: Matrix<Complex<f32>> = matrix![  Complex::from_f64(1.0), Complex::from_f64(2.0);
                                            Complex::from_f64(3.0), Complex::from_f64(4.0)];

    let v: Vector<Complex<f32>> = vector![  Complex::from_f64(2.0);
                                            Complex::from_f64(4.0)];
    let prod_ref: Vector<Complex<f32>> = vector![   Complex::from_f64(10.0);
                                                    Complex::from_f64(22.0)];

    assert_relative_eq!(prod_ref, m * v);
}