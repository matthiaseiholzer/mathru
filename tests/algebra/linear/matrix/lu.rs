use mathru::algebra::linear::{
    matrix::{Inverse, LUDec, Solve, Substitute},
    Matrix, Vector,
};
use mathru::algebra::abstr::Complex;
use crate::mathru::algebra::abstr::Zero;

#[test]
fn decompose_lu_0()
{
    let l_ref: Matrix<f64> = matrix![   1.0, 0.0, 0.0;
                                        0.0, 1.0, 0.0;
                                        0.5, 0.25, 1.0];

    let u_ref: Matrix<f64> = matrix![   2.0, -5.0, 12.0;
                                        0.0, 2.0, -10.0;
                                        0.0, 0.0, -0.5];

    let p_ref: Matrix<f64> = matrix![   0.0, 1.0, 0.0;
                                        0.0, 0.0, 1.0;
                                        1.0, 0.0, 0.0];

    let a: Matrix<f64> = matrix![   1.0, -2.0, 3.0;
                                    2.0, -5.0, 12.0;
                                    0.0, 2.0, -10.0];

    let (l, u, p): (Matrix<f64>, Matrix<f64>, Matrix<f64>) = a.dec_lu().unwrap().lup();

    assert_relative_eq!(l, l_ref, epsilon=1.0e-10);
    assert_relative_eq!(u, u_ref, epsilon=1.0e-10);
    assert_relative_eq!(p, p_ref, epsilon=1.0e-10);

    assert_relative_eq!(&p * &a, &l * &u, epsilon=1.0e-10);
}

#[test]
fn decompose_lu_1()
{
    let a: Matrix<f64> = matrix![   4.0, 1.0, -2.0, 2.0;
                                    1.0, 2.0, 0.0, -2.0;
                                    0.0, 3.0, -2.0, 2.0;
                                    2.0, 1.0, -2.0, -1.0];

    let l_ref: Matrix<f64> = matrix![   1.0, 0.0, 0.0, 0.0;
                                        0.0, 1.0, 0.0, 0.0;
                                        0.25, 0.5833333333333334, 1.0, 0.0;
                                        0.5, 0.16666666666666666, -0.4, 1.0 ];

    let u_ref: Matrix<f64> = matrix![   4.0, 1.0, -2.0, 2.0;
                                        0.0, 3.0, -2.0, 2.0;
                                        0.0, 0.0, 1.6666666666666667, -3.666666666666667;
                                        0.0, 0.0, 0.0, -3.8000000000000003 ];

    let p_ref: Matrix<f64> = matrix![   1.0, 0.0, 0.0, 0.0;
                                        0.0, 0.0, 1.0, 0.0;
                                        0.0, 1.0, 0.0, 0.0;
                                        0.0, 0.0, 0.0, 1.0];

    let (l, u, p): (Matrix<f64>, Matrix<f64>, Matrix<f64>) = a.dec_lu().unwrap().lup();

    assert_relative_eq!(l, l_ref, epsilon=1.0e-10);
    assert_relative_eq!(u, u_ref, epsilon=1.0e-10);
    assert_relative_eq!(p, p_ref, epsilon=1.0e-10);
}

#[test]
fn decompose_lu2()
{
    let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 3.0, -2.0, -7.0]);

    let l_ref: Matrix<f64> = matrix![   1.0, 0.0;
                                        1.0 / 3.0, 1.0];

    let u_ref: Matrix<f64> = matrix![   3.0, -7.0;
                                        0.0, 1.0/3.0];

    let (l, u, _p): (Matrix<f64>, Matrix<f64>, Matrix<f64>) = a.dec_lu().unwrap().lup();

    assert_relative_eq!(l, l_ref, epsilon=1.0e-10);
    assert_relative_eq!(u, u_ref, epsilon=1.0e-10);
}

#[test]
fn decompose_lu3()
{
    let a: Matrix<f64> = Matrix::new(2, 2, vec![1.0, 2.0, -3.0, -7.0]);
    let b: Vector<f64> = vector![1.0; 3.0];
    let x_ref: Vector<f64> = vector![-2.25; 8.5];
    let (l, u, p): (Matrix<f64>, Matrix<f64>, Matrix<f64>) = a.dec_lu().unwrap().lup();

    let b_hat = &p * &b;

    let y = u.substitute_backward(b_hat);

    let x = p * l.substitute_forward(y);

    assert_relative_eq!(x, x_ref, epsilon=1.0e-10);
}

#[test]
fn decompose_complex_f32()
{
    let a: Matrix<Complex<f32>> = matrix![  Complex::new(1.0, 1.0), Complex::new(2.0, 2.0);
                                            Complex::new(3.0, 3.0), Complex::new(-4.0, 4.0)];

    let l_ref: Matrix<Complex<f32>> = matrix![  Complex::new(1.0, 0.0), Complex::zero();
                                                Complex::new(1.0/3.0, 0.0), Complex::new(1.0000, 0.0)];

    let u_ref: Matrix<Complex<f32>>  = matrix![ Complex::new(3.0, 3.0), Complex::new(-4.0, 4.0);
                                                Complex::zero(), Complex::new(10.0 / 3.0, 2.0 / 3.0)];

    let p_ref: Matrix<Complex<f32>> = matrix![  Complex::zero(), Complex::new(1.0, 0.0);
                                                Complex::new(1.0, 0.0), Complex::zero()];

    let (l, u, p): (Matrix<Complex<f32>>, Matrix<Complex<f32>>, Matrix<Complex<f32>>) = a.dec_lu().unwrap().lup();

    assert_relative_eq!(l, l_ref);
    assert_relative_eq!(u, u_ref);
    assert_relative_eq!(p, p_ref);

    assert_relative_eq!(p * l * u, a);
}

#[test]
fn decompose_complex_f64()
{
    let a: Matrix<Complex<f64>> = matrix![  Complex::new(1.0, 1.0), Complex::new(2.0, 2.0);
                                            Complex::new(3.0, 3.0), Complex::new(-4.0, 4.0)];

    let l_ref: Matrix<Complex<f64>> = matrix![  Complex::new(1.0, 0.0), Complex::zero();
                                                Complex::new(1.0/3.0, 0.0), Complex::new(1.0000, 0.0)];

    let u_ref: Matrix<Complex<f64>>  = matrix![ Complex::new(3.0, 3.0), Complex::new(-4.0, 4.0);
                                                Complex::zero(), Complex::new(10.0 / 3.0, 2.0 / 3.0)];

    let p_ref: Matrix<Complex<f64>> = matrix![  Complex::zero(), Complex::new(1.0, 0.0);
                                                Complex::new(1.0, 0.0), Complex::zero()];

    let (l, u, p): (Matrix<Complex<f64>>, Matrix<Complex<f64>>, Matrix<Complex<f64>>) = a.dec_lu().unwrap().lup();

    assert_relative_eq!(l, l_ref);
    assert_relative_eq!(u, u_ref);
    assert_relative_eq!(p, p_ref);

    assert_relative_eq!(p * l * u, a);
}

#[test]
fn solve_0()
{
    let a: Matrix<f64> = matrix![6.0, 2.0, -1.0; -3.0, 5.0, 3.0; -2.0, 1.0, 3.0];
    let b: Vector<f64> = vector![48.0; 49.0; 24.0];

    let lu_dec: LUDec<f64> = a.dec_lu().unwrap();
    let x: Vector<f64> = lu_dec.solve(&b).unwrap();
    let x_ref: Vector<f64> = vector![7.0; 8.0; 10.0];

    assert_relative_eq!(x, x_ref, epsilon=10e-10);
}



#[test]
fn inv_0()
{
    let a: Matrix<f64> = matrix![   1.0, -2.0, 3.0;
                                    2.0, -5.0, 12.0;
                                    0.0, 2.0, -10.0];

    let a_inv_ref: Matrix<f64> = matrix![   -13.0, 7.0, 4.5;
                                            -10.0, 5.0, 3.0;
                                            -2.0, 1.0, 0.5];

    let a_inv: Matrix<f64> = a.dec_lu().unwrap().inv().unwrap();

    assert_relative_eq!(a_inv, a_inv_ref, epsilon=1.0e-10);
}

#[test]
fn inv_1()
{
    let a: Matrix<f64> = matrix![   1.0, 0.0, 2.0;
                                    -1.0, 5.0, 0.0;
                                    0.0, 3.0, -9.0];

    let a_inv_ref: Matrix<f64> = matrix![   0.8823529411764706, -0.11764705882352942, 0.19607843137254904;
                                            0.17647058823529413, 0.17647058823529413, 0.03921568627450981;
                                            0.05882352941176471, 0.05882352941176471, -0.09803921568627452];

    let a_inv: Matrix<f64> = a.dec_lu().unwrap().inv().unwrap();

    assert_relative_eq!(a_inv, a_inv_ref, epsilon=1.0e-10);
}

#[test]
fn inv_2()
{
    let a: Matrix<f64> = matrix![   -1.0, 2.0, 3.0, 4.0, 5.0;
                                    -6.0, -7.0, 8.0, 9.0, 10.0;
                                    -11.0, 12.0, 13.0, 14.0, 15.0;
                                    -16.0, -17.0, -18.0, -19.0, 20.0;
                                    -21.0, 22.0, -23.0, 24.0, 25.0];

    let a_inv_ref: Matrix<f64> = matrix![   0.38478669499836576, -0.03759398496240601, -0.08489293886891143, -0.006578947368421052, -0.005720823798627002;
                                            0.03571428571428603, -0.07142857142857142, 0.03571428571428571, 0.0, -0.000000000000000001734723475976807;
                                            -0.021739130434782705, 0.0, 0.04347826086956519, 0.0, -0.021739130434782608;
                                            -0.024517816279830296, 0.06390977443609022, -0.033671134357633165, -0.02631578947368421, 0.020594965675057208;
                                            0.2953293559986926, -0.03007518796992481, -0.030414351095129147, 0.019736842105263157, -0.004576659038901602];

    let a_inv: Matrix<f64> = a.dec_lu().unwrap().inv().unwrap();

    assert_relative_eq!(a_inv, a_inv_ref, epsilon=1.0e-10);
}
