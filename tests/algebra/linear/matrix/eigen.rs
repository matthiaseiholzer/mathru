use mathru::algebra::linear::{Matrix, Vector};
use mathru::algebra::abstr::Complex;
use crate::mathru::algebra::abstr::cast::FromPrimitive;

#[test]
fn eigenvalue_f32()
{
    let a: Matrix<f32> = matrix![   1.0, -3.0, 3.0;
                                    3.0, -5.0,  3.0;
                                    6.0, -6.0,  4.0];

    let eig_ref: Vector<f32> = vector![ 4.0;
                                        -2.0;
                                        -2.0];

    let (value, _vector): (Vector<f32>, Matrix<f32>) = a.dec_eigen().pair();

    assert_relative_eq!(value, eig_ref, epsilon=1.0e-5);
}

#[test]
fn eigenvalue_f64()
{
    let a: Matrix<f64> = matrix![   1.0, -3.0, 3.0;
                                    3.0, -5.0,  3.0;
                                    6.0, -6.0,  4.0];

    let eig_ref: Vector<f64> = vector![ 4.0;
                                        -2.0;
                                        -2.0];

    let (value, _vector): (Vector<f64>, Matrix<f64>) = a.dec_eigen().pair();

    assert_relative_eq!(value, eig_ref, epsilon=1.0e-10);
}

#[cfg(feature = "native")]
#[test]
fn eigenvalue_complex_f32()
{
    let a: Matrix<Complex<f32>> = matrix![  Complex::from_f32(1.0), Complex::from_f32(-3.0), Complex::from_f32(3.0);
                                            Complex::from_f32(3.0), Complex::from_f32(-5.0), Complex::from_f32(3.0);
                                            Complex::from_f32(6.0), Complex::from_f32(-6.0), Complex::from_f32(4.0)];

    let eig_value_ref: Vector<Complex<f32>> = vector![  Complex::from_f32(4.0);
                                                        Complex::from_f32(-2.0);
                                                        Complex::from_f32(-2.0)];

    let (value, _vector): (Vector<Complex<f32>>, Matrix<Complex<f32>>) = a.dec_eigen().pair();

    assert_relative_eq!(value, eig_value_ref, epsilon=Complex::new(2.0e-5, 2.0e-5));
}

#[cfg(feature = "lapack")]
#[test]
fn eigenvalue_complex_f32()
{
    let a: Matrix<Complex<f32>> = matrix![  Complex::from_f32(1.0), Complex::from_f32(-3.0), Complex::from_f32(3.0);
                                            Complex::from_f32(3.0), Complex::from_f32(-5.0), Complex::from_f32(3.0);
                                            Complex::from_f32(6.0), Complex::from_f32(-6.0), Complex::from_f32(4.0)];

    let eig_value_ref: Vector<Complex<f32>> = vector![  Complex::from_f32(-2.0);
                                                        Complex::from_f32(4.0);
                                                        Complex::from_f32(-2.0)];

    let (value, _vector): (Vector<Complex<f32>>, Matrix<Complex<f32>>) = a.dec_eigen().pair();

    assert_relative_eq!(value, eig_value_ref, epsilon=Complex::new(1.0e-5, 1.0e-5));
}

#[test]
fn eigenvalue_complex_f64()
{
    let a: Matrix<Complex<f64>> = matrix![  Complex::from_f64(1.0), Complex::from_f64(-3.0), Complex::from_f64(3.0);
                                            Complex::from_f64(3.0), Complex::from_f64(-5.0), Complex::from_f64(3.0);
                                            Complex::from_f64(6.0), Complex::from_f64(-6.0), Complex::from_f64(4.0)];

    let eig_value_ref: Vector<Complex<f64>> = vector![  Complex::from_f64(4.0);
                                                        Complex::from_f64(-2.0);
                                                        Complex::from_f64(-2.0)];

    // let eig_vector_ref = matrix![Complex::new(-0.4082, 0.0), Complex::new(0.2440, -0.4070), Complex::new(0.2440, 0.4070);
    //                             Complex::new(-0.4082, 0.0), Complex::new(-0.4162, -0.4070), Complex::new(-0.4162, 0.4070);
    //                             Complex::new(-0.8165, 0.0), Complex::new(-0.6602, 0.0), Complex::new(-0.6602, 0.0)];

    let (value, _vector): (Vector<Complex<f64>>, Matrix<Complex<f64>>) = a.dec_eigen().pair();

    assert_relative_eq!(value, eig_value_ref, epsilon=Complex::new(1.0e-10, 1.0e-10));
    // assert_relative_eq!(vector, eig_vector_ref, epsilon=Complex::new(1.0e-10, 1.0e-10));
}


