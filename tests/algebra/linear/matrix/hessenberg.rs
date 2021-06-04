use mathru::algebra::linear::Matrix;

#[cfg(feature = "native")]
use mathru::algebra::abstr::Complex;

#[cfg(feature = "native")]
#[test]
fn dec_f64()
{
    let a: Matrix<f64> = matrix![   1.0, 5.0, 3.0;
                                    1.0, 0.0, -7.0;
                                    3.0, 8.0, 9.0];

    let h_ref: Matrix<f64> = matrix![   1.0, -4.427188724235731, 3.7947331922020537;
                                        -3.162277660168379, 8.4, -5.2;
                                        0.0, 9.8, 0.6];

    let q_ref: Matrix<f64> = matrix![   1.0, 0.0, 0.0;
                                        0.0, -0.316227766016838, 0.9486832980505137;
                                        0.0, -0.9486832980505137, -0.3162277660168381];

    let (q, h): (Matrix<f64>, Matrix<f64>) = a.dec_hessenberg().qh();

    assert_relative_eq!(q, q_ref, epsilon=1.0e-10);
    assert_relative_eq!(h, h_ref, epsilon=1.0e-10);
}

#[cfg(feature = "native")]
#[test]
fn dec_complex_f64()
{
    let a: Matrix<Complex<f64>> = matrix![  Complex::new(1.0, 0.0), Complex::new(5.0, 0.0), Complex::new(3.0, 0.0);
                                            Complex::new(1.0, 0.0), Complex::new(0.0, 0.0), Complex::new(-7.0, 0.0);
                                            Complex::new(3.0, 0.0), Complex::new(8.0, 0.0), Complex::new(9.0, 0.0)];

    let q_ref: Matrix<Complex<f64>> = matrix![  Complex::new(1.0, 0.0), Complex::new(0.0, 0.0), Complex::new(0.0, 0.0);
                                                Complex::new(0.0, 0.0), Complex::new(-0.316227766016838, 0.0), Complex::new(-0.9486832980505137, 0.0);
                                                Complex::new(0.0, 0.0), Complex::new(-0.9486832980505137, 0.0), Complex::new(0.3162277660168381, 0.0)];

    let h_ref: Matrix<Complex<f64>> = matrix![  Complex::new(1.0, 0.0), Complex::new(-4.427188724235731, 0.0), Complex::new(-3.7947331922020537, 0.0);
                                                Complex::new(-3.162277660168379, 0.0), Complex::new(8.4, 0.0), Complex::new(5.2, 0.0);
                                                Complex::new(0.0, 0.0), Complex::new(-9.8, 0.0), Complex::new(0.6, 0.0)];

    let (q, h): (Matrix<Complex<f64>>, Matrix<Complex<f64>>) = a.dec_hessenberg().qh();

    assert_relative_eq!(q, q_ref, epsilon=Complex::new(1.0e-10, 1.0e-10));
    assert_relative_eq!(h, h_ref, epsilon=Complex::new(1.0e-10, 1.0e-10));
}

#[cfg(feature = "lapack")]
#[test]
fn dec_f32()
{
    let a: Matrix<f32> = matrix![   1.0, 5.0, 3.0;
                                    1.0, 0.0, -7.0;
                                    3.0, 8.0, 9.0];

    let h_ref: Matrix<f32> = matrix![   1.0, -4.427188724235731, -3.7947331922020537;
                                        -3.162277660168379, 8.4, 5.2;
                                        0.0, -9.8, 0.6];

    let q_ref: Matrix<f32> = matrix![   1.0, 0.0, 0.0;
                                        0.0, -0.316227766016838, -0.9486832980505137;
                                        0.0, -0.9486832980505137, 0.3162277660168381];

    let (q, h): (Matrix<f32>, Matrix<f32>) = a.dec_hessenberg().qh();

    assert_relative_eq!(q, q_ref, epsilon=1.0e-5);
    assert_relative_eq!(h, h_ref, epsilon=1.0e-5);
}

#[cfg(feature = "lapack")]
#[test]
fn dec_f64()
{
    let a: Matrix<f64> = matrix![   1.0, 5.0, 3.0;
                                    1.0, 0.0, -7.0;
                                    3.0, 8.0, 9.0];

    let h_ref: Matrix<f64> = matrix![   1.0, -4.427188724235731, -3.7947331922020537;
                                        -3.162277660168379, 8.4, 5.2;
                                        0.0, -9.8, 0.6];

    let q_ref: Matrix<f64> = matrix![   1.0, 0.0, 0.0;
                                        0.0, -0.316227766016838, -0.9486832980505137;
                                        0.0, -0.9486832980505137, 0.3162277660168381];

    let (q, h): (Matrix<f64>, Matrix<f64>) = a.dec_hessenberg().qh();

    assert_relative_eq!(q, q_ref, epsilon=1.0e-10);
    assert_relative_eq!(h, h_ref, epsilon=1.0e-10);
}
