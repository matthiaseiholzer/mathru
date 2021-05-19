use mathru::algebra::linear::Matrix;
use crate::mathru::algebra::abstr::cast::FromPrimitive;
use mathru::algebra::abstr::Complex;

#[test]
fn sub()
{
    let dim: usize = 5;
    let m_zero: Matrix<f32> = Matrix::zero(dim, dim);
    let m_one: Matrix<f32> = Matrix::one(dim);

    let m_res: Matrix<f32> = &m_zero - &m_one;

    assert_relative_eq!(m_one * -1.0, m_res);
}

#[test]
fn sub_1()
{
    let a: Matrix<f32> = matrix![   1.0, -2.0, -3.0;
                                    -4.0, -1.0, -2.5];

    let b: Matrix<f32> = matrix![   -2.0, -7.0, 3.0;
                                    -5.0, 3.5, 0.0];

    let sum_ref: Matrix<f32> = matrix![ 3.0, 5.0, -6.0;
                                        1.0, -4.5, -2.5];

    assert_relative_eq!(sum_ref, a - b);
}

#[test]
fn sub_2()
{
    let a: Matrix<f64> = matrix![   1.0, -2.0, -3.0;
                                    -4.0, -1.0, -2.5];

    let b: Matrix<f64> = matrix![   -2.0, -7.0, 3.0;
                                    -5.0, 3.5, 0.0];

    let sum_ref: Matrix<f64> = matrix![ 3.0, 5.0, -6.0;
                                        1.0, -4.5, -2.5];

    assert_relative_eq!(sum_ref, a - b);
}

#[test]
fn sub_3()
{
    let a: Matrix<Complex<f32>> = matrix![   Complex::from_f32(1.0), Complex::from_f32(-2.0), Complex::from_f32(-3.0);
                                    Complex::from_f32(-4.0), Complex::from_f32(-1.0), Complex::from_f32(-2.5)];

    let b: Matrix<Complex<f32>> = matrix![   Complex::from_f32(-2.0), Complex::from_f32(-7.0), Complex::from_f32(3.0);
                                    Complex::from_f32(-5.0), Complex::from_f32(3.5), Complex::from_f32(0.0)];

    let sum_ref: Matrix<Complex<f32>> = matrix![Complex::from_f32(3.0), Complex::from_f32(5.0), Complex::from_f32(-6.0);
                                        Complex::from_f32(1.0), Complex::from_f32(-4.5), Complex::from_f32(-2.5)];

    assert_relative_eq!(sum_ref, a - b);
}

#[test]
fn sub_4()
{
    let a: Matrix<Complex<f64>> = matrix![   Complex::new(1.0, 1.0), Complex::new(-2.0, 2.0) ;
                                    Complex::new(-4.0, 3.0), Complex::new(1.0, -5.0)];

    let b: Matrix<Complex<f64>> = matrix![   Complex::new(-2.0, 3.0), Complex::new(-7.0, 5.0);
                                    Complex::new(-5.0, 2.0), Complex::new(-1.0, -4.0)];

    let sum_ref: Matrix<Complex<f64>> = matrix![Complex::new(3.0, -2.0), Complex::new(5.0, -3.0);
                                        Complex::new(1.0, 1.0), Complex::new(2.0, -1.0)];

    assert_relative_eq!(sum_ref, a - b);
}

