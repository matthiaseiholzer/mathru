use mathru::algebra::linear::Matrix;
use mathru::algebra::abstr::Complex;
use mathru::algebra::linear::matrix::Transpose;

#[test]
fn transpose_f32()
{
    let uut: Matrix<f32> = matrix![ 1.0, 0.0;
                                    3.0, 0.0;
                                    1.0, -7.0;
                                   0.5, 0.25];

    let res: Matrix<f32> = uut.clone().transpose();

    let trans_ref: Matrix<f32> = matrix![   1.0, 3.0, 1.0, 0.5;
                                            0.0, 0.0, -7.0, 0.25];

    assert_relative_eq!(res.clone().transpose(), uut);
    assert_relative_eq!(trans_ref, res);
}

#[test]
fn transpose_f64()
{
    let uut: Matrix<f64> = matrix![ 1.0, 0.0;
                                    3.0, 0.0;
                                    1.0, -7.0;
                                   0.5, 0.25];

    let res: Matrix<f64> = uut.clone().transpose();

    let trans_ref: Matrix<f64> = matrix![   1.0, 3.0, 1.0, 0.5;
                                            0.0, 0.0, -7.0, 0.25];

    assert_relative_eq!(res.clone().transpose(), uut);
    assert_relative_eq!(trans_ref, res);
}

#[test]
fn transpose_square()
{
    let uut: Matrix<f32> = matrix![ 1.0, 0.0, 4.0;
                                    3.0, 2.0, 5.0;
                                    7.0, -1.0, 8.0];

    let uut_t_ref: Matrix<f32> = Matrix::new(3, 3, vec![1.0, 0.0, 4.0, 3.0, 2.0, 5.0, 7.0, -1.0, 8.0]);

    assert_relative_eq!(uut_t_ref, uut.transpose());
}

#[test]
fn transpose_row()
{
    let uut: Matrix<f32> = matrix![ 1.0;
                                    3.0;
                                    1.0;
                                    0.5];

    let uut_t_ref: Matrix<f32> = matrix![1.0, 3.0, 1.0, 0.5];

    assert_relative_eq!(uut_t_ref, uut.transpose());
}

#[test]
fn transpose_column()
{
    let uut: Matrix<f32> = matrix![1.0, 3.0, 1.0, 0.5];

    let uut_t_ref: Matrix<f32> = matrix![   1.0;
                                            3.0;
                                            1.0;
                                            0.5];

    assert_relative_eq!(uut_t_ref, uut.transpose());
}

#[test]
fn transpose_complex_f32()
{
    let uut: Matrix<Complex<f32>> = matrix![    Complex::new(1.0, 0.0);
                                                Complex::new(13.0, 0.0);
                                                Complex::new(11.0, 0.0);
                                                Complex::new(10.5, 0.0)];

    let uut_t_ref: Matrix<Complex<f32>> = matrix![Complex::new(1.0, 0.0), Complex::new(13.0, 0.0), Complex::new(11.0, 0.0), Complex::new(10.5, 0.0)];

    assert_relative_eq!(uut_t_ref, uut.transpose());
}

#[test]
fn transpose_complex_f64()
{
    let uut: Matrix<Complex<f64>> = matrix![    Complex::new(1.0, 0.0);
                                                Complex::new(13.0, 0.0);
                                                Complex::new(11.0, 0.0);
                                                Complex::new(10.5, 0.0)];

    let uut_t_ref: Matrix<Complex<f64>> = matrix![Complex::new(1.0, 0.0), Complex::new(13.0, 0.0), Complex::new(11.0, 0.0), Complex::new(10.5, 0.0)];

    assert_relative_eq!(uut_t_ref, uut.transpose());
}