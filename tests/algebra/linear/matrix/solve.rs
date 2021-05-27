use mathru::algebra::linear::{Matrix, Vector, matrix::Solve};
use mathru::algebra::abstr::Complex;

#[test]
fn solve_matrix_f32()
{
    let a: Matrix<f32> = matrix![   1.0, -2.0, 3.0;
                                    2.0, -5.0, 12.0;
                                    0.0, 2.0, -10.0];

    let x_ref: Matrix<f32> = matrix![   -13.0, 7.0, 4.5;
                                        -10.0, 5.0, 3.0;
                                        -2.0, 1.0, 0.5];

    let id: Matrix<f32> = Matrix::one(3);

    let x: Matrix<f32> = a.solve(&id).unwrap();

    assert_relative_eq!(x, x_ref);
}

#[test]
fn solve_matrix_f64()
{
    let a: Matrix<f64> = matrix![   1.0, -2.0, 3.0;
                                    2.0, -5.0, 12.0;
                                    0.0, 2.0, -10.0];

    let x_ref: Matrix<f64> = matrix![   -13.0, 7.0, 4.5;
                                        -10.0, 5.0, 3.0;
                                        -2.0, 1.0, 0.5];

    let id: Matrix<f64> = Matrix::one(3);

    let x: Matrix<f64> = a.solve(&id).unwrap();

    assert_relative_eq!(x, x_ref);
}

#[test]
fn solve_matrix_complex_f32()
{
    let a: Matrix<Complex<f32>> = matrix![  Complex::new(1.0, 0.0), Complex::new(-2.0, 0.0), Complex::new(3.0, 0.0);
                                            Complex::new(2.0, 0.0), Complex::new(-5.0, 0.0), Complex::new(12.0, 0.0);
                                            Complex::new(0.0, 0.0), Complex::new(2.0, 0.0), Complex::new(-10.0, 0.0)];

    let x_ref: Matrix<Complex<f32>> = matrix![  Complex::new(-13.0, 0.0), Complex::new(7.0, 0.0), Complex::new(4.5, 0.0);
                                                Complex::new(-10.0, 0.0), Complex::new(5.0, 0.0), Complex::new(3.0, 0.0);
                                                Complex::new(-2.0, 0.0), Complex::new(1.0, 0.0), Complex::new(0.5, 0.0)];

    let id: Matrix<Complex<f32>> = Matrix::one(3);

    let x: Matrix<Complex<f32>> = a.solve(&id).unwrap();

    assert_relative_eq!(x, x_ref);
}

#[test]
fn solve_matrix_complex_f64()
{
    let a: Matrix<Complex<f64>> = matrix![  Complex::new(1.0, 0.0), Complex::new(-2.0, 0.0), Complex::new(3.0, 0.0);
                                            Complex::new(2.0, 0.0), Complex::new(-5.0, 0.0), Complex::new(12.0, 0.0);
                                            Complex::new(0.0, 0.0), Complex::new(2.0, 0.0), Complex::new(-10.0, 0.0)];

    let x_ref: Matrix<Complex<f64>> = matrix![  Complex::new(-13.0, 0.0), Complex::new(7.0, 0.0), Complex::new(4.5, 0.0);
                                                Complex::new(-10.0, 0.0), Complex::new(5.0, 0.0), Complex::new(3.0, 0.0);
                                                Complex::new(-2.0, 0.0), Complex::new(1.0, 0.0), Complex::new(0.5, 0.0)];
    let id: Matrix<Complex<f64>> = Matrix::one(3);

    let x: Matrix<Complex<f64>> = a.solve(&id).unwrap();

    assert_relative_eq!(x, x_ref);
}

#[test]
fn solve_vector_f32()
{
    let a: Matrix<f32> = matrix![   6.0, 2.0, -1.0;
                                    -3.0, 5.0, 3.0;
                                    -2.0, 1.0, 3.0];

    let b: Vector<f32> = vector![   48.0;
                                    49.0;
                                    24.0];

    let x: Vector<f32> = a.solve(&b).unwrap();
    let x_ref: Vector<f32> = vector![   7.0;
                                        8.0;
                                        10.0];

    assert_relative_eq!(x, x_ref);
}

#[test]
fn solve_vector_f64()
{
    let a: Matrix<f64> = matrix![   6.0, 2.0, -1.0;
                                    -3.0, 5.0, 3.0;
                                    -2.0, 1.0, 3.0];

    let b: Vector<f64> = vector![   48.0;
                                    49.0;
                                    24.0];

    let x: Vector<f64> = a.solve(&b).unwrap();
    let x_ref: Vector<f64> = vector![   7.0;
                                        8.0;
                                        10.0];

    assert_relative_eq!(x, x_ref);
}

#[test]
fn solve_vector_complex_f32()
{
    let a: Matrix<Complex<f32>> = matrix![  Complex::new(1.0, 0.0), Complex::new(0.0, 1.0), Complex::new(-3.0, 1.0);
                                            Complex::new(2.0, 0.0), Complex::new(1.0, 3.0), Complex::new(-4.0, 2.0);
                                            Complex::new(0.0, 2.0), Complex::new(-2.0, 0.0), Complex::new(-2.0, -3.0)];

    let b: Vector<Complex<f32>> = vector![  Complex::new(-1.0, -1.0);
                                            Complex::new(0.0, 2.0);
                                            Complex::new(-1.0, 1.0)];

    let x: Vector<Complex<f32>> = a.solve(&b).unwrap();
    let x_ref: Vector<Complex<f32>> = vector![  Complex::new(4.0, 0.0);
                                                Complex::new(1.0, 1.0);
                                                Complex::new(1.0, 1.0)];

    assert_relative_eq!(x, x_ref);
}

#[test]
fn solve_vector_complex_f64()
{
    let a: Matrix<Complex<f64>> = matrix![  Complex::new(1.0, 0.0), Complex::new(0.0, 1.0), Complex::new(-3.0, 1.0);
                                            Complex::new(2.0, 0.0), Complex::new(1.0, 3.0), Complex::new(-4.0, 2.0);
                                            Complex::new(0.0, 2.0), Complex::new(-2.0, 0.0), Complex::new(-2.0, -3.0)];

    let b: Vector<Complex<f64>> = vector![  Complex::new(-1.0, -1.0);
                                            Complex::new(0.0, 2.0);
                                            Complex::new(-1.0, 1.0)];

    let x: Vector<Complex<f64>> = a.solve(&b).unwrap();
    let x_ref: Vector<Complex<f64>> = vector![  Complex::new(4.0, 0.0);
                                                Complex::new(1.0, 1.0);
                                                Complex::new(1.0, 1.0)];

    assert_relative_eq!(x, x_ref);
}