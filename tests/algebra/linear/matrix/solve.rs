use mathru::algebra::linear::{Matrix, Vector, matrix::Solve};

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