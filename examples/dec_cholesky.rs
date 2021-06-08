use mathru::{algebra::linear::Matrix, matrix};

fn main()
{
    let a: Matrix<f64> = matrix![   2.0, -1.0, 0.0;
                                	-1.0, 2.0, -1.0;
                                	0.0, -1.0,  2.0];

    let _l: Matrix<f64> = a.dec_cholesky().unwrap().l();
}
