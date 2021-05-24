use crate::mathru::algebra::linear::matrix::Substitute;
use mathru::algebra::linear::{matrix::{Transpose, Solve}, Matrix, Vector};

#[test]
fn macro_0()
{
    //Construct a 2x3 matrix of f32
    let mat: Matrix<f32> = matrix![ 1.0, 2.0, 3.0;
                                    4.0, 5.0, 6.0];

    let mat_ref: Matrix<f32> = Matrix::new(2, 3, vec![1.0, 4.0, 2.0, 5.0, 3.0, 6.0]);

    assert_relative_eq!(mat, mat_ref);
}

#[test]
fn macro_1()
{
    //Construct a 2x3 matrix of f32
    let mat: Matrix<f32> = matrix![1.0, 2.0, 3.0, 4.0, 5.0, 6.0];

    let mat_ref: Matrix<f32> = Matrix::new(1, 6, vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0]);

    assert_relative_eq!(mat, mat_ref);
}

#[cfg(feature = "serde")]
#[test]
fn serde_0()
{
    let mat: Matrix<f64> = matrix![1.0, 2.0; 3.0, 4.0];
    let serialized = serde_json::to_string(&mat).unwrap();

    let mat_s: Matrix<f64> = serde_json::from_str(&serialized).unwrap();

    assert_relative_eq!(mat_s, mat);
}

#[test]
fn zeros()
{
    let rows: usize = 5;
    let cols: usize = 7;
    let m_zero: Matrix<f32> = Matrix::zero(rows, cols);
    let dim: (usize, usize) = m_zero.dim();
    assert_eq!(dim, (rows, cols));

    for i in 0..rows
    {
        for k in 0..cols
        {
            assert_relative_eq!(*(m_zero.get(i, k)), 0.0);
        }
    }
}

#[test]
fn one()
{
    let rows: usize = 5;
    let m_ones: Matrix<f32> = Matrix::one(rows);
    let dim: (usize, usize) = m_ones.dim();
    assert_eq!(dim, (rows, rows));

    for i in 0..rows
    {
        for k in 0..rows
        {
            if i == k
            {
                assert_relative_eq!(*m_ones.get(i, k), 1.0);
            }
            else
            {
                assert_relative_eq!(*m_ones.get(i, k), 0.0);
            }
        }
    }
}

#[test]
fn ones()
{
    let rows: usize = 3;
    let columns: usize = 5;
    let m_ones: Matrix<f64> = Matrix::ones(rows, columns);

    let ones_ref: Matrix<f64> = matrix![    1.0, 1.0, 1.0, 1.0, 1.0;
                                1.0, 1.0, 1.0, 1.0, 1.0;
                                1.0, 1.0, 1.0, 1.0, 1.0];

    assert_eq!(m_ones, ones_ref);
}

#[test]
fn get_column()
{
    let a: Matrix<f32> = matrix![4.0, 1.0, -3.0, 2.0; 1.0, 2.0, 0.0, 1.0; -2.0, 0.0, 3.0, -2.0; 2.0, 1.0, -2.0,
    -1.0];

    let x: Vector<f32> = a.get_column(0);

    let x_ref: Vector<f32> = Vector::new_column(4, vec![4.0, 1.0, -2.0, 2.0]);

    for i in 0..4
    {
        assert_relative_eq!(*x.get(i), *x_ref.get(i));
    }
}

#[test]
fn get_row()
{
    let a: Matrix<f64> = matrix![   4.0, 1.0, -2.0, 2.0;
                                    1.0, 2.0, 3.0, 1.0;
                                    -2.0, 0.0, 3.0, -2.0;
                                    2.0, 3.0, -2.0, -1.0];

    let x: Vector<f64> = a.get_row(1);
    let x_ref: Vector<f64> = Vector::new_row(4, vec![1.0, 2.0, 3.0, 1.0]);

    for i in 0..4
    {
        assert_relative_eq!(*(x.get(i)), *(x_ref.get(i)));
    }
}

#[test]
fn givens()
{
    let m: usize = 4;
    let i: usize = 1;
    let j: usize = 2;
    let theta: f32 = 1.0;
    let c: f32 = theta.cos();
    let s: f32 = theta.sin();

    let givens: Matrix<f32> = Matrix::givens(m, i, j, c, s);

    assert_relative_eq!(*(givens.get(0, 0)), 1.0);
    assert_relative_eq!(*(givens.get(i, i)), c);
    assert_relative_eq!(*(givens.get(j, j)), c);
    assert_relative_eq!(*(givens.get(j, i)), -s);
    assert_relative_eq!(*(givens.get(i, j)), s);
}

#[cfg(feature = "native")]
#[test]
fn givens2()
{
    let (c, s): (f64, f64) = Matrix::givens_cosine_sine_pair(3.0, 5.0);

    assert_relative_eq!(-0.5144957554275266, c);
    assert_relative_eq!(0.8574929257125443, s);
}

#[test]
fn givens_3()
{
    let m: usize = 4;
    let i: usize = 1;
    let j: usize = 2;
    let theta: f64 = 1.0;
    let c: f64 = theta.cos();
    let s: f64 = theta.sin();
    let givens: Matrix<f64> = Matrix::givens(m, i, j, c, s);
    let givens_t: Matrix<f64> = givens.clone().transpose();
    let res_ref: Matrix<f64> = Matrix::one(m);
    let res: Matrix<f64> = givens_t * givens;

    assert_relative_eq!(res_ref, res, epsilon=0.00001, max_relative=1.0e-10);
}

#[test]
fn determinant_0()
{
    let a: Matrix<f64> = matrix![-2.0];
    let det: f64 = a.det();

    assert_abs_diff_eq!(-2.0, det, epsilon=1.0e-10);
}

#[test]
fn determinant_1()
{
    let a: Matrix<f64> = matrix![   1.0, -2.0;
                                    3.0, -7.0];
    let det: f64 = a.det();

    assert_abs_diff_eq!(-1.0, det, epsilon=1.0e-10);
}

#[test]
fn determinant_2()
{
    let a: Matrix<f32> = matrix![   1.0, -2.0, 3.0;
                                    2.0, -5.0, 12.0;
                                    1.0, 2.0, -10.0];
    let det: f32 = a.det();

    assert_abs_diff_eq!(-11.0, det, epsilon=1.0e-4);
}

#[test]
fn determinant_3()
{
    let a: Matrix<f64> = matrix![   4.0, 1.0, -2.0, 2.0;
                                    1.0, 2.0, 0.0, -2.0;
                                    0.0, 3.0, -2.0, 2.0;
                                    2.0, 1.0, -2.0, -1.0];

    let det: f64 = a.det();

    assert_abs_diff_eq!(76.0, det, epsilon=1.0e-10);
}

#[test]
fn determinant_4()
{
    let a: Matrix<f64> = matrix![   -9.0, -8.0, -7.0;
                                    -6.0, 5.0, -6.0;
                                    -7.0, -8.0, -9.0];
    let det: f64 = a.det();

    assert_abs_diff_eq!(352.0, det, epsilon=1.0e-10);
}

#[test]
fn trace_0()
{
    let a: Matrix<f64> = matrix![0.0];
    let tr: f64 = a.trace();

    assert_abs_diff_eq!(0.0, tr, epsilon=1.0e-10);
}

#[test]
fn trace_1()
{
    let a: Matrix<f64> = matrix![-9.0];
    let tr: f64 = a.trace();

    assert_abs_diff_eq!(-9.0, tr, epsilon=1.0e-10);
}

#[test]
fn trace_2()
{
    let a: Matrix<f64> = matrix![   1.0, -2.0;
                                    3.0, -7.0];
    let tr: f64 = a.trace();

    assert_abs_diff_eq!(-6.0, tr, epsilon=1.0e-10);
}

#[test]
fn householder_0()
{
    let v: Vector<f64> = Vector::new_column(3, vec![1.0, 2.0, 3.0]);
    let h: Matrix<f64> = Matrix::householder(&v, 0);

    let h_ref: Matrix<f64> = matrix![   -0.2672612419124243, -0.5345224838248488, -0.8017837257372731;
                                        -0.5345224838248488, 0.7745419205884382, -0.33818711911734267;
                                        -0.8017837257372731, -0.33818711911734267, 0.4927193213239861];

    assert_relative_eq!(h_ref, h, epsilon=0.00001, max_relative=1.0e-10);
}

#[test]
fn householder_1()
{
    let v: Vector<f64> = Vector::new_column(3, vec![1.0, 2.0, 3.0]);
    let h: Matrix<f64> = Matrix::householder(&v, 1);

    let h_ref: Matrix<f64> = matrix![   1.0, 0.0, 0.0;
                                        0.0, -0.5547001962252291, -0.8320502943378437;
                                        0.0, -0.8320502943378437, 0.5547001962252291];

    assert_relative_eq!(h_ref, h, epsilon=0.00001, max_relative=1.0e-10);
}

#[test]
fn householder_2()
{
    let v: Vector<f64> = Vector::new_column(3, vec![1.0, 2.0, 3.0]);
    let h: Matrix<f64> = Matrix::householder(&v, 2);

    let h_ref: Matrix<f64> = matrix![   1.0, 0.0, 0.0;
                                        0.0, 1.0, 0.0;
                                        0.0, 0.0, -1.0];

    assert_relative_eq!(h_ref, h, epsilon=0.00001, max_relative=1.0e-10);
}

#[test]
fn slice_get_0()
{
    let a: Matrix<f32> = matrix![   1.0, 2.0, 3.0;
                                    4.0, 5.0, 6.0;
                                    7.0, 8.0, 9.0];

    let a_ref: Matrix<f32> = matrix![5.0];

    let slice: Matrix<f32> = a.get_slice(1, 1, 1, 1);

    assert_relative_eq!(a_ref, slice);
}

#[test]
fn slice_get_1()
{
    let a: Matrix<f32> = matrix![   1.0, 2.0, 3.0;
                                    4.0, 5.0, 6.0;
                                    7.0, 8.0, 9.0];

    let a_ref: Matrix<f32> = matrix![   5.0, 6.0;
                                        8.0, 9.0];

    let slice: Matrix<f32> = a.get_slice(1, 2, 1, 2);

    assert_relative_eq!(a_ref, slice, epsilon=0.0000001, max_relative=1.0e-10);
}

#[test]
fn slice_set_1()
{
    let a: Matrix<f32> = matrix![   1.0, 2.0, 3.0;
                                    4.0, 5.0, 6.0;
                                    7.0, 8.0, 9.0];

    let a_ref: Matrix<f32> = matrix![   1.0, 2.0, 3.0;
                                        4.0, -5.0, -6.0;
                                        7.0, -8.0, -9.0];

    let b: Matrix<f32> = matrix![   -5.0, -6.0;
                                    -8.0, -9.0];

    let slice: Matrix<f32> = a.set_slice(&b, 1, 1);

    assert_relative_eq!(a_ref, slice);
}

#[test]
fn householder_bidiagonal_0()
{
    let a: Matrix<f32> = matrix![   1.0, 2.0, 3.0;
                                    4.0, 5.0, 6.0;
                                    7.0, 8.0, 9.0;
                                    10.0, 11.0, 12.0];

    let b_ref: Matrix<f32> = matrix![   -12.884099, 21.876434, 0.0;
                                        0.0, 2.2462382, -0.61328155;
                                        0.0, 0.0, -0.000000029802322;
                                        0.0, 0.0, 0.0];

    let v_ref: Matrix<f32> = matrix![   1.0, 0.0, 0.0;
                                        0.0, -0.6670023, -0.7450557;
                                        0.0, -0.7450557, 0.6670023];

    let (_u, b, v): (Matrix<f32>, Matrix<f32>, Matrix<f32>) = a.householder_bidiag();

    assert_relative_eq!(b_ref, b, epsilon=0.00001, max_relative=1.0e-10);
    assert_relative_eq!(v_ref, v, epsilon=0.00001, max_relative=1.0e-10);
}

#[test]
fn householder_bidiagonal_1()
{
    let a: Matrix<f64> = matrix![   1.0, 5.0, 3.0;
                                    1.0, 0.0, -7.0;
                                    3.0, 8.0, 9.0];

    let b_ref: Matrix<f64> = matrix![   -3.3166247903554, 11.15999348321739, 0.0;
                                        0.0, -8.27496123318713, 5.336122204714563;
                                        0.0, 0.0, 2.5505610873193763];

    let v_ref: Matrix<f64> = matrix![   1.0, 0.0, 0.0;
                                        0.0, -0.783497679089534, -0.6213947110020441;
                                        0.0, -0.6213947110020441, 0.7834976790895338];

    let (_u, b, v): (Matrix<f64>, Matrix<f64>, Matrix<f64>) = a.householder_bidiag();

    assert_relative_eq!(b_ref, b, epsilon=0.00001, max_relative=1.0e-10);
    assert_relative_eq!(v_ref, v, epsilon=0.00001, max_relative=1.0e-10);
}

#[test]
fn householder_bidiagonal_2()
{
    let a: Matrix<f32> = matrix![   4.0, 1.0, -2.0, 2.0;
                                    1.0, 2.0, 0.0, -2.0;
                                    0.0, 3.0, -2.0, 2.0;
                                    2.0, 1.0, -2.0, -1.0];

    let b_ref: Matrix<f32> = matrix![   -4.582576, 3.2659864, 0.0, 0.0;
                                        0.0, -3.7764935, -1.5535977, 0.0;
                                        0.0, 0.0, 1.4568509, -1.203649;
                                        0.0, 0.0, 0.0, -3.014395];

    let (_u, b, _v): (Matrix<f32>, Matrix<f32>, Matrix<f32>) = a.householder_bidiag();

    assert_relative_eq!(b_ref, b, epsilon=0.00001, max_relative=1.0e-10);
}

#[test]
fn rot_0()
{
    let f: f64 = 0.0;
    let g: f64 = -3.0;
    let (c_ref, s_ref, r_ref): (f64, f64, f64) = (0.0, 1.0, -3.0);

    let (c, s, r): (f64, f64, f64) = Matrix::rot(f, g);

    assert_relative_eq!(c_ref, c);
    assert_relative_eq!(s_ref, s);
    assert_relative_eq!(r_ref, r);
}

#[test]
fn rot_1()
{
    let f: f64 = 2.0;
    let g: f64 = 3.0;
    let (c_ref, s_ref, r_ref): (f64, f64, f64) =
        (0.554700196225229, 0.8320502943378437, 3.6055512754639896);

    let (c, s, r): (f64, f64, f64) = Matrix::rot(f, g);

    assert_relative_eq!(c_ref, c);
    assert_relative_eq!(s_ref, s);
    assert_relative_eq!(r_ref, r);
}

#[test]
fn rot_2()
{
    let f: f64 = 3.0;
    let g: f64 = -5.0;
    let (c_ref, s_ref, r_ref): (f64, f64, f64) =
        (-0.5144957554275266, 0.8574929257125443, -5.8309518948453);

    let (c, s, r): (f64, f64, f64) = Matrix::rot(f, g);

    assert_relative_eq!(c_ref, c, epsilon=0.0000001, max_relative=1.0e-10);
    assert_relative_eq!(s_ref, s, epsilon=0.0000001, max_relative=1.0e-10);
    assert_relative_eq!(r_ref, r, epsilon=0.0000001, max_relative=1.0e-10);
}

#[test]
fn svd_0()
{
    let a: Matrix<f64> = matrix![   4.0, 1.0, -2.0, 2.0;
                                    1.0, 2.0, 0.0, -2.0;
                                    0.0, 3.0, -2.0, 2.0;
                                    2.0, 1.0, -2.0, -1.0];

    let (u, s, v): (Matrix<f64>, Matrix<f64>, Matrix<f64>) = a.dec_sv();

    let u_ref: Matrix<f64> = matrix![   0.750618975420566, -0.10396419803551285, -0.5550030946996072, -0.3431241235169467;
                                        0.16660611095451502, 0.7151480228364393, 0.42071909480760994, -0.5327299057590199;
                                        0.49145768126990524, -0.4896500732330627, 0.7174571018374751, 0.06298777918484924;
                                        0.40900232890779536, 0.48784993260512777, -0.014910236788856494, 0.7710364602559809];

    let v_ref: Matrix<f64> = matrix![   0.641413939431974, 0.3770763401833474, -0.5915891468488664, -0.3105219369173845;
                                        0.4773429354587908, 0.10210161763764197, 0.7839592576263352, -0.3835711981175503;
                                        -0.5312275172316614, 0.06255922864595094, -0.09543998328421795, -0.8395087119487312;
                                        0.2802306286004916, -0.9184089700459931, -0.1622386307366823, -0.2273200062245417];

    let s_ref: Matrix<f64> = matrix![   6.216089837372844, 0.0, 0.0, 0.0;
                                        0.0, 3.3812545631600996, 0.0, 0.0;
                                        0.0, 0.0, 3.0918649662716553, 0.0;
                                        0.0, 0.0, 0.0, 1.1694937978293738];

    assert_relative_eq!(u_ref, u, epsilon=0.0000001, max_relative=1.0e-10);
    assert_relative_eq!(s_ref, s, epsilon=0.0000001, max_relative=1.0e-10);
    assert_relative_eq!(v_ref, v, epsilon=0.0000001, max_relative=1.0e-10);
}

#[test]
fn svd_1()
{
    let a: Matrix<f64> = matrix![   4.0, 1.0, -2.0, 2.0;
                                    1.0, 2.0, 0.0, -2.0;
                                    0.0, 3.0, -2.0, 2.0;
                                    2.0, 1.0, -2.0, -1.0];

    let (u, _s, v): (Matrix<f64>, Matrix<f64>, Matrix<f64>) = a.dec_sv();

    let u_ref: Matrix<f64> = matrix![   0.750618975420566, -0.10396419803551285, -0.5550030946996072, -0.3431241235169467;
                                        0.16660611095451502, 0.7151480228364393, 0.42071909480760994, -0.5327299057590199;
                                        0.49145768126990524, -0.4896500732330627, 0.7174571018374751, 0.06298777918484924;
                                        0.40900232890779536, 0.48784993260512777, -0.014910236788856494, 0.7710364602559809];

    let v_ref: Matrix<f64> = matrix![   0.641413939431974, 0.3770763401833474, -0.5915891468488664, -0.3105219369173845;
                                        0.4773429354587908, 0.10210161763764197, 0.7839592576263352, -0.3835711981175503;
                                        -0.5312275172316614, 0.06255922864595094, -0.09543998328421795, -0.8395087119487312;
                                         0.2802306286004916, -0.9184089700459931, -0.1622386307366823, -0.2273200062245417];

    assert_relative_eq!(u_ref, u, epsilon=0.0000001, max_relative=1.0e-10);
    assert_relative_eq!(v_ref, v, epsilon=0.0000001, max_relative=1.0e-10);
}

#[test]
fn apply_0()
{
    let a: Matrix<f64> = matrix![   1.0, -3.0, 3.0;
                                    3.0, -5.0, 3.0;
                                    6.0, -6.0, 4.0];

    let a_ref: Matrix<f64> = matrix![   -1.0, 3.0, -3.0;
                                        -3.0, 5.0, -3.0;
                                        -6.0, 6.0, -4.0];

    let b: Matrix<f64> = a.apply(&|x| -x);

    assert_relative_eq!(a_ref, b, epsilon=1.0e-10);
}

#[test]
fn subst_backward()
{
    let a: Matrix<f64> = matrix![   1.0, 2.0, 3.0;
                                    0.0, 4.0, 5.0;
                                    0.0, 0.0, 6.0];
    let b: Vector<f64> = vector![7.0; 8.0; 9.0];

    let c_ref: Vector<f64> = vector![2.25; 0.125; 1.5];

    let c: Vector<f64> = a.substitute_backward(b);

    assert_relative_eq!(c, c_ref, epsilon=1.0e-10);
}

#[test]
fn subst_forward()
{
    let a: Matrix<f64> = matrix![   6.0, 0.0, 0.0;
                                    5.0, 4.0, 0.0;
                                    3.0, 2.0, 1.0];

    let b: Vector<f64> = vector![9.0; 8.0; 7.0];

    let c_ref: Vector<f64> = vector![1.5; 0.125; 2.25];

    let c: Vector<f64> = a.substitute_forward(b);

    assert_relative_eq!(c, c_ref, epsilon=1.0e-10);
}

#[test]
fn solve_0()
{
    let a: Matrix<f64> = matrix![6.0, 2.0, -1.0; -3.0, 5.0, 3.0; -2.0, 1.0, 3.0];
    let b: Vector<f64> = vector![48.0; 49.0; 24.0];

    let x: Vector<f64> = a.solve(&b).unwrap();
    let x_ref: Vector<f64> = vector![7.0; 8.0; 10.0];

    assert_relative_eq!(x, x_ref, epsilon=10e-10);
}

#[test]
fn pinv_0()
{
    let a: Matrix<f64> = matrix![1.0, 2.0, 3.0; 4.0, 5.0, 6.0; 8.0, 8.0, 9.0];
    let a_pinv: Matrix<f64> = a.pinv();
    let a_pinv_ref: Matrix<f64> = matrix![  1.0, -2.0, 1.0;
                                            -4.0, 5.0, -2.0;
                                            2.6666666666666474, -2.66666666666661, 1.0];

    assert_relative_eq!(a_pinv, a_pinv_ref, epsilon=10e-10);
}
