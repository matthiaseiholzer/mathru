use mathru::{algebra::linear::Matrix};

#[cfg(feature = "native")]
#[test]
fn decompose_qr0()
{
    let a: Matrix<f64> = matrix![   6.0, 5.0, 0.0;
                                    5.0, 1.0, 4.0;
                                    0.0, 4.0, 3.0];

    let (q, r): (Matrix<f64>, Matrix<f64>) = a.dec_qr().qr();

    let q_ref: Matrix<f64> = matrix![   7.682212795973757e-01, 3.326541793600714e-01, 5.469709887444194e-01;
                                        6.401843996644797e-01, -3.991850152320858e-01, -6.563651864933034e-01;
                                        0.0, 8.543959975142890e-01, -5.196224393071984e-01];

    let r_ref: Matrix<f64> = matrix![   7.810249675906654, 4.48129079765136, 2.5607375986579197;
                                        0.0, 4.681669871625427, 0.9664479316145234;
                                        0.0, 0.0, -4.184328063894809];

    assert_relative_eq!(q, q_ref, epsilon=1.0e-10, max_relative=1.0e-10);
    assert_relative_eq!(r, r_ref, epsilon=1.0e-10, max_relative=1.0e-10);
}

#[cfg(feature = "lapack")]
#[test]
fn decompose_qr0()
{
    let a: Matrix<f64> = matrix![   6.0, 5.0, 0.0;
                                    5.0, 1.0, 4.0;
                                    0.0, 4.0, 3.0];

    let (q, r): (Matrix<f64>, Matrix<f64>) = a.dec_qr().qr();

    let q_ref: Matrix<f64> = matrix![   -7.682212795973757e-01, 3.326541793600714e-01, -5.469709887444194e-01;
                                        -6.401843996644797e-01, -3.991850152320858e-01, 6.563651864933034e-01;
                                        -0.000000000000000e+00, 8.543959975142890e-01, 5.196224393071984e-01];

    let r_ref: Matrix<f64> = matrix![  -7.810249675906654, -4.48129079765136, -2.5607375986579197;
                                        0.0, 4.681669871625427, 0.9664479316145234;
                                        0.0, 0.0, 4.184328063894809];

    assert_relative_eq!(q, q_ref, epsilon=1.0e-10, max_relative=1.0e-10);
    assert_relative_eq!(r, r_ref, epsilon=1.0e-10, max_relative=1.0e-10);
}

#[cfg(feature = "native")]
#[test]
fn decompose_qr1()
{
    let a: Matrix<f64> = matrix![   3.0, 5.0;
                                    0.0, 2.0;
                                    0.0, 0.0;
                                    4.0, 5.0];

    let (q, r): (Matrix<f64>, Matrix<f64>) = a.dec_qr().qr();

    let r_ref: Matrix<f64> = matrix![   5.0, 7.0;
                                        0.0, 2.2360679775;
                                        0.0, 0.0;
                                        0.0, 0.0];

    let q_ref: Matrix<f64> = matrix![   0.6, 0.35777087639996635, 0.0, -0.7155417527999327;
                                        0.0, 0.8944271909999159, 0.0, 0.4472135954999579;
                                        0.0, 0.0, 1.0, 0.0;
                                        0.8, -0.2683281572999747, 0.0, 0.5366563145999494 ];

    assert_relative_eq!(r, r_ref, epsilon=1.0e-10);
    assert_relative_eq!(q, q_ref, epsilon=1.0e-10);
    assert_relative_eq!(a, q * r, epsilon=1.0e-10);
}

#[cfg(feature = "lapack")]
#[test]
fn decompose_qr1()
{
    let a: Matrix<f64> = matrix![  3.0, 5.0;
                                    0.0, 2.0;
                                    0.0, 0.0;
                                    4.0, 5.0];

    let r: Matrix<f64> = a.dec_qr().r();

    let r_ref: Matrix<f64> = matrix![  -5.0, -7.0;
                                        0.0, -2.2360679775;
                                        0.0, 0.0;
                                        0.0, 0.0];

    let _q_ref: Matrix<f64> = matrix![  0.6, 0.35777087639996635, 0.0, -0.7155417527999327;
                                        0.0, 0.8944271909999159, 0.0, 0.4472135954999579;
                                        0.0, 0.0, 1.0, 0.0;
                                        0.8, -0.2683281572999747, 0.0, 0.5366563145999494 ];
    assert_relative_eq!(r, r_ref, epsilon=1.0e-10);
}

#[cfg(feature = "native")]
#[test]
fn decompose_qr2()
{
    let a: Matrix<f64> = matrix![   12.0, -51.0, 4.0;
                                    6.0, 167.0, -68.0;
                                    -4.0, 24.0, -41.0];

    let (q, r): (Matrix<f64>, Matrix<f64>) = a.dec_qr().qr();

    let r_ref: Matrix<f64> = matrix![   14.0, 21.0, -14.0;
                                        0.0, 175.0, -70.0;
                                        0.0, 0.0, -35.0];

    let q_ref: Matrix<f64> = matrix![   8.571428571428572e-01, -3.942857142857143e-01, 3.314285714285715e-01;
                                        4.285714285714286e-01, 9.028571428571428e-01, -3.428571428571425e-02;
                                        -2.857142857142858e-01, 1.714285714285714e-01,  9.428571428571428e-01];

    assert_relative_eq!(q, q_ref, epsilon=1.0e-10);
    assert_relative_eq!(r, r_ref, epsilon=1.0e-10);
    assert_relative_eq!(a, &q * &r, epsilon=1.0e-10);
}

#[cfg(feature = "lapack")]
#[test]
fn decompose_qr2()
{
    let a: Matrix<f64> = matrix![  12.0, -51.0, 4.0;
                                    6.0, 167.0, -68.0;
                                    -4.0, 24.0, -41.0];

    let (q, r): (Matrix<f64>, Matrix<f64>) = a.dec_qr().qr();

    let q_ref: Matrix<f64> = matrix![   -8.571428571428572e-01, 3.942857142857143e-01, 3.314285714285715e-01;
                                        -4.285714285714286e-01, -9.028571428571428e-01, -3.428571428571425e-02;
                                        2.857142857142858e-01, -1.714285714285714e-01,  9.428571428571428e-01];

    let r_ref = matrix![    -14.0, -21.0, 14.0;
                            0.0, -175.0, 70.0;
                            0.0, 0.0, -35.0];

    assert_relative_eq!(q, q_ref, epsilon=1.0e-10);
    assert_relative_eq!(r, r_ref, epsilon=1.0e-10);
    assert_relative_eq!(a, &q * &r, epsilon=1.0e-10);
}

// #[test]
// fn decompose_complex_f32()
// {
//     let a: Matrix<Complex<f32>> = matrix![  Complex::new(1.0, 1.0), Complex::new(2.0, 2.0);
//                                             Complex::new(3.0, 3.0), Complex::new(-4.0, 4.0)];
//
//
//     let q_ref: Matrix<Complex<f32>> = matrix![  Complex::new(-0.2236, -0.2236), Complex::new(0.9303, -0.1861);
//                                                 Complex::new(-0.6708, -0.6708), Complex::new(-0.3101, 0.0620)];
//
//     let r_ref: Matrix<Complex<f32>>  = matrix![ Complex::new(-4.4721, 0.0), Complex::new(5.3666 , 0.8944);
//                                                 Complex::zero(), Complex::new(3.2249, 0.0)];
//
//     let (q, r): (Matrix<Complex<f32>>, Matrix<Complex<f32>>) = a.dec_qr().qr();
//
//     assert_relative_eq!(q, q_ref);
//     // assert_relative_eq!(r, r_ref);
//     // assert_relative_eq!(p, p_ref);
//     //
//     assert_relative_eq!(a, &q * &r, epsilon=Complex::new(1.0e-6, 1.0e-6));
// }

// #[test]
// fn decompose_complex_f64()
// {
//     let a: Matrix<Complex<f64>> = matrix![  Complex::new(1.0, 1.0), Complex::new(2.0, 2.0);
//                                             Complex::new(3.0, 3.0), Complex::new(-4.0, 4.0)];
//
//     let l_ref: Matrix<Complex<f64>> = matrix![  Complex::new(1.0, 0.0), Complex::zero();
//                                                 Complex::new(1.0/3.0, 0.0), Complex::new(1.0000, 0.0)];
//
//     let u_ref: Matrix<Complex<f64>>  = matrix![ Complex::new(3.0, 3.0), Complex::new(-4.0, 4.0);
//                                                 Complex::zero(), Complex::new(10.0 / 3.0, 2.0 / 3.0)];
//
//     let p_ref: Matrix<Complex<f64>> = matrix![  Complex::zero(), Complex::new(1.0, 0.0);
//                                                 Complex::new(1.0, 0.0), Complex::zero()];
//
//     let (l, u, p): (Matrix<Complex<f64>>, Matrix<Complex<f64>>, Matrix<Complex<f64>>) = a.dec_lu().unwrap().lup();
//
//     assert_relative_eq!(l, l_ref);
//     assert_relative_eq!(u, u_ref);
//     assert_relative_eq!(p, p_ref);
//
//     assert_relative_eq!(p * l * u, a);
// }