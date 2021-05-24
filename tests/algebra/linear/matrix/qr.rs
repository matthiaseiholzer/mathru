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
                                        0.0, 5.0.pow(0.5);
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
