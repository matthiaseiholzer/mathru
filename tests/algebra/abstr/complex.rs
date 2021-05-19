use mathru::algebra::abstr::{One, Zero};
use mathru::algebra::abstr::{Complex};
use mathru::elementary::{Trigonometry, Exponential, Hyperbolic, Power};
use mathru::algebra::abstr::cast::{FromPrimitive};
use std::f64::consts;


#[test]
fn arg0()
{
    let a: Complex<f64> = Complex::new(-3.5, -6.0);
    let angle : Complex<f64> = Complex::new(-6.0_f64.atan2(-3.5), 0.0_f64);
    assert_eq!(angle, a.arg());
}

#[test]
fn arg1()
{
    let phi: Complex<f64> = Complex::new(0.0, 0.0_f64);
    let a: Complex<f64>= Complex::new(0.5, 0.0);
    assert_eq!(phi, a.arg());
}

#[test]
fn add()
{
    let a : Complex<f32> = Complex::new(-3.5, 6.5);
    let b : Complex<f32> = Complex::new(2.0, -3.0);
    let sum: Complex<f32> = Complex::new(-1.5, 3.5);
    assert_eq!(sum, a + b);
}

#[test]
fn mul()
{
    let a : Complex<f32> = Complex::new(-3.5, 6.5);
    let b : Complex<f32> = Complex::new(2.0, -3.0);
    let prod: Complex<f32> = Complex::new(-3.5 * 2.0 - 6.5 * -3.0, 6.5 * 2.0 + -3.5 * -3.0);
    assert_eq!(prod, a * b);
}

#[test]
fn mulassign()
{
    let mut a : Complex<f32> = Complex::new(-3.5, 6.5);
    let b : Complex<f32> = Complex::new(2.0, -3.0);
    let prod : Complex<f32> = a * b;
    a *= b;
    assert_eq!(prod, a);
}

#[test]
fn zero()
{
    let reference : Complex<f32> = Complex::new(0.0, 0.0);
    assert_eq!(reference, Complex::zero());
}

#[test]
fn one()
{
    let reference : Complex<f64> = Complex::new(1.0, 0.0);
    assert_eq!(reference, Complex::one());
}

#[test]
fn conj()
{
    let a_real: f32 = -3.5;
    let a_imag: f32 = 6.5;
    let a : Complex<f32> = Complex::new(a_real, a_imag);
    let conj : Complex<f32> = Complex::new(a_real, -a_imag);
    assert_eq!(conj, a.conj());
}

#[test]
fn div0()
{
    let a_real: f32 = -3.5;
    let a_imag: f32 = 6.5;
    let b_real: f32 = 2.0;
    let b_imag: f32 = -3.0;
    let a : Complex<f32> = Complex::new(a_real, a_imag);
    let b : Complex<f32> = Complex::new(b_real, b_imag);
    let prod: Complex<f32> = Complex::new((a_real * b_real + a_imag * b_imag)/(b_real* b_real +
        b_imag * b_imag), (a_imag * b_real - a_real * b_imag)/(b_real * b_real + b_imag * b_imag));

    assert_eq!(prod, a / b);
}

#[test]
fn div1()
{
    let a_real: f32 = -3.5;
    let a_imag: f32 = 6.5;
    let b_real: f32 = 0.0;
    let b_imag: f32 = 3.0;
    let a : Complex<f32> = Complex::new(a_real, a_imag);
    let b : Complex<f32> = Complex::new(b_real, b_imag);
    let prod: Complex<f32> = Complex::new((a_real * b_real + a_imag * b_imag)/(b_real* b_real +
        b_imag * b_imag), (a_imag * b_real - a_real * b_imag)/(b_real * b_real + b_imag * b_imag));

    assert_eq!(prod, a / b);
}

#[test]
fn sub()
{
    let a : Complex<f32> = Complex::new(-3.5, 6.5);
    let b : Complex<f32> = Complex::new(2.0, -3.0);
    let diff: Complex<f32> = Complex::new(-5.5, 9.5);
    assert_eq!(diff, a - b);
}

#[test]
fn neg0()
{
    let a_real : f32 = 1.0;
    let a_imag : f32 = 2.0;
    let uut : Complex<f32> = Complex::new(a_real, a_imag);
    assert_eq!(Complex::new(-a_real, -a_imag), -uut);
}

#[test]
fn exp()
{
    let z: Complex<f64> = Complex::new(1.0, 2.0);
    let uut: Complex<f64> = z.exp();

    let refer: Complex<f64> = Complex::new(1.0.exp() * 2.0.cos(), 1.0.exp() * 2.0.sin());

    assert_eq!(refer, uut);
}

#[test]
fn sin()
{
    let a: f64 = 1.0;
    let b: f64 = 2.0;
    let z: Complex<f64> = Complex::new(a, b);
    let real: f64 = (-(-b).exp() * a.sin() - b.exp() * a.sin()) / -2.0;
    let imag: f64 = ((-b).exp() * a.cos() - b.exp() * a.cos()) / -2.0;

    let uut: Complex<f64> = z.sin();

    let refer: Complex<f64> = Complex::new(real, imag);
    assert_eq!(refer, uut);
}

#[test]
fn cos()
{
    let a: f64 = 1.0;
    let b: f64 = 2.0;
    let z: Complex<f64> = Complex::new(a, b);
    let real: f64 = ((-b).exp() * a.cos() + b.exp() * (a).cos()) / 2.0;
    let imag: f64 = ((-b).exp() * a.sin() - b.exp() * (a).sin()) / 2.0;

    let uut: Complex<f64> = z.cos();

    let refer: Complex<f64> = Complex::new(real, imag);

    assert_eq!(refer, uut);
}

#[test]
fn tan()
{
    let a: f64 = 1.0;
    let b: f64 = 2.0;
    let z: Complex<f64> = Complex::new(a, b);
    let refer: Complex<f64> = z.sin() / z.cos();

    let uut: Complex<f64> = z.tan();

    assert_eq!(refer, uut);
}

 #[test]
fn cot()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let refer: Complex<f64> = Complex::one() / a.tan();

    assert_eq!(refer, a.cot());
}

#[test]
fn sec()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let refer: Complex<f64> = Complex::one() / a.cos();

    assert_eq!(refer, a.sec());
}

#[test]
fn csc()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let refer: Complex<f64> = Complex::one() / a.sin();

    assert_eq!(refer, a.csc());
}

#[test]
fn ln()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let refer: Complex<f64> = Complex::new(5.0_f64.powf(0.5_f64).ln(), 2.0_f64.arctan2(1.0_f64));

    assert_eq!(refer, a.ln());
}

#[test]
fn arctan0()
{
    let a: Complex<f64> = Complex::new(3.0 * std::f64::consts::PI / 2.0, 0.0);
    let refer: Complex<f64> = Complex::new(1.3616916829711636, 0.0);

    assert_eq!(refer, a.arctan());
}

#[test]
fn arctan1()
{
    let a: Complex<f64> = Complex::new(0.0_f64, 2.0_f64);
    let refer: Complex<f64> = Complex::new(std::f64::consts::PI / 2.0, (4.0_f64 / 5.0_f64).artanh() / 2.0_f64);

    assert_eq!(refer, a.arctan());
}

#[test]
fn arctan2()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let refer: Complex<f64> = Complex::new((2.0_f64.atan() + std::f64::consts::PI / 2.0_f64) / 2.0_f64, (4.0_f64 / 6.0_f64).atanh() / 2.0_f64);

    assert_eq!(refer, a.arctan());
}

#[test]
fn arccot()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let refer: Complex<f64> = (Complex::one() / a).arctan();

    assert_eq!(refer, a.arccot());
}

#[test]
fn pow0()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let b: Complex<f64> = Complex::new(-2.0_f64, -1.0_f64);
    let refer: Complex<f64> = Complex::new(-0.6006033457684014, -0.07399065302898929);

    assert_eq!(refer, a.pow(b));
}

#[test]
fn pow1()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let b: Complex<f64> = Complex::new(-2.0_f64, 1.0_f64);
    let refer: Complex<f64> = Complex::new(0.010610396101816041, -0.06524284357147048);

    assert_eq!(refer, a.pow(b));
}

#[test]
fn pow2()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let b: Complex<f64> = Complex::new(-0.5_f64, -0.8_f64);
    let refer: Complex<f64> = Complex::new(0.591571342038212, -1.5097505998220102);

    assert_eq!(refer, a.pow(b));
}

#[test]
fn arcsin()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let refer: Complex<f64> = Complex::new(0.4270785863924768, 1.5285709194809995);

    assert_eq!(refer, a.arcsin());
}

#[test]
fn arccos()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let refer: Complex<f64> = Complex::new(std::f64::consts::PI / 2.0_f64, 0.0_f64) - a.arcsin();

    assert_eq!(refer, a.arccos());
}

#[test]
fn arcsec()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let refer: Complex<f64> = (Complex::one() / a).arccos();

    assert_eq!(refer, a.arcsec());
}

#[test]
fn arccsc()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let refer: Complex<f64> = (Complex::one() / a).arcsin();

    assert_eq!(refer, a.arccsc());
}

#[test]
fn sinh()
{
    let a: f64 = 1.0;
    let b: f64 = 2.0;
    let z: Complex<f64> = Complex::new(a, b);

    let uut: Complex<f64> = z.sinh();

    let refer: Complex<f64> = Complex::new(0.0_f64, -1.0_f64) * Complex::new(-2.0_f64, 1.0_f64).sin();
    assert_eq!(refer, uut);
}

#[test]
fn cosh()
{
    let a: f64 = 1.0;
    let b: f64 = 2.0;
    let z: Complex<f64> = Complex::new(a, b);

    let uut: Complex<f64> = z.cosh();

    let refer: Complex<f64> = Complex::new(-2.0_f64, 1.0_f64).cos();

    assert_eq!(refer, uut);
}

#[test]
fn tanh()
{
    let a: f64 = 1.0;
    let b: f64 = 2.0;
    let z: Complex<f64> = Complex::new(a, b);
    let refer: Complex<f64> = z.sinh() / z.cosh();

    let uut: Complex<f64> = z.tanh();
    assert_eq!(refer, uut);
}


 #[test]
fn coth()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let refer: Complex<f64> = a.cosh() / a.sinh();

    assert_eq!(refer, a.coth());
}

#[test]
fn sech()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let refer: Complex<f64> = Complex::new(-2.0_f64, 1.0_f64).sec();

    assert_eq!(refer, a.sech());
}

#[test]
fn csch()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let refer: Complex<f64> = Complex::new(0.0_f64, -1.0_f64) * Complex::new(-2.0_f64, 1.0_f64).csc();

    assert_eq!(refer, a.csch());
}

#[test]
fn artanh()
{
    let a: Complex<f64> = Complex::new(0.5_f64, -0.4_f64);
    let f: Complex<f64> = Complex::new(0.5_f64, 0.0_f64);
    let refer: Complex<f64> = ((Complex::one() + a) / (Complex::one() - a)).ln() * f;

    assert_eq!(refer, a.artanh());
}



#[test]
fn arcoth()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let f: Complex<f64> = Complex::new(0.5_f64, 0.0_f64);
    let refer: Complex<f64> = ((a + Complex::one()) / (a - Complex::one())).ln() * f;

    assert_eq!(refer, a.arcoth());
}

#[test]
fn arsinh()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let pow: Complex<f64> = Complex::new(0.5_f64, 0.0_f64);

    let refer: Complex<f64> = (a + (a * a + Complex::one()).pow(pow)).ln();

    assert_eq!(refer, a.arsinh());
}

#[test]
fn arcosh()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let pow: Complex<f64> = Complex::new(0.5_f64, 0.0_f64);

    let refer: Complex<f64> = (a + (a * a - Complex::one()).pow(pow)).ln();

    assert_eq!(refer, a.arcosh());
}

#[test]
fn arcsech()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let refer: Complex<f64> = (Complex::one() / a).arcosh();

    assert_eq!(refer, a.arsech());
}

#[test]
fn arccsch()
{
    let a: Complex<f64> = Complex::new(1.0_f64, 2.0_f64);
    let refer: Complex<f64> = (Complex::one() / a).arsinh();

    assert_eq!(refer, a.arcsch());
}

#[test]
fn from_f64()
{
    let a: Complex<f64> = Complex::from_f64(consts::PI);
    let refer: Complex<f64> = Complex::pi();

    assert_eq!(refer, a);
}