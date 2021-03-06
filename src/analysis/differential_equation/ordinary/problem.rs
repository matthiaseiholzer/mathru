use crate::{
    algebra::{
        abstr::Real,
        linear::{Matrix, Vector},
    },
    analysis::differential_equation::ordinary::{ExplicitODE, ImplicitODE},
};

///
/// ```math
/// y_{1}^{'}(x) = (I_{2} - I_{3})/I_{1} * y_{2}(x) * y_{3}(x)
/// y_{2}^{'}(x) = (I_{3} - I_{1})/I_{2} * y_{3}(x) * y_{1}(x)
/// y_{3}^{'}(x) = (I_{1} - I_{2})/I_{3} * y_{1}(x) * y_{2}(x) + f(x)
///
/// f =
/// ```
pub struct Euler<T>
{
    i1: T,
    i2: T,
    i3: T,

    time_span: (T, T),
    init_cond: Vector<T>,
}

impl<T> Default for Euler<T> where T: Real
{
    fn default() -> Euler<T>
    {
        Euler { i1: T::from_f64(0.5),
                i2: T::from_f64(2.0),
                i3: T::from_f64(3.0),
                time_span: (T::from_f64(0.0), T::from_f64(20.0)),
                init_cond: vector![T::from_f64(1.0); T::from_f64(0.0); T::from_f64(0.9)] }
    }
}

impl<T> ExplicitODE<T> for Euler<T> where T: Real
{
    fn func(self: &Self, x: &T, y: &Vector<T>) -> Vector<T>
    {
        let y_1s: T = ((self.i2 - self.i3) / self.i1) * (*y.get(1) * *y.get(2));
        let y_2s: T = ((self.i3 - self.i1) / self.i2) * (*y.get(2) * *y.get(0));

        let f: T;
        if *x >= T::from_f64(3.0) * T::pi() && *x <= T::from_f64(4.0) * T::pi()
        {
            f = T::from_f64(0.25) * x.sin() * x.sin();
        }
        else
        {
            f = T::zero();
        }
        let y_3s: T = ((self.i1 - self.i2) / self.i3) * (*y.get(0) * *y.get(1)) + f;
        return vector![y_1s; y_2s; y_3s];
    }

    fn time_span(self: &Self) -> (T, T)
    {
        return self.time_span;
    }

    fn init_cond(self: &Self) -> Vector<T>
    {
        return self.init_cond.clone();
    }
}

impl<T> ImplicitODE<T> for Euler<T> where T: Real
{
    fn func(self: &Self, x: T, y: &Vector<T>) -> Vector<T>
    {
        let a: T = (self.i2 - self.i3) / self.i1;
        let b: T = (self.i3 - self.i1) / self.i2;
        let c: T = (self.i1 - self.i2) / self.i3;

        let y_1s: T = a * (*y.get(1) * *y.get(2));
        let y_2s: T = b * (*y.get(2) * *y.get(0));

        let f: T;
        if x >= T::from_f64(3.0) * T::pi() && x <= T::from_f64(4.0) * T::pi()
        {
            f = T::from_f64(0.25) * x.sin() * x.sin();
        }
        else
        {
            f = T::zero();
        }
        let y_3s: T = c * (*y.get(0) * *y.get(1)) + f;
        return vector![y_1s; y_2s; y_3s];
    }

    fn jacobian(self: &Self, _x: T, y: &Vector<T>) -> Matrix<T>
    {
        let a: T = (self.i2 - self.i3) / self.i1;
        let b: T = (self.i3 - self.i1) / self.i2;
        let c: T = (self.i1 - self.i2) / self.i3;

        return matrix![ T::zero(), a * *y.get(2), a * *y.get(1);
                        b * *y.get(2), T::zero(), b * *y.get(0);
                        c * *y.get(1), c * *y.get(0), T::zero()];
    }

    fn time_span(self: &Self) -> (T, T)
    {
        return self.time_span;
    }

    fn init_cond(self: &Self) -> Vector<T>
    {
        return self.init_cond.clone();
    }
}

/// Van der Pol oscillator
/// ```math
/// x_{1}^{'}(t) = x_{2}(t) \\
/// x_{2}^{'}(t) = \epsilon((1 - x_{1}(t)^{2})x_{2}(t) - x_{1}(t)) \\
/// ```
///
/// ```math
/// x_{1}(0) = 1 \\
/// x_{2}(0) = 0 \\
/// \epsilon = 0.1 \\
/// ```
pub struct VanDerPolOsc<T>
{
    epsilon: T,

    time_span: (T, T),
    init_cond: Vector<T>,
}

impl<T> Default for VanDerPolOsc<T> where T: Real
{
    fn default() -> VanDerPolOsc<T>
    {
        VanDerPolOsc { epsilon: T::from_f64(0.1),
                       time_span: (T::from_f64(0.0), T::from_f64(30.0)),
                       init_cond: vector![T::from_f64(1.0); T::from_f64(0.0)] }
    }
}
/// Implicit ordinary differential equation
///
/// $`0 = f(t, x(t), x^{'}(t), \dots, x^{n}(t))`$
impl<T> ImplicitODE<T> for VanDerPolOsc<T> where T: Real
{
    fn func(self: &Self, _t: T, x: &Vector<T>) -> Vector<T>
    {
        let x_1 = *x.get(0);
        let x_2 = *x.get(1);
        return vector![x_2; self.epsilon * (T::one() - (x_1 * x_1)) * x_2 - x_1];
    }

    fn jacobian(self: &Self, _t: T, x: &Vector<T>) -> Matrix<T>
    {
        let x_1 = *x.get(0);
        let x_2 = *x.get(1);
        return matrix![T::zero(), T::one(); -T::from_f64(2.0) * self.epsilon * x_1 * x_2  - T::one(), (T::one() - x_1 *
		x_1) * self.epsilon];
    }

    fn time_span(self: &Self) -> (T, T)
    {
        return self.time_span;
    }

    fn init_cond(self: &Self) -> Vector<T>
    {
        return self.init_cond.clone();
    }
}
