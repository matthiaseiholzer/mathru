use crate::algebra::abstr::Real;
use crate::algebra::linear::{Vector, Matrix};
use crate::algebra::linear::matrix::Solve;
use crate::optimization::{OptimResult, Jacobian};

/// Levenberg-Marquardt method
///
/// input: $` f \colon \mathbb{R}^{n} \to \mathbb{R} `$ with initial approximation $` x_{0} \in \mathbb{R}^{n} `$
///
/// output: $` x_{k} `$
///
/// 1. Initialization: 0 \leq \rho^{-} < \rho^{+} \leq 1, set $` k := 0 `$
/// 2. Solve the equation
/// 3. $`\rho = \frac{\lvert \lvert f(x_{k}) \rvert \rvert_{2}^{2} - \lvert \lvert f(x_{k+1}) \rvert \rvert_{2}^{2}}{2d_{k}^{T}(f^{'}(x_{k})
/// )^T f(x_{k})}`$
/// 4. if $`\rho_{k} > \rho^{-} `$ than $` x_{k + 1} := x_{k} + d_{k}`$, else $`\Delta_{k + 1} := \Delta_{k}/2`$ and go to 2.
/// 5. if $`\rho_{k} > \rho^{+} `$ than $` \Delta_{k + 1} := 2\Delta_{k}`$, else $`\Delta_{k + 1} := \Delta_{k}`$
/// 6. set $` k:= k + 1 `$, go to 2.
pub struct LevenbergMarquardt<T>
{
    iters: u64,
    beta_0: T,
    beta_1: T
}

impl<T> LevenbergMarquardt<T>
{
    pub fn new(iters: u64, rho_minus: T, rho_plus: T) -> LevenbergMarquardt<T>
    {
        LevenbergMarquardt
        {
            iters: iters,
            beta_0: rho_minus,
            beta_1: rho_plus
        }
    }
}

impl<T> LevenbergMarquardt<T>
    where T: Real
{
    /// Minimize function func
    ///
    /// # Arguments
    ///
    /// * 'func0': Function to be minimized
    /// * 'x_0': Initial guess for the minimum
    ///
    /// # Return
    ///
    /// local minimum
    pub fn minimize<F: Jacobian<T>>(self: &Self, func: &F, x_0: &Vector<T>) -> OptimResult<Vector<T>>
    {
        let mut x_n: Vector<T> = x_0.clone();
        let mut mu_n: T = T::from_f64(0.5).unwrap();
        //let mut lambda_n: T = T::from_f64(0.4).unwrap();
        for _n in 0..self.iters
        {

            let mut d_n: Vector<T>;
            loop
            {
                let jacobian_x_n: Matrix<T> = func.jacobian(&x_n);
                let f_x_n: Vector<T> = func.eval(&x_n);

                let p_n: Vector<T> =  -(jacobian_x_n.clone().transpose() * f_x_n.clone());
                let (_j_m, j_n) = jacobian_x_n.dim();
                let left_n: Matrix<T> = jacobian_x_n.clone().transpose() * jacobian_x_n.clone() + Matrix::one(j_n) * mu_n * mu_n;
                d_n = left_n.solve(&p_n).unwrap();

                let x_n_1 = x_n.clone() + d_n.clone();
                let f_x_n_1: Vector<T> = func.eval(&x_n_1);

                let numerator: T = f_x_n.dotp(&f_x_n) - f_x_n_1.dotp(&f_x_n_1);
                let term: Vector<T> = f_x_n.clone() + jacobian_x_n.clone() * d_n.clone();
                let denumerator: T = f_x_n.dotp(&f_x_n) - term.dotp(&term);

                let epsilon: T = numerator / denumerator;

                if epsilon < self.beta_0
                {
                    mu_n = mu_n * T::from_f64(2.0).unwrap();
                }
                else
                {
                    if epsilon > self.beta_1
                    {
                        mu_n = mu_n / T::from_f64(2.0).unwrap();
                    }
                    break;
                }
            }
            x_n = x_n + d_n;
        }

        return OptimResult::new(x_n);
    }
}