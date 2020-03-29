#[cfg(test)]
mod implicit_euler
{
	use mathru::algebra::linear::{Vector};
	use mathru::analysis::ode::{ImplicitEuler, ODEProblem};
	use super::super::problem::{ImplicitODE1};

	fn compare_epsilon(a: f64, b: f64, epsilon: f64) -> bool
    {
    	if (a - b).abs() > epsilon
        {
        	println!("|a-b|: {}", (a-b).abs());
        	return false;
        }

        return true;
    }

	#[test]
	fn fn1()
	{
		let problem: ImplicitODE1 = ImplicitODE1::default();
		let solver: ImplicitEuler<f64> = ImplicitEuler::new(0.001);

		let (t, y): (Vec<f64>, Vec<Vector<f64>>) = solver.solve(&problem).unwrap();

		let len: usize = y.len();

		let time_span: (f64, f64) = problem.time_span();
		let init_cond: Vector<f64> = problem.init_cond();

		assert!(compare_epsilon(time_span.1, t[len - 1], 0.000000001));
		assert!(compare_epsilon(2.0 - (-4.0 * time_span.1).exp() , *y[len - 1].get(0), 0.00002));
	}
}