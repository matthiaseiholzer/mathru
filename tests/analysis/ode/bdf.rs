use mathru::{
    algebra::linear::Vector,
    analysis::differential_equation::ordinary::{problem, BDF},
};

#[test]
fn fn1()
{
    let problem: problem::Euler<f64> = problem::Euler::default();
    let solver: BDF<f64> = BDF::new(6, 0.001);

    let (_x, y): (Vec<f64>, Vec<Vector<f64>>) = solver.solve(&problem).unwrap();

    assert_relative_eq!(0.988, *y.last().unwrap().get(0), epsilon=0.001);
}
