# Mathru

[![crate](https://img.shields.io/crates/v/mathru.svg)](https://crates.io/crates/mathru)
[![documentation](https://docs.rs/mathru/badge.svg)](https://docs.rs/mathru)
![minimum rustc 1.51.0](https://img.shields.io/badge/rustc-1.51.0-green.svg)
![maintenance](https://img.shields.io/badge/maintenance-actively--developed-brightgreen.svg)
[![pipeline status](https://gitlab.com/matthiaseiholzer/mathru/badges/main/pipeline.svg)](https://gitlab.com/matthiaseiholzer/mathru/-/commits/main)
------------
Mathru is a numeric library containing algorithms for linear algebra, analysis and statistics written in pure Rust with BLAS/LAPACK support.

## Features
The following features are implemented in this create:
* [Algebra](https://matthiaseiholzer.gitlab.io/mathru/documentation/algebra/)
    * [Abstract](https://matthiaseiholzer.gitlab.io/mathru/documentation/algebra/abstract/)
        * [Polynomial](https://matthiaseiholzer.gitlab.io/mathru/documentation/algebra/abstract/polynomial/)
    * [Linear algebra](https://matthiaseiholzer.gitlab.io/mathru/documentation/algebra/linear/)
        * [Vector](https://matthiaseiholzer.gitlab.io/mathru/documentation/algebra/linear/vector/)
        * [Matrix](https://matthiaseiholzer.gitlab.io/mathru/documentation/algebra/linear/matrix/)
            * Basic matrix operations(+,-,*)
            * Transposition (In-place)
            * [LU decomposition](https://matthiaseiholzer.gitlab.io/mathru/documentation/algebra/linear/matrix/#lu-with-partial-pivoting)
            * [QR decomposition](https://matthiaseiholzer.gitlab.io/mathru/documentation/algebra/linear/matrix/#qr)
            * [Hessenberg decomposition](https://matthiaseiholzer.gitlab.io/mathru/documentation/algebra/linear/matrix/#hessenberg)
            * [Cholesky decomposition](https://matthiaseiholzer.gitlab.io/mathru/documentation/algebra/linear/matrix/#cholesky)
            * Eigen decomposition
            * Singular value decomposition
            * Inverse
            * Pseudo inverse
            * Determinant
            * Trace
            * [Solve linear system](https://matthiaseiholzer.gitlab.io/mathru/documentation/algebra/linear/matrix/#linear-system-resolution)

* [Ordinary differential equation (ODE)](https://matthiaseiholzer.gitlab.io/mathru/documentation/analysis/differentialeq/)
    * [Explicit methods](https://matthiaseiholzer.gitlab.io/mathru/documentation/analysis/differentialeq/ode/explicit/)
        * Euler method
        * Heun's 2nd order method
        * Midpoint method
        * Ralston's 2nd order method
        * Kutta 3rd order
        * Heun's 3rd order method
        * Ralston's 3rd order method
        * Runge-Kutta 4th order
        * Runge-Kutta-Felhberg
        * Dormand-Prince
        * Cash-Karp
        * Tsitouras
        * Bogacki-Shampine
        * Adams-Bashforth
    * Automatic step size control with starting step size
    * [Implicit methods](https://matthiaseiholzer.gitlab.io/mathru/documentation/analysis/differentialeq/ode/implicit)
        * Implicit Euler
        * Backward differentiation formula (BDF)

* [Optimization](https://matthiaseiholzer.gitlab.io/mathru/documentation/optimization)
    * Gauss-Newton algorithm
    * Gradient descent
    * Newton method
    * Levenberg-Marquardt algorithm
    * Conjugate gradient method

* [Statistics](https://matthiaseiholzer.gitlab.io/mathru/documentation/statistics)
    * Probability distribution
        * Bernoulli
        * Beta
        * Binomial
        * Exponential
        * Gamma
        * Chi-squared
        * [Normal](https://matthiaseiholzer.gitlab.io/mathru/documentation/statistics/distribution/normal/)
        * Log-Normal
        * Poisson
        * Raised cosine
        * Student-t
        * Uniform
    * [Test](https://matthiaseiholzer.gitlab.io/mathru/documentation/statistics/test/)
        * Chi-squared
        * G
        * [Student-t](https://matthiaseiholzer.gitlab.io/mathru/documentation/statistics/test/t_test/)

* Elementary functions
    * trigonometric functions
    * hyperbolic functions
    * exponential functions

* [special functions](https://matthiaseiholzer.gitlab.io/mathru/documentation/special)
    * [gamma functions](https://matthiaseiholzer.gitlab.io/mathru/documentation/special/gamma/)
        * gamma function
        * log-gamma function
        * digamma function
        * upper incomplete gamma function
        * upper incomplete regularized gamma function
        * inverse upper incomplete regularized gamma function
        * lower incomplete gamma function
        * lower regularized incomplete gamma function
        * inverse lower regularized incomplete gamma function
    * [beta functions](https://matthiaseiholzer.gitlab.io/mathru/documentation/special/beta/)
        * beta function
        * incomplete beta function
        * incomplete regularized beta function
    * [error functions](https://matthiaseiholzer.gitlab.io/mathru/documentation/special/error/)
        * error function
        * complementary error function
        * inverse error function
        * inverse complementary error function
    * hypergeometric functions

## Usage

Add this to your `Cargo.toml` for the native Rust implementation:

```toml
[dependencies.mathru]
version = "0.10"
```
Add the following lines to 'Cargo.toml' if the openblas library should be used:

```toml
[dependencies.mathru]
version = "0.10"
default-features = false
features = "openblas"
```

One of the following implementations for linear algebra can be activated as a feature:
- native: Native Rust implementation(activated by default)
- [openblas](https://www.openblas.net): Optimized BLAS library
- [netlib](https://www.netlib.org): Collection of mathematical software, papers, and databases
- [intel-mkl](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html): Intel Math Kernel Library
- [accelerate](https://developer.apple.com/documentation/accelerate) Make large-scale mathematical computations and image calculations, optimized for high performance and low-energy consumption.(macOS only)


Then import the modules and it is ready to be used.

```rust
use mathru as mr;
```

### Solve a system of linear equations

```rust
use mathru::{
    algebra::linear::{
        matrix::{LUDec, Solve},
        Matrix, Vector,
    },
    matrix, vector,
};

/// Solves a system of linear equations
fn main()
{
    let a: Matrix<f64> = matrix![6.0, 2.0, -1.0; -3.0, 5.0, 3.0; -2.0, 1.0, 3.0];
    let b: Vector<f64> = vector![48.0; 49.0; 24.0];

    // Decompose a into a lower and upper matrix
    let lu_dec: LUDec<f64> = a.dec_lu().unwrap();

    // Solve the system of linear equations with the decomposed matrix
    let _x1: Vector<f64> = lu_dec.solve(&b).unwrap();

    // Solve it directly
    let _x2: Vector<f64> = a.solve(&b).unwrap();
}
```

### Solve ordinary differential equation with the Dormand-Prince algorithm

```rust
use mathru::{
    algebra::linear::Vector,
    analysis::differential_equation::ordinary::{problem, DormandPrince54, ExplicitODE},
};

fn main()
{
    // Create an ODE instance
    let problem: problem::Euler<f64> = problem::Euler::default();

    let (x_start, x_end) = problem.time_span();

    // Create a ODE solver instance
    let h_0: f64 = 0.0001;
    let n_max: u32 = 800;
    let abs_tol: f64 = 10e-7;

    let solver: DormandPrince54<f64> = DormandPrince54::new(abs_tol, h_0, n_max);

    // Solve ODE
    let (x, y): (Vec<f64>, Vec<Vector<f64>>) = solver.solve(&problem).unwrap();
}
```

### Further examples

For further examples, see [project page](https://matthiaseiholzer.gitlab.io/mathru)

## Documentation

See [project page](https://matthiaseiholzer.gitlab.io/mathru) for more information and examples.
The API is documented on [docs.rs](https://docs.rs/mathru).

## License

Licensed under

 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

### Contribution

Any contribution is welcome!