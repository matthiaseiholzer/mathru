//! Linear algebra
//!

pub use self::{matrix::Matrix, vector::Vector};

#[macro_use]
pub mod vector;
#[macro_use]
pub mod matrix;
#[cfg(feature = "lapack")]
pub mod lapack;
#[cfg(feature = "lapack")]
pub mod blas;
