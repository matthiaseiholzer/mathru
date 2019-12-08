#[macro_use]
pub mod matrix;
mod matrixintoiterator;
mod matrixiterator;
mod matrixiteratormut;
mod matrixrowiterator;
mod matrixrowiteratormut;
mod matrixcolumniterator;
mod matrixcolumniteratormut;

mod eigenvalue;
mod hessenberg;
mod lu;
mod qr;
mod inverse;
mod mul;
mod add;
mod sub;
mod div;
mod cholesky;
mod solve;
mod substitute;


pub use self::substitute::Substitute;
pub use self::cholesky::CholeskyDec;
pub use self::lu::LUDec;
pub use self::solve::Solve;
pub use self::inverse::Inverse;
pub use self::matrix::Matrix;
pub use self::matrixintoiterator::MatrixIntoIterator;
pub use self::matrixiterator::MatrixIterator;
pub use self::matrixiteratormut::MatrixIteratorMut;
pub use self::matrixrowiterator::MatrixRowIterator;
pub use self::matrixrowiteratormut::MatrixRowIteratorMut;
pub use self::matrixcolumniterator::MatrixColumnIterator;
pub use self::matrixcolumniteratormut::MatrixColumnIteratorMut;