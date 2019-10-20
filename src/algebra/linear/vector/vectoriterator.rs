use crate::algebra::abstr::Real;
use crate::algebra::linear::matrix::MatrixIterator;

pub struct VectorIterator<'a, T>
{
    pub iter: MatrixIterator<'a, T>
}

impl<'a, T> Iterator for VectorIterator<'a, T>
    where T: Real
{
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item>
    {
        self.iter.next()
    }
}