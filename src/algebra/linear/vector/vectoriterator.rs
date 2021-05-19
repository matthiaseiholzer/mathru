use crate::algebra::{
    linear::matrix::MatrixIterator,
};

pub struct VectorIterator<'a, T>
{
    iter: MatrixIterator<'a, T>,
}

impl<'a, T> VectorIterator<'a, T>
{
    pub fn new(iter: MatrixIterator<'a, T>) -> VectorIterator<'a, T>
    {
        VectorIterator{iter}
    }
}

impl<'a, T> Iterator for VectorIterator<'a, T>
{
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item>
    {
        self.iter.next()
    }
}
