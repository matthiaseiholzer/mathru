use crate::algebra::abstr::{Field, Scalar, Zero};
use crate::algebra::linear::{Vector};
use crate::vector;
use std::slice::IterMut;

pub struct MatrixRowIteratorMut<'a, T>
{
    iter: IterMut<'a, T>,
}

impl<'a, T> MatrixRowIteratorMut<'a, T>
    where T: Zero
{
    pub fn new(iter: IterMut<'a, T>) -> MatrixRowIteratorMut<'a, T>
    {
        return MatrixRowIteratorMut{iter};
    }
}

impl<'a, T> Iterator for MatrixRowIteratorMut<'a, T>
{
    type Item = &'a mut Vector<T>;

    fn next(&'a mut self) -> Option<Self::Item>
    {
        None
    }
}
