use crate::algebra::abstr::{Field, Scalar};
use std::{iter::Iterator, vec::IntoIter};

pub struct MatrixIntoIterator<T>
{
    pub iter: IntoIter<T>,
}

impl<T> Iterator for MatrixIntoIterator<T> where T: Field + Scalar
{
    type Item = T;

    fn next(self: &mut Self) -> Option<Self::Item>
    {
        return self.iter.next();
    }
}
