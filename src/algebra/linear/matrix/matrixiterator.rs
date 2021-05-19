use std::slice::Iter;

pub struct MatrixIterator<'a, T>
{
    iter: Iter<'a, T>,
}

impl<'a, T> MatrixIterator<'a, T>
{
    pub fn new(iter: Iter<'a, T>) -> MatrixIterator<'a, T>
    {
        MatrixIterator{ iter}
    }
}

impl<'a, T> Iterator for MatrixIterator<'a, T>
{
    type Item = &'a T;

    fn next(&mut self) -> Option<Self::Item>
    {
        return self.iter.next();
    }
}
