use std::slice::IterMut;

pub struct MatrixIteratorMut<'a, T>
{
    iter: IterMut<'a, T>,
}

impl<'a, T> MatrixIteratorMut<'a, T>
{
    pub fn new(iter: IterMut<'a, T>) -> MatrixIteratorMut<'a, T>
    {
        MatrixIteratorMut{ iter }
    }
}
impl<'a, T> Iterator for MatrixIteratorMut<'a, T>
{
    type Item = &'a mut T;

    fn next(&mut self) -> Option<Self::Item>
    {
        return self.iter.next();
    }
}
