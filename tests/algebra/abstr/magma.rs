macro_rules! test_magma
{
    ($a:expr, $b:expr, $(($id:ident, $s:ty)),*) =>
    {
        $(
            mod $id
            {
                use mathru::algebra::abstr::{Magma, Addition, Multiplication};

                #[test]
                fn magma_addition()
                {
                    let a: $s = $a;
                    let b: $s = $b;

                    assert_eq!(a + b, Magma::<Addition>::operate(a, b));
                }

                #[test]
                fn magma_multiplication()
                {
                    let a: $s = $a;
                    let b: $s = $b;

                    assert_eq!(a * b, Magma::<Multiplication>::operate(a, b));
                }
            }
        )*
    };
}

test_magma!(5,
            2,
            (u8, u8),
            (u16, u16),
            (u32, u32),
            (u64, u64),
            (u128, u128));
test_magma!(-5,
            2,
            (i8, i8),
            (i16, i16),
            (i32, i32),
            (i64, i64),
            (i128, i128));
test_magma!(5.0, 2.0, (f32, f32), (f64, f64));
