macro_rules! test_semigroup
{
    ($a:expr, $b:expr, $c:expr, $(($id:ident, $s:ty)),*) =>
    {
        $(
            mod $id
            {
                use mathru::algebra::abstr::{Semigroup, Addition, Multiplication};

                #[test]
                fn semigroup_addition()
                {
                    assert!(Semigroup::<Addition>::is_associative($a, $b, $c));
                }

                #[test]
                fn semigroup_multiplication()
                {
                    assert!(Semigroup::<Multiplication>::is_associative($a, $b, $c));
                }
            }
        )*
    };
}

test_semigroup!(5,
            2,
            4,
            (u8, u8),
            (u16, u16),
            (u32, u32),
            (u64, u64),
            (u128, u128));
test_semigroup!(-5,
            2,
            3,
            (i8, i8),
            (i16, i16),
            (i32, i32),
            (i64, i64),
            (i128, i128));
test_semigroup!(5.0, 2.0, 3.0, (f32, f32), (f64, f64));
