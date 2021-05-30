
// // macro_rules! test_monoid
// // {
// //     ($a:expr, $b:expr, $(($id:ident, $s:ty)),*) =>
// //     {
// //         $(
// //             mod $id
// //             {
// //                 use mathru::algebra::abstr::{Monoid, Addition, Multiplication};
// //
// //                 #[test]
// //                 fn test_monoid_addition()
// //                 {
// //                     let a: $s = $a;
// //                     let b: $s = $b;
// //
// //                     assert_eq!(a + Monoid::<Addition>::id(), Monoid::<Addition>::operate(a, b));
// //                 }
// //
// //                 #[test]
// //                 fn test_monoid_multiplication()
// //                 {
// //                     let a: $s = $a;
// //                     let b: $s = $b;
// //
// //                     assert_eq!(a * Monoid::<Multiplication>::id(), Monoid::<Multiplication>::operate(a, b));
// //                 }
// //             }
// //         )*
// //     };
// // }
// //
// // test_monoid!(5,
// //             0,
// //             (u8, u8),
// //             (u16, u16),
// //             (u32, u32),
// //             (u64, u64),
// //             (u128, u128));
// // test_monoid!(-5,
// //             0,
// //             (i8, i8),
// //             (i16, i16),
// //             (i32, i32),
// //             (i64, i64),
// //             (i128, i128));
// // test_monoid!(5.0, 0.0, (f32, f32), (f64, f64));
//
// use mathru::algebra::abstr::{Monoid, Addition};
//
// #[test]
// fn test_monoid_addition()
// {
//     let a: u8 = 2;
//     let b: u8 = 0;
//
//     assert_eq!(a + b, <Semigroup<Addition> as Magma::<Addition>>::operate(a, b));
// }