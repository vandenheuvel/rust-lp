macro_rules! test_rational {
    ($t:ident, $test_module_name:ident, $in_t:ident) => {
        #[allow(unused_imports)]
        mod $test_module_name {
            use num::{One, Zero};
            use num::FromPrimitive;
            use num::BigInt;

            use crate::data::number_types::rational::{Rational32, Rational64, Rational128, RationalBig};
            use crate::{R32, R64, R128, RB};

            #[test]
            fn field_identities() {
                for i in -10..0 {
                    assert_eq!($t!(0, i), $in_t::zero());
                }
                for i in 1..10 {
                    assert_eq!($t!(0, i), $in_t::zero());
                }
                for i in -10..0 {
                    assert_eq!($t!(i, i), $in_t::one());
                }
                for i in 1..10 {
                    assert_eq!($t!(i, i), $in_t::one());
                }
            }

            #[test]
            #[should_panic]
            fn panic_divide_zero_by_zero() {
                let _result = $t!(0, 0);
            }

            #[test]
            #[should_panic]
            fn panic_divide_nonzero_by_zero() {
                let _result = $t!(3, 0);
            }

            #[test]
            fn eq() {
                assert_eq!($t!(3, 2), $t!(6, 4));
                assert_eq!($t!(0, 2), $t!(0, 5));
                assert_eq!($t!(0, 2), $t!(0));
            }

            #[test]
            fn add() {
                assert_eq!($t!(3, 2) + $t!(6, 4), $t!(3));
                assert_eq!($t!(0, 2) + $t!(0, 5), $t!(0, 3));

                let mut x = $t!(0);
                for _ in 0..1000 {
                    x = x + $t!(1);
                }
                assert_eq!(x, $t!(1000));
            }

            #[test]
            fn sub() {
                assert_eq!($t!(3, 2) - $t!(6, 4), $t!(0, 9));
                assert_eq!($t!(0, 2) - $t!(0, 5), $t!(0, 3));
            }

            #[test]
            fn mul() {
                assert_eq!($t!(3, 2) * $t!(6, 4), $t!(9, 4));
                assert_eq!($t!(0, 2) * $t!(0, 5), $t!(0, 3));
            }

            #[test]
            fn div() {
                assert_eq!($t!(3, 2) / $t!(6, 4), $in_t::one());
                assert_eq!($t!(0, 2) / $t!(2, 5), $in_t::zero());
            }

            #[test]
            #[should_panic]
            fn div_zero() {
                let _result = $t!(4564, 65468) / $t!(0, 654654);
            }
        }
    };
}

test_rational!(R32, test_32, Rational32);
test_rational!(R64, test_64, Rational64);
test_rational!(R128, test_128, Rational128);
test_rational!(RB, test_big, RationalBig);