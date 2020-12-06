//! # Rational numbers
//!
//! Useful for testing of methods and overall correctness of implementations.
//!
//! It appears that in practice, no fast solvers use rational numbers.
use num::BigInt;
use num::rational::Ratio;

use crate::data::number_types::traits::{Field, FieldRef, OrderedField, OrderedFieldRef};

macro_rules! impl_rational {
    ($in_t:ident) => {
        impl OrderedField for Ratio<$in_t> {
        }

        /// Except for the finiteness of the type, this is correct.
        impl Field for Ratio<$in_t> {
        }

        impl OrderedFieldRef<Ratio<$in_t>> for &Ratio<$in_t> {
        }

        impl FieldRef<Ratio<$in_t>> for &Ratio<$in_t> {
        }
    }
}
impl_rational!(i32);
impl_rational!(i64);
impl_rational!(i128);
impl_rational!(BigInt);

/// Shorthand for creating a rational number in tests.
#[macro_export]
macro_rules! R32 {
    ($value:expr) => {
        Ratio::<i32>::from_f64($value as f64).unwrap()
    };
    ($numer:expr, $denom:expr) => {
        Ratio::<i32>::new($numer, $denom)
    };
}
/// Shorthand for creating a rational number in tests.
#[macro_export]
macro_rules! R64 {
    ($value:expr) => {
        Ratio::<i64>::from_f64($value as f64).unwrap()
    };
    ($numer:expr, $denom:expr) => {
        Ratio::<i64>::new($numer, $denom)
    };
}
/// Shorthand for creating a rational number in tests.
#[macro_export]
macro_rules! R128 {
    ($value:expr) => {
        Ratio::<i128>::from_f64($value as f64).unwrap()
    };
    ($numer:expr, $denom:expr) => {
        Ratio::<i128>::new($numer, $denom)
    };
}
/// Shorthand for creating a rational number in tests.
#[macro_export]
macro_rules! BR {
    ($value:expr) => {
        Ratio::<BigInt>::from_f64($value as f64).unwrap()
    };
    ($numer:expr, $denom:expr) => {
        Ratio::<BigInt>::new(BigInt::from_i64($numer).unwrap(), BigInt::from_i64($denom).unwrap())
    };
}

#[cfg(test)]
mod test {
    macro_rules! test_rational {
        ($t:ident, $test_module_name:ident, $in_t:ident) => {
            #[allow(unused_imports)]
            mod $test_module_name {
                use num::rational::Ratio;
                use num::{One, Zero};
                use num::FromPrimitive;
                use num::BigInt;

                #[test]
                fn field_identities() {
                    for i in -10..0 {
                        assert_eq!($t!(0, i), Ratio::<$in_t>::zero());
                    }
                    for i in 1..10 {
                        assert_eq!($t!(0, i), Ratio::<$in_t>::zero());
                    }
                    for i in -10..0 {
                        assert_eq!($t!(i, i), Ratio::<$in_t>::one());
                    }
                    for i in 1..10 {
                        assert_eq!($t!(i, i), Ratio::<$in_t>::one());
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
                    assert_eq!($t!(3, 2) / $t!(6, 4), Ratio::<$in_t>::one());

                    assert_eq!($t!(0, 2) / $t!(2, 5), Ratio::<$in_t>::zero());
                }

                #[test]
                #[should_panic]
                fn div_zero() {
                    let _result = $t!(4564, 65468) / $t!(0, 654654);
                }
            }
        };
    }
    test_rational!(R32, test_i32, i32);
    test_rational!(R64, test_i64, i64);
    test_rational!(R128, test_i128, i128);
    test_rational!(BR, test_br, BigInt);
}
