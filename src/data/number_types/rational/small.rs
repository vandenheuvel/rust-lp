//! # Rational types of fixed size
use std::cmp::Ordering;
use std::num::{NonZeroI128, NonZeroI16, NonZeroI32, NonZeroI64, NonZeroI8, NonZeroU128, NonZeroU16,
               NonZeroU32, NonZeroU64, NonZeroU8};
use std::ops::Add;

use num::Zero;

use crate::algorithm::two_phase::tableau::kind::artificial::Cost;
use crate::data::number_types::nonzero::{Nonzero, NonzeroSign};
use crate::data::number_types::nonzero::sign::{NonzeroSigned, Sign};
use crate::data::number_types::rational::Ratio;

/// Aliased type to ease a possible transition to own variant in the future.
pub type Rational32 = num::rational::Rational32;
/// Aliased type to ease a possible transition to own variant in the future.
pub type Rational64 = num::rational::Rational64;
/// Aliased type to ease a possible transition to own variant in the future.
pub type Rational128 = num::rational::Ratio<i128>;

macro_rules! rational {
    ($name:ident, $numerator:ty, $denominator:ty) => {
        pub type $name = Ratio<$numerator, $denominator>;

        impl Nonzero for $name {
            fn is_not_zero(&self) -> bool {
                self.numerator != 0
            }
        }

        impl NonzeroSigned for $name {
            fn signum(&self) -> NonzeroSign {
                match self.numerator.cmp(&0) {
                    Ordering::Less => NonzeroSign::Negative,
                    Ordering::Equal => panic!("NonzeroSigned::signum can't be called on a zero value"),
                    Ordering::Greater => NonzeroSign::Positive,
                }
            }
        }

        impl Add<Cost> for $name {
            type Output = Self;

            fn add(mut self, rhs: Cost) -> Self::Output {
                match rhs {
                    Cost::Zero => self,
                    Cost::One => {
                        self.numerator += self.denominator as $numerator;
                        self
                    }
                }
            }
        }
    }
}

impl Zero for XRational32 {
    fn zero() -> Self {
        Self { numerator: 0, denominator: 1 }
    }

    fn set_zero(&mut self) {
        self.numerator = 0;
        self.denominator = 1;
    }

    fn is_zero(&self) -> bool {
        self.numerator == 0
    }
}

impl Add for XRational32 {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {

    }
}

rational!(XRational8, i8, NonZeroU8);
rational!(XRational16, i16, NonZeroU16);
rational!(XRational32, i32, NonZeroU32);
rational!(XRational64, i64, NonZeroU64);
rational!(XRational128, i128, NonZeroU128);

rational!(XNonzeroRational8, NonZeroI8, NonZeroU8);
rational!(XNonnegativeRational8, u8, NonZeroU8);
rational!(XPositiveRational8, NonZeroU8, NonZeroU8);
