//! # Rational types of fixed size
use std::ops::Add;

use crate::algorithm::two_phase::tableau::kind::artificial::Cost;
use crate::data::number_types::rational::{Ratio};
use crate::data::number_types::nonzero::Nonzero;
use crate::data::number_types::nonzero::sign::{NonzeroSigned, Sign};
use num::Zero;

/// Aliased type to ease a possible transition to own variant in the future.
pub type Rational32 = num::rational::Rational32;
/// Aliased type to ease a possible transition to own variant in the future.
pub type Rational64 = num::rational::Rational64;
/// Aliased type to ease a possible transition to own variant in the future.
pub type Rational128 = num::rational::Ratio<i128>;


impl Nonzero for Rational32 {
    fn is_not_zero(&self) -> bool {
        !self.is_zero()
    }
}

impl Nonzero for Rational64 {
    fn is_not_zero(&self) -> bool {
        !self.is_zero()
    }
}

impl Nonzero for Rational128 {
    fn is_not_zero(&self) -> bool {
        !self.is_zero()
    }
}

pub type R32 = Ratio<i32, u32>;

impl NonzeroSigned for R32 {
    fn signum(&self) -> Sign {
        NonzeroSigned::signum(&self.numerator)
    }
}

impl Add<Cost> for Rational64 {
    type Output = Self;

    fn add(self, rhs: Cost) -> Self::Output {
        // TODO(PERFORMANCE): Ensure that this is fast
        match rhs {
            Cost::Zero => self,
            Cost::One => {
                Self::new(self.numer() + self.denom(), *self.denom())
            }
        }
    }
}
