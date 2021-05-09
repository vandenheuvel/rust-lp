use std::ops::{Neg, Mul};

use num::Zero;

use crate::data::number_types::nonzero::Nonzero;
use std::cmp::Ordering;

/// A signed number that can have a nonzero value.
pub trait NonzeroSigned: Nonzero + Clone {
    /// Whether the value is positive or negative.
    fn signum(&self) -> Sign;
    /// Whether `x > 0`.
    fn is_positive(&self) -> bool {
        self.signum() == Sign::Positive
    }
    /// Whether `x < 0`.
    fn is_negative(&self) -> bool {
        self.signum() == Sign::Negative
    }
}

/// Sign of a nonzero value.
///
/// Existing `Sign` traits, such in `num`, typically have a third value for the sign of 0. Working
/// with that trait creates many branches or match cases that should never be possible.
#[derive(Eq, PartialEq, Copy, Clone, Debug)]
pub enum Sign {
    /// `x > 0`
    Positive,
    /// `x < 0`
    Negative,
}

impl<T: Zero + Nonzero + PartialOrd<Self> + Clone> NonzeroSigned for T {
    default fn signum(&self) -> Sign {
        debug_assert!(self.is_not_zero());

        match self.partial_cmp(&Self::zero()) {
            Some(Ordering::Less) => Sign::Negative,
            Some(Ordering::Greater) => Sign::Positive,
            Some(Ordering::Equal) | None => unreachable!("\
                Should only be used on nonzero values, and those should always be comparable with \
                the zero value of the type.\
            "),
        }
    }
}

impl PartialOrd for Sign {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (self, other) {
            (Sign::Positive, Sign::Positive) => None,
            (Sign::Positive, Sign::Negative) => Some(Ordering::Greater),
            (Sign::Negative, Sign::Positive) => Some(Ordering::Less),
            (Sign::Negative, Sign::Negative) => None,
        }
    }
}

impl Mul for Sign {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        match (self, rhs) {
            (Sign::Positive, Sign::Positive) => Sign::Positive,
            (Sign::Positive, Sign::Negative) => Sign::Negative,
            (Sign::Negative, Sign::Positive) => Sign::Negative,
            (Sign::Negative, Sign::Negative) => Sign::Positive,
        }
    }
}
