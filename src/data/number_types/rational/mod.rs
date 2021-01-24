//! # Rational numbers
//!
//! Primary way to do arbitrary precision computation.
pub use big::Big as RationalBig;
pub use small::Rational128 as Rational128;
pub use small::Rational32 as Rational32;
pub use small::Rational64 as Rational64;

use crate::data::number_types::nonzero::Nonzero;

mod small;
mod big;
mod macros;
mod factorization;

pub trait Rational: Nonzero {
    type Numerator;
    type Denominator;

    fn numerator(&self) -> &Self::Numerator;
    fn denominator(&self) -> &Self::Denominator;

    fn numerator_mut(&mut self) -> &mut Self::Numerator;
    fn denominator_mut(&mut self) -> &mut Self::Denominator;
}

#[cfg(test)]
mod test;
