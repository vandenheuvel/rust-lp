//! # An arbitrary precision rational type
//!
//! At the moment, this is just wrapping the `num::BigRational` type, following the newtype pattern.
//! This is needed because some of the impl's in this module are not provided by `num`. Methods on
//! this type can be modified and specialized as needed.
use std::fmt::{Debug, Formatter};
use std::fmt;

use num::Zero;

use crate::data::number_types::nonzero::Nonzero;

mod with_fixed_size_integer;
mod with_binary;
mod wrapping;
mod special;
#[cfg(test)]
mod test;

/// A big rational type that currently completely relies on the methods of `num::BigRational`.
#[derive(
    Clone,
    Ord, PartialOrd, Eq, PartialEq,
)]
pub struct Big(pub num::BigRational);

impl Big {
    /// Create a new instance by converting the two provided numbers into arbitrary size ints.
    #[must_use]
    pub fn new(numer: i64, denom: i64) -> Self {
        Self(num::BigRational::new(numer.into(), denom.into()))
    }
}

impl Nonzero for Big {
    fn is_not_zero(&self) -> bool {
        !self.is_zero()
    }
}

impl Debug for Big {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let Self(inner) = self;

        write!(f, "{}", inner.numer().to_string())?;
        write!(f, "/")?;
        write!(f, "{}", inner.denom().to_string())
    }
}