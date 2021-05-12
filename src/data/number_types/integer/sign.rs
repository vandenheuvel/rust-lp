//! # Nonzero signs of integers
use std::cmp::Ordering;

use num::Zero;

use crate::data::number_types::nonzero::NonzeroSigned;
use crate::data::number_types::nonzero::sign::Sign;

impl NonzeroSigned for i32 {
    fn signum(&self) -> Sign {
        match self.cmp(&i32::zero()) {
            Ordering::Less => Sign::Negative,
            Ordering::Equal => panic!(),
            Ordering::Greater => Sign::Positive,
        }
    }
}
