use crate::data::number_types::nonzero::NonzeroSigned;
use crate::data::number_types::nonzero::sign::Sign;
use num::Zero;
use std::cmp::Ordering;

impl NonzeroSigned for i32 {
    fn signum(&self) -> Sign {
        match self.cmp(&i32::zero()) {
            Ordering::Less => Sign::Negative,
            Ordering::Equal => panic!(),
            Ordering::Greater => Sign::Positive,
        }
    }
}
