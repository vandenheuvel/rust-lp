use smallvec::SmallVec;
use crate::data::number_types::rational::big::integer::{UnsignedDigit, SignedDigit};
use num_traits::Zero;
use std::ops::Add;

#[derive(Debug, Clone)]
pub struct Signed {
    first: SignedDigit,
    data: SmallVec<[UnsignedDigit; 4 - 1]>, // TODO(ENHANCEMENT): Generics over inline size
}

impl Add for Signed {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        if self.is_zero() {
            return rhs;
        }
        if rhs.is_zero() {}
    }
}

impl Zero for Signed {
    fn zero() -> Self {
        Self {
            first: 0,
            data: SmallVec::new(),
        }
    }

    fn set_zero(&mut self) {
        self.first = 0;
        self.data.clear();
    }

    fn is_zero(&self) -> bool {
        self.first == 0 && self.data.is_empty()
    }
}
