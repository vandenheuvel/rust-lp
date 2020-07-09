use std::ops::{Add, AddAssign, Mul, MulAssign};

use num_traits::One;
use num_traits::one;
use smallvec::SmallVec;
use crate::data::number_types::rational::big::integer::{NonZeroUnsignedDigit, UnsignedDigit};

#[derive(Eq, PartialEq, Clone, Debug)]
pub struct NonZeroUnsigned {
    lowest: NonZeroUnsignedDigit,
    low_to_high: SmallVec<[UnsignedDigit; 4 - 1]>, // minus one, because lowest already is inline
}

impl NonZeroUnsigned {
    fn gcd(&self, other: &Self) -> Self {

    }
}

impl Mul for NonZeroUnsigned {
    type Output = Self;

    fn mul(mut self, rhs: Self) -> Self::Output {
        self *= rhs;
        self
    }
}
impl MulAssign for NonZeroUnsigned {
    fn mul_assign(&mut self, rhs: Self) {
        unimplemented!()
    }
}

impl One for NonZeroUnsigned {
    fn one() -> Self {
        Self { lowest: one(), low_to_high: Default::default(), }
    }
}
