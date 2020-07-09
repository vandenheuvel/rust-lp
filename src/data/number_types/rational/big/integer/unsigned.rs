use std::cmp::Ordering;
use std::fmt;
use std::ops::{Add, Mul};

use num_traits::{One, Zero};
use smallvec::SmallVec;

use crate::data::number_types::rational::big::integer::{UnsignedDigit, cmp_slice};

/// A big unsigned integer type.
#[derive(Debug)]
pub struct Unsigned {
    data: SmallVec<[UnsignedDigit; 4]>, // TODO(ENHANCEMENT): Generics over inline size
}

impl PartialEq<Unsigned> for Unsigned {
    fn eq(&self, other: &Unsigned) -> bool {
        // TODO: Why these debug asserts?
        debug_assert!(self.data.last() != Some(&0));
        debug_assert!(other.data.last() != Some(&0));

        self.data == other.data
    }
}
impl Eq for Unsigned {}

impl PartialOrd for Unsigned {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for Unsigned {
    fn cmp(&self, other: &Self) -> Ordering {
        cmp_slice(&self.data, &other.data)
    }
}

impl fmt::Display for Unsigned {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        unimplemented!()
    }
}

impl Add for Unsigned {
    type Output = Self;

    fn add(self, rhs: Self) -> Self::Output {
        unimplemented!()
    }
}

impl Mul for Unsigned {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        unimplemented!()
    }
}

impl Zero for Unsigned {
    fn zero() -> Self {
        Self { data: Default::default(), }
    }

    fn set_zero(&mut self) {
        self.data.clear()
    }

    fn is_zero(&self) -> bool {
        self.data.is_empty()
    }
}

impl One for Unsigned {
    fn one() -> Self {
        Self { data: smallvec![1], }
    }

    fn set_one(&mut self) {
        self.data.clear();
        self.data.push(1);
    }

    fn is_one(&self) -> bool {
        self.data[..] == [1]
    }
}
