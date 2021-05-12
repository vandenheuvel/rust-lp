//! # Number factorization
//!
//! Factorize integers and rational numbers into numbers that are often primes.
use std::convert::identity;
use std::ops::{Add, Mul};

use num::One;

use crate::algorithm::utilities::merge_sparse_indices;
use crate::data::number_types::nonzero::{Nonzero, NonzeroSign};

/// Creating a factorization of an integer or rational number.
pub trait NonzeroFactorizable: Nonzero {
    /// Some prime greater than 1.
    type Factor: Nonzero + Ord + Clone;
    /// How often the factor appears in the number.
    ///
    /// This is marked Copy, because a 64-bit power already allows for values up to 2^(2^64), which
    /// has about 5.6 * 10^18 decimal digits.
    type Power: Nonzero + Copy + Clone;

    /// Decompose into factors.
    ///
    /// Note that these factors will often be, but are not guaranteed to be, primes.
    fn factorize(&self) -> NonzeroFactorization<Self::Factor, Self::Power>;
}

/// Prime factorization representation of a nonzero rational number.
///
/// Includes a sign.
#[derive(Eq, PartialEq, Clone, Debug)]
pub struct NonzeroFactorization<Factor, Power> {
    /// Whether the number is negative.
    pub sign: NonzeroSign,
    /// `(prime factor, power)` tuples.
    ///
    /// The factors should all be smaller than 64 bits and can have negative powers; that is, appear
    /// in the denominator. The powers can't be zero, as this is a sparse representation.
    ///
    /// When this field is empty, the value `1` or `-1` is represented, depending on `sign`.
    pub factors: Vec<(Factor, Power)>,
}

impl<Factor: Ord, Power: Nonzero + Eq + Add<Output=Power>> One for NonzeroFactorization<Factor, Power> {
    fn one() -> Self {
        Self {
            sign: NonzeroSign::Positive,
            factors: vec![],
        }
    }
}

impl<Factor: Ord, Power: Nonzero + Eq + Add<Output=Power>> Mul for NonzeroFactorization<Factor, Power> {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let sign = self.sign * rhs.sign;
        let factors = merge_sparse_indices(
            self.factors.into_iter(), rhs.factors.into_iter(),
            Add::add, identity, identity,
        );

        Self { sign, factors }
    }
}

#[cfg(test)]
mod test {
    use crate::data::number_types::traits::factorization::NonzeroFactorizable;

    #[test]
    fn test_multiply() {
        let a = 6;
        let b = 36;

        let a_factorized = a.factorize();
        let b_factorized = b.factorize();

        let c = a * b;
        let c_factorized = c.factorize();

        assert_eq!(a_factorized * b_factorized, c_factorized);
    }
}
