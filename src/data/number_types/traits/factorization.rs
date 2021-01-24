use crate::data::number_types::nonzero::{NonzeroSigned, Nonzero, NonzeroSign};

pub trait NonzeroFactorizable: NonzeroSigned {
    /// Some prime greater than 1.
    type Factor: Nonzero + Ord;
    /// How often the factor appears in the number.
    type Power: Nonzero;

    fn factorize(&self) -> NonzeroFactorization<Self::Factor, Self::Power>;
}

/// Prime factorization representation of a nonzero rational number.
///
/// Includes a sign.
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
