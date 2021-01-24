//! # Traits
//!
//! A hierarchy of number types is defined. The hierarchy is "mathematically exact", but the
//! implementations aren't. That is, the contracts that these traits define, or their names imply,
//! may not be kept precisely. This is due to finite representation of these numbers and is a
//! fundamental problem that cannot be avoided, but perhaps be dealt with differently.
use std::fmt::{Debug, Display};
use std::iter::Sum;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use num::{One, Zero};
use crate::data::number_types::nonzero::NonzeroSigned;

pub mod factorization;

/// The simplex algorithm is defined over the ordered fields. All methods containing algorithmic
/// logic should be defined to work an ordered field (or a field, if they don't need the ordering).
/// All methods representing a matrix should be defined over a field, because they don't need the
/// additional ordering.
pub trait OrderedField =
    Ord +
    NonzeroSigned +
    Field +
    Sized +
;

/// Absolute value of a number.
///
/// Automatically implemented for all types satisfying the trait's bounds.
pub trait Abs: Neg<Output=Self> + Ord + Zero {
    /// The absolute value of a number.
    ///
    /// Compute the additive inverse if the number is smaller than the additive identity.
    fn abs(self) -> Self {
        if self < Self::zero() {
            -self
        } else {
            self
        }
    }
}
impl<T: Neg<Output=Self> + Ord + Zero> Abs for T {
}

/// A reference to an ordered field.
pub trait OrderedFieldRef<Deref> = Ord + FieldRef<Deref>;

/// Basic field operations with Self and with references to Self.
pub trait Field =
    PartialEq + // Equivalence relation
    Eq +
    PartialOrd +
    Zero + // Additive identity
    Neg<Output=Self> + // Additive inverse
    One + // Multiplicative identity
    // First operation
    Add<Self, Output=Self> +
    for<'r> Add<&'r Self, Output=Self> +
    AddAssign<Self> +
    for<'r> AddAssign<&'r Self> +
    Sum +
    // First operation inverse
    Sub<Self, Output=Self> +
    for<'r> Sub<&'r Self, Output=Self> +
    SubAssign<Self> +
    for<'r> SubAssign<&'r Self> +
    // Second operation
    Mul<Self, Output=Self> +
    for<'r> Mul<&'r Self, Output=Self> +
    MulAssign<Self> +
    for<'r> MulAssign<&'r Self> +
    // Second operation inverse
    Div<Self, Output=Self> +
    for<'r> Div<&'r Self, Output=Self> +
    DivAssign<Self> +
    for<'r> DivAssign<&'r Self> +
    // TODO: MulAdd should be possible. Only in specialization?
    //  + MulAdd

    // Practicalities
    Clone +
    Display +
    ToString +
    Debug +
;

/// A reference to a variable that is in a `Field`.
///
/// TODO: Can less HRTB be used? Can the be written down less often? Can this trait be integrated
///  with the `Field` trait?
pub trait FieldRef<Deref> =
    // Equivalence relation
    PartialEq<Self> +
    Neg<Output=Deref> +  // Additive inverse
    // First operation
    Add<Deref, Output=Deref> +
    Add<Output=Deref> +
    // First operation inverse
    Sub<Deref, Output=Deref> +
    Sub<Output=Deref> +
    // Second operation
    Mul<Deref, Output=Deref> +
    Mul<Output=Deref> +
    // Second operation inverse
    Div<Deref, Output=Deref> +
    Div<Output=Deref> +
    // TODO: MulAdd should be possible. Only in specialization?
    //  + MulAdd

    // Practicalities
    Copy +
    Clone +
    Display +
    Debug +
    // Necessary for the Add, Sub, Mul and Div traits. References are sized anyways.
    Sized +
;

/// Helper macro for tests.
#[macro_export]
macro_rules! F {
    ($value:expr) => {
        {
            F::from_f64($value as f64).unwrap()
        }
    };
}

#[cfg(test)]
mod test {
    use std::fmt::Debug;

    use num::One;

    use crate::data::number_types::rational::{Rational128, Rational32, Rational64, RationalBig};
    use crate::data::number_types::traits::Abs;

    /// Testing the `RealField::abs` method.
    #[test]
    fn real_field_abs() {
        fn test<F: Abs + One + Clone + Debug> () {
            let v = F::one();
            let w = -v.clone();
            assert_eq!(w.abs(), v);

            let v = F::one();
            assert_eq!(v.clone().abs() ,v);

            let v = F::zero();
            assert_eq!(v.clone().abs(), v)
        }

        test::<Rational32>();
        test::<Rational64>();
        test::<Rational128>();
        test::<RationalBig>();
    }
}
