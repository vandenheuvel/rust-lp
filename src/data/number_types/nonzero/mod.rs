use std::num::{NonZeroI128, NonZeroI16, NonZeroI32, NonZeroI64, NonZeroI8};
use std::num::{NonZeroU128, NonZeroU16, NonZeroU32, NonZeroU64, NonZeroU8};

use num::Zero;

pub use sign::Sign as NonzeroSign;
pub use sign::NonzeroSigned as NonzeroSigned;

pub mod sign;

/// Implementors should not be zero.
///
/// This trait is used for debug asserts. Values in sparse data structures should never be zero, and
/// requiring that they implement `num::Zero` prohibits writing number types that can't represent
/// the value 0.
///
/// The `num::Zero` trait is for types that can be zero, this trait is for types that can be a value
/// other than zero. They may or may not be able to represent zero.
pub trait Nonzero {
    /// Whether the value is not equal to zero.
    ///
    /// Should always be `true`.
    fn is_not_zero(&self) -> bool;
}

macro_rules! could_be_zero {
    ($t: ident) => {
        impl Nonzero for $t {
            fn is_not_zero(&self) -> bool {
                !Zero::is_zero(self)
            }
        }
    }
}

could_be_zero!(i8);
could_be_zero!(u8);
could_be_zero!(i16);
could_be_zero!(u16);
could_be_zero!(i32);
could_be_zero!(u32);
could_be_zero!(i64);
could_be_zero!(u64);
could_be_zero!(i128);
could_be_zero!(u128);

macro_rules! can_not_be_zero {
    ($t: ident) => {
        impl Nonzero for $t {
            fn is_not_zero(&self) -> bool {
                true
            }
        }
    }
}

can_not_be_zero!(NonZeroI8);
can_not_be_zero!(NonZeroI16);
can_not_be_zero!(NonZeroI32);
can_not_be_zero!(NonZeroI64);
can_not_be_zero!(NonZeroI128);
can_not_be_zero!(NonZeroU8);
can_not_be_zero!(NonZeroU16);
can_not_be_zero!(NonZeroU32);
can_not_be_zero!(NonZeroU64);
can_not_be_zero!(NonZeroU128);
