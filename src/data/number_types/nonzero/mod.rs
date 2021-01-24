use std::num::{NonZeroI128, NonZeroI16, NonZeroI32, NonZeroI64, NonZeroI8};
use std::num::{NonZeroU128, NonZeroU16, NonZeroU32, NonZeroU64, NonZeroU8};

use num::Zero;

pub use sign::Sign as NonzeroSign;
pub use sign::Signed as NonzeroSigned;

pub mod sign;

/// Implementors should not be zero.
///
/// This trait is used for debug asserts. Values in sparse data structures should never be zero, and
/// requiring that they implement `num::Zero` prohibits writing number types that can't represent
/// the value 0.
///
/// The `num::Zero` trait is for types that *can* be zero, this trait is for values that *are* not
/// zero (independent of whether the type could represent the value zero).
pub trait Nonzero {
    /// Whether the value is not equal to zero.
    ///
    /// Should always be `true`.
    fn is_not_zero(&self) -> bool;
}

impl<T: Zero> Nonzero for T {
    default fn is_not_zero(&self) -> bool {
        !self.is_zero()
    }
}

// macro_rules! can_not_be_zero {
//     ($t: ident) => {
//         impl Nonzero for $t {
//             fn is_not_zero(&self) -> bool {
//                 true
//             }
//         }
//     }
// }
// impl Nonzero for NonZeroI8 {
//     fn is_not_zero(&self) -> bool {
//         true
//     }
// }
// can_not_be_zero!(NonZeroI8);
// can_not_be_zero!(NonZeroI16);
// can_not_be_zero!(NonZeroI32);
// can_not_be_zero!(NonZeroI64);
// can_not_be_zero!(NonZeroI128);
// can_not_be_zero!(NonZeroU8);
// can_not_be_zero!(NonZeroU16);
// can_not_be_zero!(NonZeroU32);
// can_not_be_zero!(NonZeroU64);
// can_not_be_zero!(NonZeroU128);
