use std::cmp::Ordering;
use std::num::{NonZeroU64, NonZeroI64};

pub mod signed;
pub mod unsigned;
pub mod non_zero_signed;
pub mod non_zero_unsigned;

pub type UnsignedDigit = u64;
pub type SignedDigit = i64;
pub type NonZeroUnsignedDigit = NonZeroU64;
pub type NonZeroSignedDigit = NonZeroI64;

pub(crate) fn cmp_slice(a: &[UnsignedDigit], b: &[UnsignedDigit]) -> Ordering {
    debug_assert!(a.last() != Some(&0));
    debug_assert!(b.last() != Some(&0));

    match Ord::cmp(&a.len(), &b.len()) {
        Ordering::Equal => Iterator::cmp(a.iter().rev(), b.iter().rev()),
        other => other,
    }
}
