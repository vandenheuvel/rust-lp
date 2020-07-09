use smallvec::SmallVec;

use crate::data::number_types::rational::big::integer::{UnsignedDigit, NonZeroSignedDigit};

#[derive(Eq, PartialEq, Clone, Debug)]
pub struct NonZeroSigned {
    lowest: NonZeroSignedDigit,
    data: SmallVec<[UnsignedDigit; 4 - 1]>,
}
