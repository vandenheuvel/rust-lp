use crate::data::number_types::rational::Ratio;
use std::num::NonZeroU64;

pub type Rational64 = Ratio<i64, NonZeroU64>;
pub type RationalNonNegative64 = Ratio<u64, NonZeroU64>;
