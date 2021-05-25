//! # Properties of number types
use relp_num::NonZero;

use crate::data::linear_algebra::traits::Element;

pub trait Rhs =
    NonZero +
    Eq +
    Element +
;
