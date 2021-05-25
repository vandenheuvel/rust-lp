//! # Number types
//!
//! Aliases for the number types needed when reading and converting MPS structures.
use std::fmt::{Debug, Display};
use std::ops::{Neg, Sub};

use num_traits::{One, Zero};
use relp_num::Abs;

pub trait Field =
    Zero +
    One +

    Neg<Output=Self> +
    Sub<Output=Self> +
    Abs +

    Eq +
    Ord +

    Debug +
    Display +
    Clone +
;
