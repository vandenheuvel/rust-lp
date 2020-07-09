use std::fmt::{Debug, Display};

pub trait Field =
    Ord +
    Debug +
    Display +
    Clone +
;

pub trait NonZero =
    Ord +
    Debug +
    Display +
    Clone +
;
