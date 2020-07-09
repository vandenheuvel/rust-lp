use std::fmt::Display;

/// A non zero number should be able to interact with
/// TODO: Write docs
pub trait NonZero =
where
    for<'r> &'r Self: NonZeroRef,
;

pub trait NonZeroRef =
;

pub trait NonZeroOption =
    Clone +
    Display +
;

pub trait NonNegative =
    Clone +
;
