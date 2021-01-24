//! # Number types
//!
//! Different data structures and different algorithms are defined over different spaces. This
//! module defines those spaces and provides implementations.
//!
//! A benefit of this approach is that the algorithms can be tested well for correctness using
//! fractional numbers.
//! A downside is that the algorithm doesn't have access to the finite representations of e.g.
//! floats, which is potentially limiting.
pub mod traits;

pub mod binary;
pub mod integer;
pub mod nonzero;
pub mod rational;
