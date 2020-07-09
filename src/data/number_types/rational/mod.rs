//! # Rational numbers
//!
//! Useful for testing of methods and overall correctness of implementations.
//!
//! It appears that in practice, no fast solvers use rational numbers.
use std::fmt::Debug;

pub mod big;
pub mod small;
pub mod macros;

#[derive(Eq, PartialEq, Debug, Clone)]
struct Ratio<N, D> {
    n: N,
    d: D,
}

impl<N: Clone + Debug, D: Clone + Debug> Ratio<N, D> {
    pub fn new(numerator: N, denominator: D) -> Self {
        Self { n: numerator, d: denominator, }
    }
}
