//! # Building blocks to describe linear programs.
use std::ops::{BitXor, Neg};
use std::ops::Not;

use num::One;

use crate::data::linear_program::solution::Solution;

/// A `Constraint` is a type of (in)equality.
#[allow(missing_docs)]
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ConstraintType {
    Equal,
    Greater,
    Less,
}

/// Direction of a bound.
///
/// Is used more generally in the case where the three variants of the `ConstraintType` don't suit
/// the needs.
#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
pub enum BoundDirection {
    /// In the case of a variable, x >= b >= 0 (currently, variables are rewritten to be nonnegative
    /// and branching might only bring that bound higher).
    Lower,
    /// In the case of a variable, 0 <= x <= b.
    Upper,
}

impl BoundDirection {
    /// Convert a bound direction into a positive or negative one.
    /// 
    /// Upper bounds needs positive slacks, lower bounds need negative slacks.
    #[must_use]
    pub fn into<F: One + Neg<Output = F>>(self) -> F {
        match self {
            BoundDirection::Lower => -F::one(),
            BoundDirection::Upper => F::one(),
        }
    }
}

impl Not for BoundDirection {
    type Output = Self;

    fn not(self) -> Self::Output {
        match self {
            Self::Lower => Self::Upper,
            Self::Upper => Self::Lower,
        }
    }
}

/// Analogue to multiplying signs of values.
///
/// Used mostly in presolving.
impl BitXor for BoundDirection {
    type Output = Self;

    fn bitxor(self, other: Self) -> Self::Output {
        match (self, other) {
            (Self::Lower, Self::Upper) => Self::Upper,
            (Self::Upper, Self::Lower) => Self::Upper,
            (Self::Lower, Self::Lower) => Self::Lower,
            (Self::Upper, Self::Upper) => Self::Lower,
        }
    }
}

/// A variable is either continuous or integer.
#[allow(missing_docs)]
#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub enum VariableType {
    Continuous,
    Integer,
}

impl Not for VariableType {
    type Output = VariableType;

    fn not(self) -> VariableType {
        match self {
            VariableType::Continuous => VariableType::Integer,
            VariableType::Integer => VariableType::Continuous,
        }
    }
}

/// After the second phase, either an optimum is found or the problem is determined to be unbounded.
#[allow(missing_docs)]
#[derive(Debug, Eq, PartialEq)]
pub enum LinearProgramType<F> {
    FiniteOptimum(Solution<F>),
    Infeasible,
    Unbounded,
}

/// Direction of optimization.
#[allow(missing_docs)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Objective {
    Maximize,
    Minimize,
}
impl Default for Objective {
    fn default() -> Self {
        Objective::Minimize
    }
}
