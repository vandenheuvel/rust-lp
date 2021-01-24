//! # Building blocks to describe linear programs.
use std::fmt;
use std::ops::Mul;
use std::ops::Not;

use crate::data::linear_program::solution::Solution;
use crate::data::number_types::nonzero::NonzeroSign;

/// When a constraint can not be an equality constraint.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum InequalityRelation {
    /// <a, x> <= b
    Less,
    /// <a, x> >= b
    Greater,
}

impl From<InequalityRelation> for ConstraintRelation {
    fn from(relation: InequalityRelation) -> Self {
        match relation {
            InequalityRelation::Less => Self::Less,
            InequalityRelation::Greater => Self::Greater
        }
    }
}

/// A constraint type describes the bound an equation implies.
///
/// Can be equality, but not a range.
#[allow(missing_docs)]
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ConstraintRelation {
    /// <a, x> <= b
    Less,
    /// <a, x> == b
    Equal,
    /// <a, x> >= b
    Greater,
}

/// Constraint relation which can be a range.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum RangedConstraintRelation<F> {
    /// <a, x> == b
    Equal,
    /// Given as the size of the range.
    ///
    /// The constraint is of the form `l <= <a, x> <= u`. This is typically represented with
    /// `<a, x> + s == b` and `0 <= s <= u_s`, with `b = u` and `u_s = u - l =: r`. We have
    /// `r >= 0`. The `r` is stored in this variant.
    ///
    /// TODO(CORRECTNESS): Should the value stored be strictly positive?
    Range(F),
    /// <a, x> <= b
    Less,
    /// <a, x> >= b
    Greater,
}

impl<F> From<InequalityRelation> for RangedConstraintRelation<F> {
    fn from(relation: InequalityRelation) -> Self {
        match relation {
            InequalityRelation::Less => Self::Less,
            InequalityRelation::Greater => Self::Greater,
        }
    }
}

impl<F> From<ConstraintRelation> for RangedConstraintRelation<F> {
    fn from(relation: ConstraintRelation) -> Self {
        match relation {
            ConstraintRelation::Less => Self::Less,
            ConstraintRelation::Equal => Self::Equal,
            ConstraintRelation::Greater => Self::Greater,
        }
    }
}

impl<F: fmt::Display> fmt::Display for RangedConstraintRelation<F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            RangedConstraintRelation::Less => f.write_str("<="),
            RangedConstraintRelation::Equal => f.write_str("=="),
            RangedConstraintRelation::Greater => f.write_str(">="),
            RangedConstraintRelation::Range(r) => write!(f, "={}=", r),
        }
    }
}

/// Constraint relation which can be a range.
///
/// This is essentially the `RangedConstraintRelation` enum, but without the range value stored such
/// that it can be cheaply copied.
pub enum RangedConstraintRelationKind {
    /// <a, x> == b
    Equal,
    /// b_l <= <a, x> <= b_u
    Range,
    /// <a, x> <= b
    Less,
    /// <a, x> >= b
    Greater,
}

impl<F> From<&RangedConstraintRelation<F>> for RangedConstraintRelationKind {
    fn from(constraint: &RangedConstraintRelation<F>) -> Self {
        match constraint {
            RangedConstraintRelation::Equal => Self::Equal,
            RangedConstraintRelation::Range(_) => Self::Range,
            RangedConstraintRelation::Less => Self::Less,
            RangedConstraintRelation::Greater => Self::Greater,
        }
    }
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

impl Not for BoundDirection {
    type Output = Self;

    fn not(self) -> Self::Output {
        match self {
            Self::Lower => Self::Upper,
            Self::Upper => Self::Lower,
        }
    }
}

/// This is like multiplying an equation with -1.
///
/// If the sign of a coefficient is negative, you need to often flip a bound.
///
/// Example:
///
/// * `3 x >= 2 <=> x >= 2 / 3`
/// * `-3 x >= -2 <=> x <= 2 / 3`
impl Mul<NonzeroSign> for BoundDirection {
    type Output = Self;

    fn mul(self, other: NonzeroSign) -> Self::Output {
        match (self, other) {
            (Self::Upper, NonzeroSign::Positive) | (Self::Lower, NonzeroSign::Negative) => Self::Upper,
            (Self::Lower, NonzeroSign::Positive) | (Self::Upper, NonzeroSign::Negative) => Self::Lower,
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
