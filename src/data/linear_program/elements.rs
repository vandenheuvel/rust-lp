use std::ops::Not;

/// A `Constraint` is a type of (in)equality.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum ConstraintType {
    Equal,
    Greater,
    Less,
}

/// A variable is either continuous or integer.
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

/// An LP either has a finite optimum, is unbounded or has no basic feasible solution.
#[derive(Debug, Copy, Clone)]
pub enum LPCategory {
    FiniteOptimum(f64),
    Unbounded,
    Infeasible,
}