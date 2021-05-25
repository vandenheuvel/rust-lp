//! # Importing MPS files
//!
//! Reading of `.mps` files, or files of the Mathematical Programming System format.
//!
//! See <http://lpsolve.sourceforge.net/5.5/mps-format.htm> for a specification.
use std::fmt::{Display, Formatter};
use std::fmt;

use relp_num::Rational64;

use crate::data::linear_algebra::{SparseTuple, SparseTupleVec};
use crate::data::linear_program::elements::{ConstraintRelation, Objective};
use crate::data::linear_program::elements::VariableType;
use crate::io::error::Import;
use crate::io::mps::parse::{fixed, free};

#[allow(clippy::type_complexity)]
mod convert;
pub mod number;
pub mod parse;
mod token;

/// Parse an MPS program, in string form, to a MPS.
///
/// # Arguments
///
/// * `program`: The input in [MPS format](https://en.wikipedia.org/wiki/MPS_(format)).
///
/// # Return value
///
/// A `Result<MPS, ImportError>` instance.
///
/// # Errors
///
/// An Import error, wrapping either a parse error indicating that the file was syntactically
/// incorrect, or an Inconsistency error indicating that the file is "logically" incorrect.
pub fn parse(
    program: &impl AsRef<str>,
) -> Result<MPS<Rational64>, Import> {
    free::parse(program.as_ref())
}

/// Parse an MPS program, in string form, to a MPS with struct assumptions on file layout.
///
/// # Arguments
///
/// * `program`: The input in [MPS format](https://en.wikipedia.org/wiki/MPS_(format)).
///
/// # Return value
///
/// A `Result<MPS, ImportError>` instance.
///
/// # Errors
///
/// An Import error, wrapping either a parse error indicating that the file was syntactically
/// incorrect, or an Inconsistency error indicating that the file is "logically" incorrect.
pub fn parse_fixed(
    program: &impl AsRef<str>,
) -> Result<MPS<Rational64>, Import> {
    fixed::parse(program.as_ref())
}

/// Represents the contents of a MPS file in a structured manner.
///
/// `usize` variables in contained structs refer to the index of the cost and row names.
#[derive(Debug, PartialEq)]
pub struct MPS<F> {
    /// Name of the linear program.
    name: String,
    /// Whether this is a minimization or maximization.
    objective: Objective,

    /// Name of the cost row.
    cost_row_name: String,
    /// Variable index and value tuples, describing how the variables appear in the objective
    /// function.
    ///
    /// Column (by index) and coefficient combinations for the objective function.
    cost_values: SparseTupleVec<F>,

    /// All named constraint types (see the ConstraintType enum).
    ///
    /// Ordering corresponds to the row_names field.
    rows: Vec<Row>,
    /// Constraint name and variable name combinations.
    ///
    /// Ordering in each variable corresponds to the row_names field.
    columns: Vec<Column<F>>,
    /// Right-hand side constraint values.
    ///
    /// Ordering in each right hand side corresponds to the row_names field.
    rhss: Vec<Rhs<F>>,
    /// Limiting constraint activations two-sidedly.
    ranges: Vec<Range<F>>,
    /// Bounds on variables.
    bounds: Vec<Bound<F>>,
}

#[allow(clippy::too_many_arguments)]
impl<F> MPS<F> {
    /// Simple constructor without any logic.
    #[must_use]
    pub fn new(
        name: String,
        objective: Objective,
        cost_row_name: String,
        cost_values: SparseTupleVec<F>,
        rows: Vec<Row>,
        columns: Vec<Column<F>>,
        rhss: Vec<Rhs<F>>,
        ranges: Vec<Range<F>>,
        bounds: Vec<Bound<F>>,
    ) -> Self {
        Self {
            name,
            objective,
            cost_row_name,
            cost_values,
            rows,
            columns,
            rhss,
            ranges,
            bounds,
        }
    }
}

/// MPS files are divided into sections.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum Section {
    /// After the name, the first section.
    Rows,
    /// Follows `Rows` section.
    Columns,
    /// Right-hand-sides, constraints, i.e. the b of Ax >=< b.
    Rhs,
    /// Bounds on variables.
    Bounds,
    /// Bounds on rows.
    ///
    /// Note: This section is not used in the MIPLIB 2010 benchmark.
    Ranges,
    /// The `Endata` variant (notice the odd spelling) denotes the end of the file.
    Endata,
}

impl Display for Section {
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            Section::Rows => "ROWS",
            Section::Columns => "COLUMNS",
            Section::Rhs => "RHS",
            Section::Bounds => "BOUNDS",
            Section::Ranges => "RANGES",
            Section::Endata => "ENDATA",
        })
    }
}

/// Every row is either a cost row or some constraint.
#[derive(Debug, Eq, PartialEq)]
enum RowType {
    /// Typically, there is only one cost row per problem. Otherwise, it's a multi-objective problem
    /// and it's not clear how they should be parsed.
    Cost,
    /// Other rows are constraints.
    Constraint(ConstraintRelation),
}

/// The MPS format defines the `BoundType`s described in this enum.
///
/// # Note
///
/// Not all `BoundType` variants are currently supported.
#[derive(Copy, Clone, Debug, Eq, PartialEq)]
pub enum BoundType<F> {
    /// b <- x (< +inf)
    LowerContinuous(F),
    /// (0 <=) x <= b
    UpperContinuous(F),
    /// x = b
    Fixed(F),
    /// -inf < x < +inf
    Free,
    /// -inf < x (<= 0)
    LowerMinusInfinity,
    /// (0 <=) x < +inf
    UpperInfinity,
    /// x = 0 or 1
    Binary,
    /// b <= x (< +inf)
    LowerInteger(F),
    /// (0 <=) x <= b
    UpperInteger(F),
    /// x = 0 or l =< x =< b
    ///
    /// Note: appears only very rarely in the MIPLIB benchmark set.
    SemiContinuous(F, F),
}

/// Every `Row` has a name and a `RowType`.
#[derive(Debug, Eq, PartialEq, Clone)]
pub struct Row {
    /// Each row has a name. Is not really used after parsing, but is stored for writing.
    pub name: String,
    /// Direction of the constraint.
    pub constraint_type: ConstraintRelation,
}

/// Is either continuous or integer, and has for some rows a coefficient.
///
/// Note that all values should be non zero, as this is a sparse representation.
#[derive(Debug, PartialEq)]
pub struct Column<F> {
    /// Name of the variable. Used for outputting the solution.
    pub name: String,
    /// Whether the variable is integer or not.
    pub variable_type: VariableType,
    /// Sparse representation of non-zero values in the column.
    pub values: Vec<SparseTuple<F>>,
}

/// The right-hand side of Ax = b.
///
/// A single linear program defined in MPS can have multiple right-hand sides. It relates a row name
/// to a real constant.
#[derive(Debug, PartialEq)]
pub struct Rhs<F> {
    /// Name of the right-hand side. Stored only for writing the problem to disk.
    pub name: String,
    /// Sparse representation of the constraint values.
    ///
    /// Note that there could be multiple per index, although this is rather rare.
    pub values: Vec<SparseTuple<F>>,
}

/// Specifies a bound on constraint activation.
///
/// Overview of how the range of a constraint is defined, depending on the constraint type:
///
/// row type | sign of r |    h    |    u
/// ---------|-----------|---------|---------
/// G        |  + or -   |    b    | b + |r|
/// L        |  + or -   | b - |r| |   b
/// E        |     +     |    b    | b + |r|
/// E        |     -     | b - |r| |   b
#[derive(Debug, PartialEq)]
pub struct Range<F> {
    /// Name of the range. Stored only for writing the problem to disk.
    pub name: String,
    /// Sorted constraint indices and their 'r' value.
    ///
    /// Typically, only few rows would have a range.
    pub values: Vec<SparseTuple<F>>,
}

/// Specifies a bound on a variable. The variable can either be continuous or integer, while the
/// bound can have any direction.
#[derive(Debug, PartialEq)]
pub struct Bound<F> {
    /// Name of the bound. Stored only for writing the problem to disk.
    pub name: String,
    /// Collection of at most one bound per row. Note that across different bounds, multiple values
    /// per row may be specified.
    pub values: Vec<SparseTuple<BoundType<F>>>,
}
