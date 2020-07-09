//! # Importing MPS files
//!
//! Reading of `.mps` files, or files of the Mathematical Programming System format.
//!
//! See http://lpsolve.sourceforge.net/5.5/mps-format.htm for a specification.
//!
//! TODO:
//!     * Support all `BoundType` variants
use std::fmt::{Display, Formatter};
use std::fmt;

use crate::data::linear_algebra::{SparseTuple, SparseTupleVec};
use crate::data::linear_program::elements::ConstraintType;
use crate::data::linear_program::elements::VariableType;
use crate::io::error::Import;
use crate::io::mps::parse::free;
use crate::data::number_types::rational::small::Rational64;

mod convert;
pub mod number;
pub mod parse;
mod parsing;
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
) -> Result<MPS<Rational64, Rational64>, Import> {
    free::parse(program.as_ref())
}

/// Represents the contents of a MPS file in a structured manner.
///
/// `usize` variables in contained structs refer to the index of the cost and row names.
#[derive(Debug, PartialEq)]
pub struct MPS<F, NZ> {
    /// Name of the linear program.
    name: String,

    /// Name of the cost row.
    cost_row_name: String,
    /// Variable index and value tuples, describing how the variables appear in the objective
    /// function.
    ///
    /// Column (by index) and coefficient combinations for the objective function.
    cost_values: SparseTupleVec<NZ>,

    /// All named constraint types (see the ConstraintType enum).
    ///
    /// Ordering corresponds to the row_names field.
    rows: Vec<Row>,
    /// Constraint name and variable name combinations.
    ///
    /// Ordering in each variable corresponds to the row_names field.
    columns: Vec<Column<NZ>>,
    /// Right-hand side constraint values.
    ///
    /// Ordering in each right hand side corresponds to the row_names field.
    rhss: Vec<Rhs<F>>,
    /// Limiting constraint activations two-sidedly.
    ranges: Vec<Range<F>>,
    /// Bounds on variables.
    bounds: Vec<Bound<F>>,
}

impl<F, NZ> MPS<F, NZ> {
    pub fn new(
        name: String,
        cost_row_name: String,
        cost_values: SparseTupleVec<NZ>,
        rows: Vec<Row>,
        columns: Vec<Column<NZ>>,
        rhss: Vec<Rhs<F>>,
        ranges: Vec<Range<F>>,
        bounds: Vec<Bound<F>>,
    ) -> Self {
        Self {
            name,
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
    Rows,
    Columns,
    Rhs,
    Bounds,
    /// This section is not used in the MIPLIB 2010 benchmark.
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
    Cost,
    Constraint(ConstraintType),
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
    pub name: String,
    pub constraint_type: ConstraintType,
}

/// Is either continuous or integer, and has for some rows a coefficient.
///
/// Note that all values should be non zero, as this is a sparse representation.
#[derive(Debug, PartialEq)]
pub struct Column<NZ> {
    pub name: String,
    pub variable_type: VariableType,
    pub values: Vec<SparseTuple<NZ>>,
}

/// The right-hand side of Ax = b.
///
/// A single linear program defined in MPS can have multiple right-hand sides. It relates a row name
/// to a real constant.
#[derive(Debug, PartialEq)]
pub struct Rhs<F> {
    pub name: String,
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
pub struct Range<F> {  // Type parameter naming: can take all values, not just be nonzero
    pub name: String,
    /// Sorted constraint indices and their 'r' value.
    pub values: Vec<SparseTuple<F>>,
}

/// Specifies a bound on a variable. The variable can either be continuous or integer, while the
/// bound can have any direction.
#[derive(Debug, PartialEq)]
pub struct Bound<F> {
    pub name: String,
    pub values: Vec<SparseTuple<BoundType<F>>>,
}
