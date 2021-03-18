//! # Representing linear programs for rapid read access
//!
//! The Simplex method algorithms work on a tableau. Because this tableau is very sparse in
//! practice, we store in a matrix that describes the current basis together with the original
//! (also sparse) matrix data. This module contains structures that can provide a matrix.
use crate::algorithm::two_phase::matrix_provider::column::{Column, OrderedColumn};
use crate::data::linear_algebra::traits::{Element, SparseComparator, SparseElement};
use crate::data::linear_algebra::vector::{DenseVector, SparseVector};
use crate::data::linear_program::elements::BoundDirection;

pub mod column;
pub mod matrix_data;
pub mod filter;
pub mod variable;

/// Abstract interface for a matrix and constraint vector.
///
/// This is the data of the "problem relative to the initial basis"; that is, nothing in
/// data structures implementing this trait determines a basis. The implementors of this trait
/// should be primarily read-only, with basis changes, the `Carry` fields of the `Tableau`
/// should change instead.
///
/// Note that a this trait doesn't have to be implemented by a (sparse) matrix data structure per
/// se; it could also be implemented by a graph, which lets itself be represented by data in a
/// matrix.
/// The indexing for the variables and constraints is as follows:
///
/// /                 || Vars of which we want a solution | Constraint slack vars | Bound slack vars |
/// ==================||==================================|=======================|==================|-----
/// Constraints       ||            constants             |       constants       |         0        || b |
/// ------------------||----------------------------------|-----------------------|------------------||---|
///                   ||                                  |                       | +/- 1            |
/// Bound constraints ||    constants (one 1 per row)     |           0           |       +/- 1      |
///                   ||                                  |                       |            +/- 1 |
/// --------------------------------------------------------------------------------------------------
pub trait MatrixProvider {
    /// Representation of a column of the matrix.
    ///
    /// TODO(ARCHITECTURE): When GATs are working, cloning can be avoided in some implementations,
    ///  such as the ones that explicitly store the column data, by giving this associated type a
    ///  lifetime parameter. Keep an eye on <https://github.com/rust-lang/rust/issues/44265>.
    /// TODO(ARCHITECTURE): When specializing on the generic arguments of trait methods is possible,
    ///  the columns no longer need to be ordered necessarily and the bound can be removed here.
    type Column<'a>: Column + OrderedColumn;
    /// Cost row type.
    ///
    /// This type will often be of the form `Option<_>` so to not have to store any zero values, the
    /// inner type would never be zero in that case.
    type Cost<'a>;

    /// Right hand side type.
    type Rhs: Element;

    /// Column of the problem.
    ///
    /// # Arguments
    ///
    /// * `j`: Column index.
    ///
    /// # Return value
    ///
    /// A sparse vector.
    fn column(&self, j: usize) -> Self::Column<'_>;

    /// Cost of a variable.
    ///
    /// # Arguments
    ///
    /// * `j`: Column index.
    ///
    /// # Return value
    ///
    /// Cost value.
    fn cost_value(&self, j: usize) -> Self::Cost<'_>;

    /// Constraint values.
    ///
    /// Note: constraint values of both the constraints and bounds. Lengths should be
    /// `self.nr_rows()`.
    ///
    /// TODO(OPTIMIZATION): Can this clone be avoided?
    ///
    /// # Return value
    ///
    /// A dense vector of constraint values, often called `b` in mathematical notation.
    fn right_hand_side(&self) -> DenseVector<Self::Rhs>;

    /// Index of the row of a virtual bound, if any.
    ///
    /// TODO(ARCHITECTURE): Currently, the return value is a row index. Make this relative to
    ///  `self.nr_constraints`?
    ///
    /// # Arguments
    ///
    /// * `j`: Variable index for the bound, if it exists.
    /// * `bound_type`: Whether it concerns a lower or upper bound.
    ///
    /// # Return value
    ///
    /// The index of the row in which the bound is virtually represented, if the bound exists.
    fn bound_row_index(&self, j: usize, bound_type: BoundDirection) -> Option<usize>;

    /// The number of constraints in the problem. This excludes simple variable bounds.
    fn nr_constraints(&self) -> usize;

    /// The number of simple variable bounds in the problem. This excludes more complicated
    /// constraints.
    fn nr_variable_bounds(&self) -> usize;

    /// The total number of rows in the provided virtual matrix.
    fn nr_rows(&self) -> usize {
        self.nr_constraints() + self.nr_variable_bounds()
    }

    /// The total number of columns in the provided virtual matrix. This does not include artificial
    /// variables; those are virtually represented by the `Artificial` `TableauType`.
    fn nr_columns(&self) -> usize;

    /// Reconstruct a solution.
    ///
    /// Not all variables that a provider presents to the solution algorithms might be relevant for
    /// the final solution. Free variables that are split, for example, could be recombined here.
    ///
    /// # Arguments
    ///
    /// * `column_values`: A solution for each of the variables that this provider presents.
    ///
    /// # Return value
    ///
    /// A solution that might be smaller than the number of variables in this problem.
    fn reconstruct_solution<G>(&self, column_values: SparseVector<G, G>) -> SparseVector<G, G>
    where
        G: SparseElement<G> + SparseComparator,
    ;
}
