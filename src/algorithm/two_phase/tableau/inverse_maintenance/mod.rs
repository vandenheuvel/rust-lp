//! # Maintaining a basis inverse
//!
//! The simplex method requires us to keep track of the basis inverse. The latest constraint values
//! and a variable called `minus_pi` are also tracked, as is the objective value.
//!
//! It is currently only possible to track this inverse using a sparse, row-major matrix (no
//! factorization).
//!
//! TODO(ENHANCEMENT): A better inverse maintenance algorithm. Start with factorization?
use std::fmt::Display;

use crate::algorithm::two_phase::matrix_provider::{Column, MatrixProvider, OrderedColumn};
use crate::algorithm::two_phase::matrix_provider::filter::Filtered;
use crate::data::linear_algebra::SparseTuple;
use crate::data::linear_algebra::vector::{Dense as DenseVector, Sparse as SparseVector};
use crate::data::number_types::traits::Field;

pub mod carry;

/// Maintain a basis inverse.
///
/// Should facilitate quick solving of a linear system.
pub trait InverseMaintenance: Display {
    /// Type used for computations inside the instance.
    ///
    /// Because the algorithm works with results from this object, many other parts of the code use
    /// the same number type. Examples are the tableau and pivot rules.
    type F: 'static + Field;

    /// Create a `Carry` for a tableau with only artificial variables.
    ///
    /// Only used if it is known in advance that there will be no positive slacks.
    ///
    /// # Arguments
    ///
    /// * `b`: Constraint values of the original problem, i.e. the original `b` with respect to the
    /// unit basis.
    ///
    /// # Return value
    ///
    /// `Carry` with a `minus_pi` equal to -1's and the standard basis.
    fn create_for_fully_artificial(b: DenseVector<Self::F>) -> Self;

    /// Create a `Carry` for a tableau with some artificial variables.
    ///
    /// Only used if unknown in advance how many positive slacks there will be.
    ///
    /// # Arguments
    ///
    /// * `artificial`: Indices of rows where an artificial variable is needed.
    /// * `basis`: (row index, column index) tuples of given basis variables.
    /// * `b`: Constraint values of the original problem, i.e. the original `b` with respect to the
    /// unit basis.
    ///
    /// # Return value
    ///
    /// `Carry` with a `minus_pi` equal to -1's and the standard basis.
    fn create_for_partially_artificial(
        artificial: &[usize],
        basis: &Vec<(usize, usize)>,
        b: DenseVector<Self::F>,
    ) -> Self;

    /// Create a basis inverse when only the basis indices are known.
    ///
    /// Often, one would invert a matrix right at the start here to build up the basis inverse
    /// representation.
    ///
    /// # Arguments
    ///
    /// * `basis`: Indices of columns that are to be in the basis. Should match the number of rows
    /// of the provider. Values should be unique, could have been a set.
    /// * `provider`: Problem representation.
    fn from_basis(basis: &[usize], provider: &impl MatrixProvider) -> Self;

    /// Create a basis inverse when the basis indices and their pivot rows are known.
    ///
    /// Often, one would invert a matrix right at the start here to build up the basis inverse
    /// representation.
    ///
    /// # Arguments
    ///
    /// * `basis`: Indices of columns that are to be in the basis. Should match the number of rows
    /// of the provider. Values should be unique, could have been a set.
    /// * `provider`: Problem representation.
    fn from_basis_pivots(
        basis: &Vec<(usize, usize)>,
        provider: &impl MatrixProvider,
    ) -> Self;

    /// When a previous basis inverse representation was used to find a basic feasible solution.
    ///
    /// Artificial variables need to be removed somewhere in the process that calls this method.
    ///
    /// # Arguments
    ///
    /// * `artificial`: Indices of rows where an artificial variable is needed.
    /// * `provider`: Original problem representation.
    /// * `basis`: (row index, column index) tuples of given basis variables.
    fn from_artificial(
        artificial: Self,
        // TODO(ENHANCEMENT): Decouple these two `F` types
        provider: &impl MatrixProvider<Column: Column<F=Self::F>>,
        basis: &[usize],
    ) -> Self;

    /// When a previous basis inverse representation was used to find a basic feasible solution.
    ///
    /// Artificial variables need to be removed somewhere in the process that calls this method.
    ///
    /// # Arguments
    ///
    /// * `artificial`: Indices of rows where an artificial variable is needed.
    /// * `provider`: Original problem representation.
    /// * `basis`: (row index, column index) tuples of given basis variables.
    fn from_artificial_remove_rows(
        artificial: Self,
        rows_removed: &impl Filtered<Column: Column<F=Self::F>>,
        basis: &[usize],
    ) -> Self;

    /// Update the basis by representing one row reduction operation.
    ///
    /// Supply a column, that should become part of the basis, to this function. Then this
    /// `Carry` gets updated such that if one were to matrix-multiply the concatenation column
    /// `[cost, c for c in column]` with this `Carry`, an (m + 1)-dimensional unitvector would
    /// be the result.
    ///
    /// # Arguments
    ///
    /// * `pivot_row_index`: Index of the pivot row.
    /// * `column`: Column relative to the current basis to be entered into that basis.
    /// * `cost`: Relative cost of that column. The objective function value will change by this
    /// amount.
    fn change_basis(&mut self, pivot_row_index: usize, column: &SparseVector<Self::F, Self::F>, cost: Self::F);

    /// Calculates the cost difference `c_j`.
    ///
    /// This cost difference is the inner product of `minus_pi` and the column.
    // TODO(ENHANCEMENT): Drop the OrderedColumn trait bound once it is possible to specialize on
    //  it.
    fn cost_difference<C: Column<F=Self::F> + OrderedColumn>(&self, original_column: &C) -> Self::F;

    /// Multiplies the submatrix consisting of `minus_pi` and B^-1 by a original_column.
    ///
    /// # Arguments
    ///
    /// * `original_column`: A `SparseVector<T>` of length `m`.
    ///
    /// # Return value
    ///
    /// A `SparseVector<T>` of length `m`.
    /// TODO(ENHANCEMENT): Drop the OrderedColumn trait bound once it is possible to specialize on
    ///  it.
    fn generate_column<C: Column<F=Self::F> + OrderedColumn>(&self, original_column: C) -> SparseVector<Self::F, Self::F>;

    /// Generate a single element in the tableau with respect to the current basis.
    ///
    /// # Arguments
    ///
    /// * `i`: Row index
    /// * `original_column`: Column with respect to the original basis.
    /// TODO(ENHANCEMENT): Drop the OrderedColumn trait bound once it is possible to specialize on
    ///  it.
    fn generate_element<'a, I: Iterator<Item=&'a SparseTuple<Self::F>>>(
        &self,
        i: usize,
        original_column: I,
    ) -> Self::F
    where
        Self::F: 'a,
    ;

    /// Clone the latest constraint vector.
    ///
    /// TODO: Can this cloning be avoided?
    fn b(&self) -> DenseVector<Self::F>;

    /// Get the objective function value for the current basis.
    ///
    /// # Return value
    ///
    /// The objective value.
    fn get_objective_function_value(&self) -> Self::F;

    /// Get the `i`th constraint value relative to the current basis.
    ///
    /// # Arguments
    ///
    /// * `i`: Index of the constraint to retrieve value from, in range `0` until `m`.
    ///
    /// # Return value
    ///
    /// The constraint value.
    fn get_constraint_value(&self, i: usize) -> &Self::F;
}
