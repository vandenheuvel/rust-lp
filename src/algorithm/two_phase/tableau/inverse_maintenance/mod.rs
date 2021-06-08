//! # Maintaining a basis inverse
//!
//! The simplex method requires us to keep track of the basis inverse. The latest constraint values
//! and a variable called `minus_pi` are also tracked, as is the objective value.
//!
//! It is currently only possible to track this inverse using a sparse, row-major matrix (no
//! factorization).
//!
//! TODO(ENHANCEMENT): A better inverse maintenance algorithm. Start with factorization?
use std::fmt::{Debug, Display};

use crate::algorithm::two_phase::matrix_provider::column::{Column, ColumnNumber, OrderedColumn};
use crate::algorithm::two_phase::matrix_provider::filter::Filtered;
use crate::algorithm::two_phase::matrix_provider::MatrixProvider;
use crate::algorithm::two_phase::tableau::kind::Kind;
use crate::data::linear_algebra::SparseTuple;
use crate::data::linear_algebra::traits::Element;
use crate::data::linear_algebra::vector::{DenseVector, SparseVector};

pub mod carry;
pub mod ops;

/// Maintain a basis inverse.
///
/// Facilitates quick solving of a linear system.
///
/// TODO(ARCHITECTURE): Consider making this trait generic over `MatrixProvider`'s associated types.
///  See issue #14.
pub trait InverseMaintener: Display + Sized {
    /// Type used for computations inside the instance.
    ///
    /// Because the algorithm works with results from this object, many other parts of the code use
    /// the same number type. Examples are the tableau and pivot rules.
    type F: ops::Field;

    /// Contains the computed column and potentially other information that can be reused from that
    /// process.
    type ColumnComputationInfo: ColumnComputationInfo<Self::F>;

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
    fn create_for_fully_artificial<Rhs: Element>(rhs: DenseVector<Rhs>) -> Self
    where
        Self::F: ops::Rhs<Rhs>,
    ;

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
    fn create_for_partially_artificial<G: Element>(
        artificial: &[usize],
        basis: &[(usize, usize)],
        b: DenseVector<G>,
        basis_indices: Vec<usize>,
    ) -> Self
    where
        Self::F: ops::Rhs<G>,
    ;

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
    fn from_basis<'a, MP: MatrixProvider>(basis: &[usize], provider: &'a MP) -> Self
    where
        Self::F:
            ops::Column<<MP::Column as Column>::F> +
            ops::Rhs<MP::Rhs> +
            ops::Cost<MP::Cost<'a>> +
            ops::Column<MP::Rhs> +
        ,
        MP::Rhs: 'static,
    ;

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
    fn from_basis_pivots<'a, MP: MatrixProvider>(basis: &[(usize, usize)], provider: &'a MP) -> Self
    where
        Self::F:
            ops::Column<<MP::Column as Column>::F> +
            ops::Rhs<MP::Rhs> +
            ops::Cost<MP::Cost<'a>> +
            ops::Column<MP::Rhs> +
        ,
        MP::Rhs: 'static + ColumnNumber,
    ;

    /// When a previous basis inverse representation was used to find a basic feasible solution.
    ///
    /// Artificial variables need to be removed somewhere in the process that calls this method.
    ///
    /// # Arguments
    ///
    /// * `artificial`: Indices of rows where an artificial variable is needed.
    /// * `provider`: Original problem representation.
    /// * `basis`: (row index, column index) tuples of given basis variables.
    // TODO(ARCHITECTURE): Specialize this method with the "remove_rows" version below.
    fn from_artificial<'provider, MP: MatrixProvider>(
        artificial: Self,
        provider: &'provider MP,
        nr_artificial: usize,
    ) -> Self
    where
        Self::F: ops::FieldHR + ops::Column<<MP::Column as Column>::F> + ops::Cost<MP::Cost<'provider>>,
    ;

    /// When a previous basis inverse representation was used to find a basic feasible solution.
    ///
    /// Artificial variables need to be removed somewhere in the process that calls this method.
    ///
    /// # Arguments
    ///
    /// * `artificial`: Indices of rows where an artificial variable is needed.
    /// * `provider`: Original problem representation.
    /// * `basis`: (row index, column index) tuples of given basis variables.
    fn from_artificial_remove_rows<'a, MP: Filtered>(
        artificial: Self,
        rows_removed: &'a MP,
        nr_artificial: usize,
    ) -> Self
    where
        Self::F: ops::Column<<<MP as MatrixProvider>::Column as Column>::F> + ops::Cost<MP::Cost<'a>>,
    ;

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
    ///
    /// # Return value
    ///
    /// The index of the column being removed from the basis.
    fn change_basis(
        &mut self,
        pivot_row_index: usize,
        pivot_column_index: usize,
        column: Self::ColumnComputationInfo,
        cost: Self::F,
    ) -> usize;

    /// Calculates the cost difference `c_j`.
    ///
    /// This cost difference is the inner product of `minus_pi` and the column.
    // TODO(ENHANCEMENT): Drop the OrderedColumn trait bound once it is possible to specialize on
    //  it.
    fn cost_difference<G, C: Column<F=G> + OrderedColumn>(&self, original_column: &C) -> Self::F
    where
        Self::F: ops::Column<G>,
        G: Display + Debug,
    ;

    /// Multiplies the submatrix consisting of `minus_pi` and B^-1 by an `original_column`.
    ///
    /// # Arguments
    ///
    /// * `original_column`: A `SparseVector<T>` of length `m`.
    ///
    /// # Return value
    ///
    /// A `SparseVector<T>` of length `m`.
    /// TODO(ENHANCEMENT): Drop the `OrderedColumn` trait bound once it is possible to specialize on
    ///  it.
    fn generate_column<G, C: Column<F=G> + OrderedColumn>(
        &self,
        original_column: C,
    ) -> Self::ColumnComputationInfo
    where
        Self::F: ops::Column<G>,
    ;

    /// Generate a single element in the tableau with respect to the current basis.
    ///
    /// # Arguments
    ///
    /// * `i`: Row index
    /// * `original_column`: Column with respect to the original basis.
    /// TODO(ENHANCEMENT): Drop the `OrderedColumn` trait bound once it is possible to specialize on
    ///  it.
    fn generate_element<C: Column + OrderedColumn>(
        &self,
        i: usize,
        original_column: C,
    ) -> Option<Self::F>
    where
        Self::F: ops::Column<C::F>,
    ;

    /// Housekeeping that needs to happen after a basis change.
    ///
    /// # Arguments
    ///
    /// * `kind`: Provides access to the matrix being solved, can be used to retrieve columns.
    fn after_basis_change<K: Kind>(&mut self, kind: &K)
    where
        Self::F: ops::Column<<<K as Kind>::Column as Column>::F>,
    ;

    /// Extract the current basic feasible solution.
    fn current_bfs(&self) -> Vec<SparseTuple<Self::F>>;

    /// Which basis column has a pivot in the provided row index?
    fn basis_column_index_for_row(&self, row: usize) -> usize;

    /// Clone the latest constraint vector.
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

/// Allows additional values to be computed when the column is computed.
pub trait ColumnComputationInfo<F>: Debug {
    /// Get a reference to the column.
    fn column(&self) -> &SparseVector<F, F>;
    /// Consume this value, leaving only the inner column.
    fn into_column(self) -> SparseVector<F, F>;
}
