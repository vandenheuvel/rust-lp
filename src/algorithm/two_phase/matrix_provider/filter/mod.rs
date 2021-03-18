//! # Removing rows from a matrix provider
//!
//! If a problem turns out to have redundant constraints after the first phase, we need to remove
//! those rows from the problem for the rest of the implementation to work.
use crate::algorithm::two_phase::matrix_provider::column::Column;
use crate::algorithm::two_phase::matrix_provider::MatrixProvider;

pub mod generic_wrapper;

/// A filtered matrix provider had rows removed from it.
pub trait Filtered: MatrixProvider {
    /// The rows that were removed.
    ///
    /// Indexes are relevant to the original problem.
    fn filtered_rows(&self) -> &[usize];
}

/// Derive a variant of the matrix provider that has rows removed from it.
pub trait ToFiltered: MatrixProvider {
    /// The resulting matrix provider type.
    type Filtered<'provider>: Filtered<Column: Column<F=<Self::Column<'provider> as Column>::F>>;

    /// Derive a variant of the matrix provider that has rows removed from it.
    ///
    /// # Arguments
    ///
    /// * `rows_to_filter`: Indices of rows to remove from the problem. Indices are relative to the
    /// original problem.
    fn to_filtered(&self, rows_to_filter: Vec<usize>) -> Self::Filtered<'_>;
}
