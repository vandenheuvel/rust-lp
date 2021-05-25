//! When a matrix is represented with an unusual backend, like a network, it might be practical to
//! remove rows from the matrix it represents without having to adapt the underlying implementation.
//! This module provides a wrapper around any matrix provider that removes rows from it.
//!
//! Note that this implementation is there only as a fall-back option; really, implementors of the
//! `MatrixProvider` trait should either ensure that the problem is full rank, or implement one of
//! the traits in the parent module and write code that is efficiently converts their specific
//! instance.
use std::collections::HashSet;
use std::fmt::{Display, Formatter};
use std::fmt;

use crate::algorithm::two_phase::matrix_provider::column::{Column, OrderedColumn};
use crate::algorithm::two_phase::matrix_provider::filter::{Filtered, ToFiltered};
use crate::algorithm::two_phase::matrix_provider::MatrixProvider;
use crate::algorithm::two_phase::matrix_provider::variable::FeasibilityLogic;
use crate::algorithm::two_phase::phase_one::PartialInitialBasis;
use crate::algorithm::utilities::remove_sparse_indices;
use crate::data::linear_algebra::traits::{SparseComparator, SparseElement};
use crate::data::linear_algebra::vector::{DenseVector, SparseVector, Vector};
use crate::data::linear_program::elements::BoundDirection;

/// Remove a set of rows from a column.
///
/// Note that only constraint rows should be removed.
pub trait IntoFilteredColumn: Column {
    /// The type used to represent the filtered version of the column.
    ///
    /// This will often be `Self`.
    ///
    /// This type has no lifetime attached to it, because it is likely better to filter once and
    /// then iterate over references of the filtered object, rather than to filter repeatedly while
    /// iterating.
    type Filtered: Column<F=Self::F> + OrderedColumn;

    /// Filter the column.
    ///
    /// # Arguments
    ///
    /// * `to_remove`: Indices of rows to remove.
    fn into_filtered(self, to_remove: &[usize]) -> Self::Filtered;
}

/// Wraps a `MatrixProvider`, acting as if rows were removed from the inner matrix provider.
///
/// Used for deleting duplicate constraints after finding primal feasibility.
///
/// Only constraint rows should be deleted, variable bounds should remain intact.
/// TODO(DOCUMENTATION): Check the above property.
#[derive(PartialEq, Debug)]
pub struct RemoveRows<'a, MP> {
    provider: &'a MP,
    /// List of rows that this method removes.
    ///
    /// Sorted at all times.
    // TODO(OPTIMIZATION): Consider using a `HashMap`.
    pub rows_to_skip: Vec<usize>,
}

impl<'provider, MP> Filtered for RemoveRows<'provider, MP>
where
    MP: MatrixProvider<Column: IntoFilteredColumn>,
{
    fn filtered_rows(&self) -> &[usize] {
        &self.rows_to_skip
    }
}

impl<MP: 'static> ToFiltered for MP
where
    MP: MatrixProvider<Column: IntoFilteredColumn>,
    <<MP as MatrixProvider>::Column as IntoFilteredColumn>::Filtered: OrderedColumn,
{
    type Filtered<'provider> = RemoveRows<'provider, MP>;

    /// Create a new `RemoveRows` instance.
    ///
    /// # Arguments
    ///
    /// * `provider`: Reference to an instance implementing the `MatrixProvider` trait. Rows from
    /// this provider will be removed.
    /// * `rows_to_skip`: A **sorted** list of rows that are skipped.
    ///
    /// # Return value
    ///
    /// A new `RemoveRows` instance.
    default fn to_filtered(&self, rows_to_skip: Vec<usize>) -> Self::Filtered<'_> {
        debug_assert!(rows_to_skip.is_sorted());
        debug_assert_eq!(rows_to_skip.iter().collect::<HashSet<_>>().len(), rows_to_skip.len());

        RemoveRows {
            provider: self,
            rows_to_skip,
        }
    }
}

impl<'provider, MP> RemoveRows<'provider, MP>
where
    MP: MatrixProvider,
{
    /// Create a new `RemoveRows` instance.
    ///
    /// # Arguments
    ///
    /// * `provider`: Reference to an instance implementing the `MatrixProvider` trait. Rows from
    /// this provider will be removed.
    /// * `rows_to_skip`: A **sorted** list of rows that are skipped.
    ///
    /// # Return value
    ///
    /// A new `RemoveRows` instance.
    pub fn new(provider: &'provider MP, rows_to_skip: Vec<usize>) -> Self {
        debug_assert!(rows_to_skip.is_sorted());
        debug_assert_eq!(rows_to_skip.iter().collect::<HashSet<_>>().len(), rows_to_skip.len());

        RemoveRows {
            provider,
            rows_to_skip,
        }
    }

    /// Get the index of the same row in the original `MatrixProvider`.
    ///
    /// # Arguments
    ///
    /// * `i`: Index of row in the version of the problem from which rows were removed (this
    /// struct).
    ///
    /// # Return value
    ///
    /// Index of row in the original problem.
    #[must_use]
    pub fn get_underlying_row_index(&self, i: usize) -> usize {
        debug_assert!(i < self.provider.nr_rows() - self.rows_to_skip.len());

        Self::get_underlying_index(&self.rows_to_skip, i)
    }

    /// Delete a row
    ///
    /// # Arguments
    ///
    /// * `i`: Index of row in the version of the problem from which rows were removed (this
    /// struct).
    pub fn delete_row(&mut self, i: usize) {
        debug_assert!(i < self.provider.nr_rows() - self.rows_to_skip.len());

        let in_original_problem = self.get_underlying_row_index(i);
        debug_assert!(self.rows_to_skip.contains(&in_original_problem));
        let insertion_index = match self.rows_to_skip.binary_search(&in_original_problem) {
            Ok(_) => panic!("Deleting a row that already was deleted!"),
            Err(nr) => nr,
        };
        self.rows_to_skip.insert(insertion_index, in_original_problem);
    }

    /// How many constraints were removed.
    #[must_use]
    pub fn nr_constraints_deleted(&self) -> usize {
        self.rows_to_skip.len()
    }

    /// Method abstracting over the row and column getter methods.
    ///
    /// # Arguments
    ///
    /// * `i`: Index in the reduced version of the problem.
    ///
    /// # Return value
    ///
    /// Index in the original problem.
    fn get_underlying_index(skip_indices_array: &[usize], i: usize) -> usize {
        if skip_indices_array.is_empty() {
            // If no indices have been deleted, it's just the original value
            i
        } else if skip_indices_array.len() == 1 {
            // If one index has been deleted, see if that was before of after the value tested
            if i < skip_indices_array[0] {
                i
            } else {
                i + skip_indices_array.len()
            }
        } else {
            // skip_indices_array.len() >= 2
            if i < skip_indices_array[0] {
                i
            } else {
                // Binary search with invariants:
                //   1. skip_indices_array[lower_bound] - lower_bound <= i
                //   2. skip_indices_array[upper_bound] - upper_bound > i
                let (mut lower_bound, mut upper_bound) = (0, skip_indices_array.len());
                while upper_bound - lower_bound != 1 {
                    let middle = (lower_bound + upper_bound) / 2;
                    if skip_indices_array[middle] - middle <= i {
                        lower_bound = middle
                    } else {
                        upper_bound = middle
                    }
                }

                i + upper_bound
            }
        }
    }

    /// Method abstracting over the row and column deletion methods.
    ///
    /// # Arguments
    ///
    /// * `i`: Index in the reduced version of the problem to be deleted from the original problem.
    fn delete_index(skip_indices_array: &mut Vec<usize>, i: usize) {
        let in_original_problem = RemoveRows::<MP>::get_underlying_index(skip_indices_array, i);
        debug_assert!(skip_indices_array.contains(&in_original_problem));

        let insertion_index = match skip_indices_array.binary_search(&in_original_problem) {
            Ok(_) => panic!("Deleting an index that already was deleted."),
            Err(nr) => nr,
        };
        skip_indices_array.insert(insertion_index, in_original_problem);
    }
}

impl<'provider, MP> MatrixProvider for RemoveRows<'provider, MP>
where
    MP: MatrixProvider<Column: IntoFilteredColumn>,
{
    type Column = <MP::Column as IntoFilteredColumn>::Filtered;
    type Cost<'a> = MP::Cost<'a>;
    type Rhs = MP::Rhs;

    fn column(&self, j: usize) -> Self::Column {
        debug_assert!(j < self.nr_columns());

        self.provider.column(j).into_filtered(&self.rows_to_skip)
    }

    fn cost_value(&self, j: usize) -> Self::Cost<'_> {
        debug_assert!(j < self.nr_columns());

        self.provider.cost_value(j)
    }

    fn right_hand_side(&self) -> DenseVector<Self::Rhs> {
        let mut all = self.provider.right_hand_side();
        all.remove_indices(&self.rows_to_skip);
        all
    }

    fn bound_row_index(&self, j: usize, bound_type: BoundDirection) -> Option<usize> {
        debug_assert!(j < self.nr_columns());
        debug_assert!(self.rows_to_skip.iter().all(|&i| i < self.provider.nr_constraints()));

        // Only constraint rows are deleted,
        self.provider.bound_row_index(j, bound_type).map(|nr| nr - self.nr_constraints_deleted())
    }

    /// This implementation assumes that only constraint rows are removed from the `MatrixProvider`.
    fn nr_constraints(&self) -> usize {
        self.provider.nr_constraints() - self.nr_constraints_deleted()
    }

    /// This implementation assumes that only constraint rows are removed from the `MatrixProvider`.
    fn nr_variable_bounds(&self) -> usize {
        debug_assert!(self.rows_to_skip.iter().all(|&i| i < self.provider.nr_constraints()));

        self.provider.nr_variable_bounds()
    }

    fn nr_rows(&self) -> usize {
        self.provider.nr_rows() - self.nr_constraints_deleted()
    }

    fn nr_columns(&self) -> usize {
        self.provider.nr_columns()
    }

    fn reconstruct_solution<G>(&self, column_values: SparseVector<G, G>) -> SparseVector<G, G>
    where
        G: SparseElement<G> + SparseComparator,
    {
        self.provider.reconstruct_solution(column_values)
    }
}

impl<'provider, MP> PartialInitialBasis for RemoveRows<'provider, MP>
where
    MP: MatrixProvider<Column: IntoFilteredColumn> + PartialInitialBasis,
{
    fn pivot_element_indices(&self) -> Vec<(usize, usize)> {
        let mut from_parent = self.provider.pivot_element_indices();
        remove_sparse_indices(&mut from_parent, &self.rows_to_skip);
        from_parent
    }

    fn nr_initial_elements(&self) -> usize {
        // TODO(PERFORMANCE): Avoid the recomputation of the elements.
        self.pivot_element_indices().len()
    }
}

impl<'provider, MP> FeasibilityLogic for RemoveRows<'provider, MP>
where
    MP: MatrixProvider<Column: IntoFilteredColumn> + FeasibilityLogic,
{
    fn is_feasible(&self, j: usize, value: <MP::Column as Column>::F) -> bool {
        debug_assert!(j < self.nr_columns());

        self.provider.is_feasible(j, value)
    }

    fn closest_feasible(
        &self,
        j: usize,
        value: <MP::Column as Column>::F,
    ) -> (Option<<MP::Column as Column>::F>, Option<<MP::Column as Column>::F>) {
        debug_assert!(j < self.nr_columns());

        self.provider.closest_feasible(j, value)
    }
}

impl<'provider, MP> Display for RemoveRows<'provider, MP>
where
    MP: MatrixProvider<Column: IntoFilteredColumn>,
    <<MP as MatrixProvider>::Column as IntoFilteredColumn>::Filtered: OrderedColumn,
{
    fn fmt(&self, f: &mut Formatter<'_>) -> fmt::Result {
        let width = 8;

        f.write_str(&" ".repeat(width))?;
        for column in 0..self.nr_columns() {
            write!(f, "{:^width$}", column, width = width)?;
        }
        writeln!(f)?;
        f.write_str(&"-".repeat((1 + self.nr_columns()) * width))?;

        for row in 0..self.nr_rows() {
            write!(f, "{:>width$}", format!("{} |", row), width = width)?;
            for j in 0..self.nr_columns() {
                 write!(f, "{:^width$}", self.column(j).index_to_string(row), width = width)?;
            }
            writeln!(f)?;
        }
        writeln!(f)
    }
}

#[cfg(test)]
mod test {
    use relp_num::Rational32;

    use crate::algorithm::two_phase::matrix_provider::filter::generic_wrapper::RemoveRows;
    use crate::algorithm::two_phase::matrix_provider::matrix_data::MatrixData;

    #[test]
    fn get_underlying_index() {
        type T = Rational32;

        for (deleted, size) in vec![
            (vec![2, 5, 7, 9, 12, 15, 16, 19, 20, 21], 25),
            (vec![2], 5),
            (vec![2, 3], 6),
        ] {
            let left_from_original = (0..size).filter(|i| !deleted.contains(i)).collect::<Vec<_>>();
            for (in_reduced, in_original) in (0..left_from_original.len()).zip(left_from_original.into_iter()) {
                assert_eq!(RemoveRows::<MatrixData<T>>::get_underlying_index(&deleted, in_reduced), in_original);
            }
        }
    }
}
