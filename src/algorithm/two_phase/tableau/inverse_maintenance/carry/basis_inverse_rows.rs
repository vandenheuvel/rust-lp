//! # Basis inverse rows
//!
//! Explicit row-major representation of the basis inverse B^-1. The inverse of a sparse matrix is
//! not generally sparse, so this is not a scalable algorithm. It is however useful for debugging
//! purposes to have an explicit representation of the basis inverse at hand.
use std::cmp::{Ordering, max};
use std::fmt;

use relp_num::One;

use crate::algorithm::two_phase::matrix_provider::column::{Column, OrderedColumn};
use crate::algorithm::two_phase::matrix_provider::column::identity::IdentityColumnStruct;
use crate::algorithm::two_phase::tableau::inverse_maintenance::{ColumnComputationInfo, ops};
use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::{BasisInverse, RemoveBasisPart};
use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::lower_upper::LUDecomposition;
use crate::algorithm::utilities::remove_indices;
use crate::data::linear_algebra::traits::SparseElement;
use crate::data::linear_algebra::vector::{SparseVector, Vector};

/// Explicit row-major sparse representation of the basis inverse.
#[derive(Eq, PartialEq, Debug)]
pub struct BasisInverseRows<F> {
    rows: Vec<SparseVector<F, F>>,
}

impl<F> BasisInverseRows<F>
where
    F: ops::Field + ops::FieldHR,
{
    /// Create a new instance by wrapping sparse vectors.
    #[must_use]
    pub fn new(rows: Vec<SparseVector<F, F>>) -> Self {
        Self { rows }
    }

    /// Normalize the pivot row.
    ///
    /// That is, the pivot value will be set to `1`.
    ///
    /// # Arguments
    ///
    /// * `pivot_row_index`: Index of the pivot row.
    /// * `column`: Column relative to the current basis to be entered into that basis.
    fn normalize_pivot_row(
        &mut self,
        pivot_row_index: usize,
        column: &SparseVector<F, F>,
    ) {
        let pivot_value = column.get(pivot_row_index)
            .expect("Pivot value can't be zero.");

        self.rows[pivot_row_index].element_wise_divide(pivot_value);
    }

    /// Normalize the pivot row and row reduce the other basis inverse rows.
    ///
    /// # Arguments
    ///
    /// * `pivot_row_index`: Index of the pivot row.
    /// * `column`: Column relative to the current basis to be entered into that basis.
    ///
    /// # Note
    ///
    /// This method requires a normalized pivot element.
    fn row_reduce(
        &mut self,
        pivot_row_index: usize,
        column: &SparseVector<F, F>,
    ) {
        debug_assert!(pivot_row_index < self.m());

        // TODO(OPTIMIZATION): Improve the below algorithm; when does SIMD kick in?
        let (rows_left, rows_right) = self.rows.split_at_mut(pivot_row_index);
        let (rows_middle, rows_right) = rows_right.split_first_mut().unwrap();

        for (edit_row_index, column_value) in column.iter_values() {
             match edit_row_index.cmp(&pivot_row_index) {
                Ordering::Less => rows_left[*edit_row_index]
                    .add_multiple_of_row(&-column_value, &rows_middle),
                Ordering::Equal => {},
                Ordering::Greater => rows_right[*edit_row_index - (pivot_row_index + 1)]
                    .add_multiple_of_row(&-column_value, &rows_middle),
            }
        }
    }

    fn m(&self) -> usize {
        self.rows.len()
    }
}

impl<F> BasisInverse for BasisInverseRows<F>
where
    F: ops::Field + ops::FieldHR,
{
    type F = F;
    type ColumnComputationInfo = SparseVector<Self::F, Self::F>;

    fn identity(m: usize) -> Self {
        Self {
            rows: (0..m).map(|i| SparseVector::standard_basis_vector(i, m)).collect(),
        }
    }

    fn invert<C: Column + OrderedColumn>(columns: Vec<C>) -> Self
    where
        Self::F: ops::Column<C::F>,
    {
        let m = columns.len();
        let lower_upper = LUDecomposition::invert(columns);

        let inverted_columns = (0..m)
            .map(|i| IdentityColumnStruct((i, One)))
            .map(|column| lower_upper.generate_column(column))
            .collect::<Vec<_>>();

        let mut row_major = vec![Vec::new(); m];
        for (j, column) in inverted_columns.into_iter().enumerate() {
            for (i, value) in column.into_column().into_iter() {
                row_major[i].push((j, value));
            }
        }

        let rows = row_major.into_iter()
            .map(|tuples| SparseVector::new(tuples, m))
            .collect();

        Self {
            rows,
        }
    }

    fn change_basis(
        &mut self,
        pivot_row_index: usize,
        column: Self::ColumnComputationInfo,
    ) {
        debug_assert!(pivot_row_index < self.m());
        debug_assert_eq!(column.column().len(), self.m());

        // The order of these calls matters: the first of the two normalizes the pivot row
        self.normalize_pivot_row(pivot_row_index, column.column());
        self.row_reduce(pivot_row_index, column.column());
    }

    fn generate_column<C: Column + OrderedColumn>(&self, original_column: C) -> Self::ColumnComputationInfo
    where
        Self::F: ops::Column<C::F>,
    {
        let tuples = (0..self.m())
            .map(|i| self.generate_element(i, original_column.clone()))
            .enumerate()
            .filter_map(|(i, v)| v.map(|inner| (i, inner)))
            .collect();

        SparseVector::new(tuples, self.m())
    }

    fn generate_element<C: Column + OrderedColumn>(
        &self,
        i: usize,
        original_column: C,
    ) -> Option<Self::F>
    where
        Self::F: ops::Column<C::F>,
    {
        debug_assert!(i < self.m());

        let element = self.rows[i].sparse_inner_product::<Self::F, _, _>(original_column.iter());
        if element.is_not_zero() {
            Some(element)
        } else {
            None
        }
    }

    fn should_refactor(&self) -> bool {
        // Refactoring is not a concept in this implementation, it doesn't degenerate through
        // updating, as there is no updating.
        false
    }

    fn basis_inverse_row(&self, row: usize) -> SparseVector<Self::F, Self::F> {
        self.rows[row].clone()
    }

    fn m(&self) -> usize {
        self.rows.len()
    }
}

impl<F> RemoveBasisPart for BasisInverseRows<F>
where
    F: SparseElement<F> + ops::Field + ops::FieldHR,
{
    fn remove_basis_part(&mut self, indices: &[usize]) {
        let old_m = self.m();
        debug_assert!(indices.len() < old_m);

        remove_indices(&mut self.rows, indices);
        // Remove the columns
        for element in &mut self.rows {
            element.remove_indices(indices);
        }

        debug_assert_eq!(self.m(), old_m - indices.len());
    }
}

impl<F: SparseElement<F>> ColumnComputationInfo<F> for SparseVector<F, F> {
    fn column(&self) -> &SparseVector<F, F> {
        &self
    }

    fn into_column(self) -> SparseVector<F, F> {
        self
    }
}

impl<F> fmt::Display for BasisInverseRows<F>
where
    F: ops::Field + ops::FieldHR,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let rows = (0..self.m()).map(|i| {
            (0..self.m()).map(|j| match self.rows[i].get(j) {
                None => "0".to_string(),
                Some(value) => value.to_string(),
            }).collect::<Vec<_>>()
        }).collect::<Vec<_>>();

        let row_counter_width = (self.m() - 1).to_string().len();
        let column_width = (0..self.m()).map(|j| {
            max(j.to_string().len(), (0..self.m()).map(|i| rows[i][j].len()).max().unwrap())
        }).collect::<Vec<_>>();

        // Column counters
        write!(f, "{0:>width$} |", "", width = row_counter_width)?;
        for (j, width) in column_width.iter().enumerate() {
            write!(f, " {0:^width$}", j, width = width)?;
        }
        writeln!(f)?;

        // Separator
        let total_width = (row_counter_width + 1) + 1 +
            column_width.iter().map(|l| 1 + l).sum::<usize>();
        writeln!(f, "{}", "-".repeat(total_width))?;

        // Row counter and row data
        for (i, row) in rows.into_iter().enumerate() {
            write!(f, "{0:>width$} |", i, width = row_counter_width)?;
            for (width, value) in column_width.iter().zip(row) {
                write!(f, " {0:^width$}", value, width = width)?;
            }
            writeln!(f)?;
        }
        writeln!(f)
    }
}

#[cfg(test)]
mod test {
    use relp_num::RationalBig;
    use relp_num::RB;

    use crate::algorithm::two_phase::matrix_provider::matrix_data;
    use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::basis_inverse_rows::BasisInverseRows;
    use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::BasisInverse;
    use crate::data::linear_algebra::vector::{SparseVector, Vector};

    #[test]
    fn invert_identity() {
        let columns = vec![
            matrix_data::Column::Slack([(0, RB!(1))], []),
            matrix_data::Column::Slack([(1, RB!(1))], []),
        ];

        let result = BasisInverseRows::<RationalBig>::invert(columns);
        let expected = BasisInverseRows {
            rows: vec![
                SparseVector::standard_basis_vector(0, 2),
                SparseVector::standard_basis_vector(1, 2),
            ],
        };
        assert_eq!(result, expected);
    }

    #[test]
    fn invert_non_identity() {
        let columns = vec![
            matrix_data::Column::TwoSlack([(0, RB!(1)), (1, RB!(1))], []),
            matrix_data::Column::Slack([(1, RB!(1))], []),
        ];

        let result = BasisInverseRows::<RationalBig>::invert(columns);
        let expected = BasisInverseRows {
            rows: vec![
                SparseVector::standard_basis_vector(0, 2),
                SparseVector::new(vec![(0, RB!(-1)), (1, RB!(1))], 2),
            ],
        };
        assert_eq!(result, expected);
    }
}
