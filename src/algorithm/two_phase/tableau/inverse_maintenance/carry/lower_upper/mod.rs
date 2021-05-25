//! # LU decomposition
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fmt;
use std::fmt::Display;
use std::iter;

use num_traits::Zero;

use crate::algorithm::two_phase::matrix_provider::column::{Column, OrderedColumn};
use crate::algorithm::two_phase::tableau::inverse_maintenance::{ColumnComputationInfo, ops};
use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::BasisInverse;
use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::lower_upper::eta_file::EtaFile;
use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::lower_upper::permutation::{FullPermutation, Permutation, RotateToBackPermutation};
use crate::data::linear_algebra::traits::SparseElement;
use crate::data::linear_algebra::vector::{SparseVector, Vector};

mod decomposition;
mod permutation;
mod eta_file;

/// Decompose a matrix `B` into `PBQ = LU` where
///
/// * `P` is a row permutation
/// * `Q` is a column permutation
/// * `L` is lower triangular with `1`'s on the diagonal
/// * `U` is upper triangular
///
/// Note that permutations `P` and `Q` have the transpose equal to their inverse as they are
/// orthogonal matrices.
///
/// `P` and `Q` are "full" permutations, not to be confused with the simpler "rotating" permutations
/// in the `updates` field.
#[derive(Eq, PartialEq, Clone, Debug)]
pub struct LUDecomposition<F> {
    /// Row permutation `P`.
    ///
    /// The `forward` application of the permutation to rows of `M` corresponds to `PM`.
    row_permutation: FullPermutation,
    /// Column permutation `Q`.
    ///
    /// The
    column_permutation: FullPermutation,
    /// Lower triangular matrix `L`.
    ///
    /// Column major, one's on diagonal implied. So the length of the vector is `m - 1`.
    lower_triangular: Vec<Vec<(usize, F)>>,
    /// Upper triangular matrix `U`.
    ///
    /// Column major.
    // TODO(PERFORMANCE): Consider storing the diagonal separately.
    upper_triangular: Vec<Vec<(usize, F)>>,

    // TODO(PERFORMANCE): Consider an Option<Permutation> instead to avoid doing identities
    // TODO(ARCHITECTURE): Consider grouping pivot index and eta file in a single update struct.
    updates: Vec<(EtaFile<F>, RotateToBackPermutation)>,
}

impl<F> BasisInverse for LUDecomposition<F>
where
    F: ops::Field + ops::FieldHR,
{
    type F = F;
    type ColumnComputationInfo = ColumnAndSpike<Self::F>;

    fn identity(m: usize) -> Self {
        Self {
            row_permutation: FullPermutation::identity(m),
            column_permutation: FullPermutation::identity(m),
            lower_triangular: vec![vec![]; m - 1],
            upper_triangular: (0..m).map(|i| vec![(i, F::one())]).collect(),
            updates: vec![],
        }
    }

    fn invert<C: Column + OrderedColumn>(columns: Vec<C>) -> Self
    where
        Self::F: ops::Column<C::F>,
    {
        let m = columns.len();
        let mut rows = vec![Vec::new(); m];
        for (j, column) in columns.into_iter().enumerate() {
            for (i, value) in column.iter() {
                rows[*i].push((j, value.into()));
            }
        }
        debug_assert!(rows.iter().all(|row| row.is_sorted_by_key(|&(j, _)| j)));

        Self::rows(rows)
    }

    fn change_basis(
        &mut self,
        pivot_row_index: usize,
        column: Self::ColumnComputationInfo,
    ) {
        let m = self.m();

        let pivot_column_index = {
            // Column with a pivot in `pivot_row_index` is leaving
            let mut pivot_column_index = pivot_row_index;
            self.column_permutation.forward(&mut pivot_column_index);
            for (_, q) in &self.updates {
                Permutation::forward(q, &mut pivot_column_index);
            }
            pivot_column_index
        };

        // Compute r
        let (u_bar, indices_to_zero): (Vec<_>, Vec<_>) = ((pivot_column_index + 1)..m)
            .filter_map(|j| {
                self.upper_triangular[j]
                    .binary_search_by_key(&pivot_column_index, |&(i, _)| i)
                    .ok()
                    .map(|data_index| (
                        (j, self.upper_triangular[j][data_index].1.clone()),
                        (j, data_index),
                    ))
            })
            .unzip();
        let u_bar = u_bar.into_iter().collect();
        let r = self.invert_upper_left(u_bar);
        let eta_file = EtaFile::new(r, pivot_column_index, self.m());

        // We now have all information needed to recover the triangle, start modifying it
        // Zero out part of a row that will be rotated to the bottom
        for (j, data_index) in indices_to_zero {
            self.upper_triangular[j].remove(data_index);
        }
        let Self::ColumnComputationInfo { column: _, mut spike } = column;
        debug_assert!(spike.iter().is_sorted_by_key(|&(i, _)| i));

        eta_file.update_spike_pivot_value(&mut spike);
        debug_assert!(
            spike.binary_search_by_key(&pivot_column_index, |&(i, _)| i).is_ok(),
            "This value should be present to avoid singularity because it will be the bottom corner value.",
        );
        // Insert the spike for upper
        self.upper_triangular[pivot_column_index] = spike;

        // Move the spike to the end
        self.upper_triangular[pivot_column_index..].rotate_left(1);

        // Rotate the "empty except for one value" pivot row to the bottom
        let q = RotateToBackPermutation::new(pivot_column_index, self.m());
        for j in pivot_column_index..m {
            q.forward_sorted(&mut self.upper_triangular[j])
        }

        // `upper_triangular` is diagonal again
        debug_assert!(self.upper_triangular.iter().enumerate().all(|(j, column)| {
            column.last().unwrap().0 == j && column.is_sorted_by_key(|&(i, _)| i)
        }));
        self.updates.push((eta_file, q));
    }

    fn generate_column<C: Column>(&self, original_column: C) -> Self::ColumnComputationInfo
    where
        Self::F: ops::Column<C::F>,
    {
        let rhs = original_column.iter()
            .map(|(mut i, v)| {
                self.row_permutation.forward(&mut i);
                (i, v.into())
            })
            // Also sorts after the row permutation
            .collect::<BTreeMap<_, _>>();
        let mut w = self.invert_lower_right(rhs);

        for (eta, q) in &self.updates {
            eta.apply_right(&mut w);
            // Q^-1 c = (c^T Q)^T because Q orthogonal matrix. Here we
            q.forward_sorted(&mut w);
        }

        let spike = w.clone();

        let mut column = self.invert_upper_right(w.into_iter().collect());

        for (_, q) in self.updates.iter().rev() {
            q.backward_unsorted(&mut column);
        }
        self.column_permutation.backward_unsorted(&mut column);
        column.sort_unstable_by_key(|&(i, _)| i);

        Self::ColumnComputationInfo {
            column: SparseVector::new(column, self.m()),
            spike,
        }
    }

    fn generate_element<C: Column + OrderedColumn>(&self, i: usize, original_column: C) -> Option<Self::F>
    where
        Self::F: ops::Column<C::F>,
    {
        self.generate_column(original_column).into_column().get(i).cloned()
    }

    fn should_refactor(&self) -> bool {
        // TODO(ENHANCEMENT): What would be a good decision rule?
        self.updates.len() > 10
    }

    fn basis_inverse_row(&self, mut row: usize) -> SparseVector<Self::F, Self::F> {
        self.column_permutation.forward(&mut row);
        for (_, q) in &self.updates {
            q.forward(&mut row);
        }

        let initial_rhs = iter::once((row, Self::F::one())).collect();
        let mut w = self.invert_upper_left(initial_rhs);

        for (eta, q) in self.updates.iter().rev() {
            q.backward_sorted(&mut w);
            eta.apply_left(&mut w);
        }

        let mut tuples = self.invert_lower_left(w.into_iter().collect());
        self.row_permutation.backward_sorted(&mut tuples);

        SparseVector::new(tuples, self.m())
    }

    fn m(&self) -> usize {
        self.row_permutation.len()
        // == self.column_permutation.len()
        // == self.upper_triangular.len()
        // == self.lower_triangular.len() + 1
    }
}

impl<F> LUDecomposition<F>
where
    F: ops::Field + ops::FieldHR,
{
    fn invert_lower_right(&self, mut rhs: BTreeMap<usize, F>) -> Vec<(usize, F)> {
        let mut result = Vec::new();

        while let Some((row, rhs_value)) = rhs.pop_first() {
            let column = row;

            let result_value = rhs_value; // Implicit 1 on the diagonal.

            if column != self.m() - 1 {
                // Introduce new nonzeros
                for &(row, ref value) in &self.lower_triangular[column] {
                    insert_or_shift_maybe_remove(row, &result_value * value, &mut rhs);
                }
            }

            result.push((row, result_value));
        }

        result
    }

    fn invert_upper_right(&self, mut rhs: BTreeMap<usize, F>) -> Vec<(usize, F)> {
        let mut result = Vec::new();

        while let Some((row, rhs_value)) = rhs.pop_last() {
            let column = row;
            let result_value = self.compute_result(column, rhs_value);
            self.update_rhs(column, &result_value, &mut rhs);
            result.push((row, result_value));
        }

        result.reverse();
        debug_assert!(result.is_sorted_by_key(|&(i, _)| i));

        result
    }

    fn invert_upper_right_row(&self, mut rhs: BTreeMap<usize, F>, target_row: usize) -> Option<F> {
        while let Some((row, rhs_value)) = rhs.pop_last() {
            let column = row;
            match row.cmp(&target_row) {
                Ordering::Less => return None,
                Ordering::Equal => return Some(self.compute_result(column, rhs_value)),
                Ordering::Greater => {
                    let result_value = self.compute_result(column, rhs_value);
                    self.update_rhs(column, &result_value, &mut rhs);
                },
            }
        }

        None
    }

    fn compute_result(&self, column: usize, rhs_value: F) -> F {
        debug_assert_eq!(
            self.upper_triangular[column].last().unwrap().0,
            column, "Needs to have a diagonal element",
        );

        let diagonal_item = &self.upper_triangular[column].last().unwrap().1;
        rhs_value / diagonal_item
    }

    fn update_rhs(&self, column: usize, result_value: &F, rhs: &mut BTreeMap<usize, F>) {
        let nr_column_items = self.upper_triangular[column].len();
        for &(row, ref value) in &self.upper_triangular[column][..(nr_column_items - 1)] {
            insert_or_shift_maybe_remove(row, result_value * value, rhs);
        }
    }

    fn invert_lower_left(&self, mut rhs: BTreeMap<usize, F>) -> Vec<(usize, F)> {
        let mut result = Vec::new();

        while let Some((column, rhs_value)) = rhs.pop_last() {
            let row = column;
            let result_value = rhs_value;

            // Introduce new nonzeros, by rows
            for j in 0..column {
                // TODO(PERFORMANCE): Avoid scanning all columns for row values
                let has_row = self.lower_triangular[j].binary_search_by_key(&row, |&(i, _)| i);
                if let Ok(data_index) = has_row {
                    let value = &self.lower_triangular[j][data_index].1;
                    insert_or_shift_maybe_remove(j, &result_value * value, &mut rhs);
                }
            }

            result.push((row, result_value));
        }

        result.reverse();
        debug_assert!(result.is_sorted_by_key(|&(i, _)| i));

        result
    }

    fn invert_upper_left(&self, mut rhs: BTreeMap<usize, F>) -> Vec<(usize, F)> {
        let mut result = Vec::new();

        while let Some((column, rhs_value)) = rhs.pop_first() {
            let row = column;
            let diagonal_item = &self.upper_triangular[column].last().unwrap().1;
            let result_value = rhs_value / diagonal_item;

            // Introduce new nonzeros, by rows
            for j in (column + 1)..self.m() {
                // TODO(PERFORMANCE): Avoid scanning all columns for row values
                let has_row = self.upper_triangular[j].binary_search_by_key(&row, |&(i, _)| i);
                if let Ok(data_index) = has_row {
                    let value = &self.upper_triangular[j][data_index].1;
                    insert_or_shift_maybe_remove(j, &result_value * value, &mut rhs);
                }
            }

            result.push((column, result_value));
        }

        debug_assert!(result.windows(2).all(|w| w[0].0 < w[1].0));

        result
    }
}

fn insert_or_shift_maybe_remove<F>(index: usize, change: F, rhs: &mut BTreeMap<usize, F>)
where
    F: ops::Field + ops::FieldHR + Zero,
{
    match rhs.get_mut(&index) {
        None => {
            rhs.insert(index, -change);
        }
        Some(existing) => {
            *existing -= change;
            if existing.is_zero() {
                rhs.remove(&index);
            }
        }
    }
}

/// The generated column and the spike that was generated in the process wrapped together.
#[derive(Debug)]
pub struct ColumnAndSpike<F> {
    column: SparseVector<F, F>,
    spike: Vec<(usize, F)>,
}

impl<F: SparseElement<F>> ColumnComputationInfo<F> for ColumnAndSpike<F> {
    fn column(&self) -> &SparseVector<F, F> {
        &self.column
    }

    fn into_column(self) -> SparseVector<F, F> {
        self.column
    }
}

impl<F> Display for LUDecomposition<F>
where
    F: ops::Field + ops::FieldHR,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let width = 10;
        let column_width = 3;

        writeln!(f, "Lower:")?;
        write!(f, "{:>width$} |", "", width = column_width)?;
        for j in 0..self.m() {
            write!(f, "{0:^width$}", j, width = width)?;
        }
        writeln!(f)?;
        let total_width = column_width + 1 + 1 + self.m() * width;
        writeln!(f, "{}", "-".repeat(total_width))?;

        for i in 0..self.m() {
            write!(f, "{0:>width$} |", i, width = column_width)?;
            for j in 0..self.m() {
                let value = match j.cmp(&i) {
                    Ordering::Greater => "".to_string(),
                    Ordering::Equal => "1".to_string(),
                    Ordering::Less => {
                        match self.lower_triangular[j].binary_search_by_key(&i, |&(i, _)| i) {
                            Ok(index) => self.lower_triangular[j][index].1.to_string(),
                            Err(_) => "0".to_string(),
                        }
                    }
                };
                write!(f, "{0:^width$}", value, width = width)?;
            }
            writeln!(f)?;
        }
        writeln!(f)?;

        writeln!(f, "Upper:")?;
        write!(f, "{:>width$} |", "", width = column_width)?;
        for j in 0..self.m() {
            write!(f, "{0:^width$}", j, width = width)?;
        }
        writeln!(f)?;
        writeln!(f, "{}", "-".repeat(total_width))?;

        for i in 0..self.m() {
            write!(f, "{0:>width$} |", i, width = column_width)?;
            for j in 0..self.m() {
                let value = match j.cmp(&i) {
                    Ordering::Equal | Ordering::Greater => {
                        match self.upper_triangular[j].binary_search_by_key(&i, |&(i, _)| i) {
                            Ok(index) => self.upper_triangular[j][index].1.to_string(),
                            Err(_) => "0".to_string(),
                        }
                    },
                    Ordering::Less => "".to_string(),
                };
                write!(f, "{0:^width$}", value, width = width)?;
            }
            writeln!(f)?;
        }
        writeln!(f)?;

        writeln!(f, "Row permutation:")?;
        writeln!(f, "{:?}", self.row_permutation)?;
        writeln!(f, "Column permutation:")?;
        writeln!(f, "{:?}", self.column_permutation)?;

        writeln!(f, "Updates:")?;
        for (i, (eta, t)) in self.updates.iter().enumerate() {
            writeln!(f, "Update {}: ", i)?;
            writeln!(f, "R: {:?}", eta)?;
            writeln!(f, "pivot index: {}", t)?;
        }
        writeln!(f)
    }
}

#[cfg(test)]
mod test {
    use std::collections::BTreeMap;

    use relp_num::{R64, RB};
    use relp_num::{Rational64, RationalBig};

    use crate::algorithm::two_phase::matrix_provider::column::identity::IdentityColumnStruct;
    use crate::algorithm::two_phase::matrix_provider::matrix_data;
    use crate::algorithm::two_phase::matrix_provider::matrix_data::Column;
    use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::BasisInverse;
    use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::lower_upper::ColumnAndSpike;
    use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::lower_upper::LUDecomposition;
    use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::lower_upper::permutation::FullPermutation;
    use crate::algorithm::two_phase::tableau::inverse_maintenance::ColumnComputationInfo;
    use crate::data::linear_algebra::vector::{SparseVector, Vector};

    mod matmul {
        use super::*;

        #[test]
        fn identity_empty() {
            let identity = LUDecomposition::<RationalBig>::identity(2);

            let column = BTreeMap::new();
            let result = identity.invert_upper_right(column);
            assert!(result.is_empty());
            let column = BTreeMap::new();
            let result = identity.invert_upper_left(column);
            assert!(result.is_empty());
            let column = BTreeMap::new();
            let result = identity.invert_lower_right(column);
            assert!(result.is_empty());
            let column = BTreeMap::new();
            let result = identity.invert_lower_left(column);
            assert!(result.is_empty());
        }

        #[test]
        fn identity_single() {
            let identity = LUDecomposition::<RationalBig>::identity(2);

            let column = vec![(0, RB!(1))];
            let result = identity.invert_upper_right(column.clone().into_iter().collect());
            assert_eq!(result, column);
            let column = vec![(0, RB!(1))];
            let result = identity.invert_upper_left(column.clone().into_iter().collect());
            assert_eq!(result, column);
            let column = vec![(0, RB!(1))];
            let result = identity.invert_lower_right(column.clone().into_iter().collect());
            assert_eq!(result, column);
            let column = vec![(0, RB!(1))];
            let result = identity.invert_lower_left(column.clone().into_iter().collect());
            assert_eq!(result, column);

            let column = vec![(1, RB!(1))];
            let result = identity.invert_upper_right(column.clone().into_iter().collect());
            assert_eq!(result, column);
            let column = vec![(1, RB!(1))];
            let result = identity.invert_upper_left(column.clone().into_iter().collect());
            assert_eq!(result, column);
            let column = vec![(1, RB!(1))];
            let result = identity.invert_lower_right(column.clone().into_iter().collect());
            assert_eq!(result, column);
            let column = vec![(1, RB!(1))];
            let result = identity.invert_lower_left(column.clone().into_iter().collect());
            assert_eq!(result, column);
        }

        #[test]
        fn identity_double() {
            let identity = LUDecomposition::<RationalBig>::identity(2);

            let column = vec![(0, RB!(1)), (1, RB!(1))];
            let result = identity.invert_upper_right(column.clone().into_iter().collect());
            assert_eq!(result, column);
            let column = vec![(0, RB!(1)), (1, RB!(1))];
            let result = identity.invert_upper_left(column.clone().into_iter().collect());
            assert_eq!(result, column);
            let column = vec![(0, RB!(1)), (1, RB!(1))];
            let result = identity.invert_lower_right(column.clone().into_iter().collect());
            assert_eq!(result, column);
            let column = vec![(0, RB!(1)), (1, RB!(1))];
            let result = identity.invert_lower_left(column.clone().into_iter().collect());
            assert_eq!(result, column);
        }

        #[test]
        fn offdiagonal_empty() {
            let offdiag = LUDecomposition {
                row_permutation: FullPermutation::identity(2),
                column_permutation: FullPermutation::identity(2),
                lower_triangular: vec![vec![(1, RB!(1))]],
                upper_triangular: vec![vec![(0, RB!(1))], vec![(1, RB!(1))]],
                updates: vec![],
            };

            let column: Column<Rational64> = Column::Sparse {
                constraint_values: vec![],
                slack: None,
                mock_array: [],
            };
            let result = offdiag.generate_column(column);
            assert_eq!(result.column, SparseVector::new(vec![], 2));
        }

        #[test]
        fn offdiagonal_single() {
            let offdiag = LUDecomposition {
                row_permutation: FullPermutation::identity(2),
                column_permutation: FullPermutation::identity(2),
                lower_triangular: vec![vec![(1, RB!(1))]],
                upper_triangular: vec![vec![(0, RB!(1))], vec![(1, RB!(1))]],
                updates: vec![],
            };

            let column = Column::Sparse {
                constraint_values: vec![(0, R64!(1))],
                slack: None,
                mock_array: [],
            };
            let result = offdiag.generate_column(column);
            assert_eq!(result.column, SparseVector::new(vec![(0, RB!(1)), (1, RB!(-1))], 2));

            let column = Column::Sparse {
                constraint_values: vec![(1, R64!(1))],
                slack: None,
                mock_array: [],
            };
            let result = offdiag.generate_column(column);
            assert_eq!(result.column, SparseVector::new(vec![(1, RB!(1))], 2));
        }
    }

    mod change_basis {
        use relp_num::One;

        use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::lower_upper::eta_file::EtaFile;
        use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::lower_upper::permutation::RotateToBackPermutation;

        use super::*;

        /// Spike is the column which is already there.
        #[test]
        fn no_change() {
            let mut initial = LUDecomposition::<RationalBig>::identity(3);

            let spike = vec![(1, RB!(1))];
            let column_computation_info = ColumnAndSpike {
                column: SparseVector::new(spike.clone(), 3),
                spike,
            };
            initial.change_basis(1, column_computation_info);
            let modified = initial;

            let mut expected = LUDecomposition::<RationalBig>::identity(3);
            expected.updates.push((
                EtaFile::new(vec![], 1, 3),
                RotateToBackPermutation::new(1, 3),
            ));
            assert_eq!(modified, expected);
        }

        #[test]
        fn from_identity_2() {
            let mut identity = LUDecomposition::<RationalBig>::identity(2);

            let spike = vec![(0, RB!(1)), (1, RB!(1))];
            let column_computation_info = ColumnAndSpike {
                column: SparseVector::new(spike.clone(), 2),
                spike,
            };
            identity.change_basis(0, column_computation_info);
            let modified = identity;

            let expected = LUDecomposition {
                row_permutation: FullPermutation::identity(2),
                column_permutation: FullPermutation::identity(2),
                lower_triangular: vec![vec![]],
                upper_triangular: vec![vec![(0, RB!(1))], vec![(0, RB!(1)), (1, RB!(1))]],
                updates: vec![(
                    EtaFile::new(vec![], 0, 2),
                    RotateToBackPermutation::new(0, 2),
                )],
            };
            assert_eq!(modified, expected);
        }

        /// Doesn't require any `r`, permutations only are sufficient
        #[test]
        fn from_5x5_identity_no_r() {
            let m = 5;
            let mut initial = LUDecomposition::<RationalBig>::identity(5);

            let spike = vec![(0, RB!(2)), (1, RB!(3)), (2, RB!(5)), (3, RB!(7))];
            let column_computation_info = ColumnAndSpike {
                column: SparseVector::new(spike.clone(), m),
                spike,
            };
            initial.change_basis(1, column_computation_info);
            let modified = initial;

            let expected = LUDecomposition {
                row_permutation: FullPermutation::identity(m),
                column_permutation: FullPermutation::identity(m),
                lower_triangular: vec![vec![]; m - 1],
                upper_triangular: vec![
                    vec![(0, RB!(1))],
                    vec![(1, RB!(1))],
                    vec![(2, RB!(1))],
                    vec![(3, RB!(1))],
                    vec![(0, RB!(2)), (1, RB!(5)), (2, RB!(7)), (4, RB!(3))],
                ],
                updates: vec![(
                    EtaFile::new(vec![], 1, m),
                    RotateToBackPermutation::new(1, m),
                )],
            };
            assert_eq!(modified, expected);
        }

        /// Does require an `r`, permutations are not sufficient.
        #[test]
        fn from_4x4_identity() {
            let m = 4;
            let mut initial = LUDecomposition {
                row_permutation: FullPermutation::identity(m),
                column_permutation: FullPermutation::identity(m),
                lower_triangular: vec![vec![]; m - 1],
                upper_triangular: vec![
                    vec![(0, RB!(1))],
                    vec![(1, RB!(1))],
                    vec![(2, RB!(4))],
                    vec![(1, RB!(5)), (3, RB!(6))],
                ],
                updates: vec![],
            };

            let spike = vec![(1, RB!(2)), (2, RB!(3)), (3, RB!(4))];
            let column_computation_info = ColumnAndSpike {
                // They are the same in this case, row permutation and lower are identity
                column: SparseVector::new(spike.clone(), m),
                spike,
            };
            initial.change_basis(1, column_computation_info);
            let modified = initial;

            let expected = LUDecomposition {
                row_permutation: FullPermutation::identity(m),
                column_permutation: FullPermutation::identity(m),
                lower_triangular: vec![vec![]; m - 1],
                upper_triangular: vec![
                    vec![(0, RB!(1))],
                    vec![(1, RB!(4))],
                    vec![(2, RB!(6))],
                    vec![(1, RB!(3)), (2, RB!(4)), (3, -RB!(8, 6))],
                ],
                updates: vec![(
                    EtaFile::new(vec![(3, RB!(5, 6))], 1, m),
                    RotateToBackPermutation::new(1, m),
                )],
            };
            assert_eq!(modified, expected);

            // Columns
            assert_eq!(
                modified.generate_column(IdentityColumnStruct((0, One))).into_column(),
                SparseVector::standard_basis_vector(0, m),
            );
            assert_eq!(
                modified.generate_column(IdentityColumnStruct((1, One))).into_column(),
                SparseVector::new(vec![(1, RB!(-3, 4)), (2, RB!(9, 16)), (3, RB!(1, 2))], m),
            );
            assert_eq!(
                modified.generate_column(IdentityColumnStruct((2, One))).into_column(),
                SparseVector::new(vec![(2, RB!(1, 4))], m),
            );
            assert_eq!(
                modified.generate_column(IdentityColumnStruct((3, One))).into_column(),
                SparseVector::new(vec![(1, RB!(5, 8)), (2, RB!(-15, 32)), (3, RB!(-1, 4))], m),
            );

            // Rows
            assert_eq!(
                modified.basis_inverse_row(0),
                SparseVector::standard_basis_vector(0, m),
            );
            assert_eq!(
                modified.basis_inverse_row(1),
                SparseVector::new(vec![(1, RB!(-3, 4)), (3, RB!(5, 8))], m),
            );
            assert_eq!(
                modified.basis_inverse_row(2),
                SparseVector::new(vec![(1, RB!(9, 16)), (2, RB!(1, 4)), (3, RB!(-15, 32))], m),
            );
            assert_eq!(
                modified.basis_inverse_row(3),
                SparseVector::new(vec![(1, RB!(1, 2)), (3, RB!(-1, 4))], m),
            );
        }

        /// From "A review of the LU update in the simplex algorithm" by Joseph M. Elble and
        /// Nikolaes V. Sahinidis, Int. J. Mathematics in Operational Research, Vol. 4, No. 4, 2012.
        #[test]
        fn from_5x5_identity() {
            let m = 5;
            let mut initial = LUDecomposition {
                row_permutation: FullPermutation::identity(m),
                column_permutation: FullPermutation::identity(m),
                lower_triangular: vec![vec![]; m - 1],
                upper_triangular: vec![
                    vec![(0, RB!(11))],
                    vec![(0, RB!(12)), (1, RB!(22))],
                    vec![(0, RB!(13)), (1, RB!(23)), (2, RB!(33))],
                    vec![(0, RB!(14)), (1, RB!(24)), (2, RB!(34)), (3, RB!(44))],
                    vec![(0, RB!(15)), (1, RB!(25)), (2, RB!(35)), (3, RB!(45)), (4, RB!(55))],
                ],
                updates: vec![],
            };

            let spike = vec![(0, RB!(12)), (1, RB!(22)), (2, RB!(32)), (3, RB!(42))];
            let column_computation_info = ColumnAndSpike {
                column: SparseVector::new(spike.clone(), m),
                spike,
            };
            initial.change_basis(1, column_computation_info);
            let modified = initial;

            let expected = LUDecomposition {
                row_permutation: FullPermutation::identity(m),
                column_permutation: FullPermutation::identity(m),
                lower_triangular: vec![vec![]; m - 1],
                upper_triangular: vec![
                    vec![(0, RB!(11))],
                    vec![(0, RB!(13)), (1, RB!(33))],
                    vec![(0, RB!(14)), (1, RB!(34)), (2, RB!(44))],
                    vec![(0, RB!(15)), (1, RB!(35)), (2, RB!(45)), (3, RB!(55))],
                    vec![(0, RB!(12)), (1, RB!(32)), (2, RB!(42)), (4, RB!(-215, 363))],
                ],
                updates: vec![(
                    EtaFile::new(vec![
                            (2, RB!(23, 33)),
                            (3, RB!(24 * 33 - 34 * 23, 33 * 44)),
                            (4, RB!(43, 7986)),
                        ], 1, m),
                    RotateToBackPermutation::new(1, m),
                )],
            };
            assert_eq!(modified, expected);

            // Columns
            assert_eq!(
                modified.generate_column(IdentityColumnStruct((0, One))).into_column(),
                SparseVector::new(vec![(0, RB!(1, 11))], m),
            );
            assert_eq!(
                modified.generate_column(IdentityColumnStruct((1, One))).into_column(),
                SparseVector::new(vec![(0, RB!(-2, 11)), (1, RB!(-363, 215)), (2, RB!(-1, 43)), (3, RB!(693, 430))], m),
            );
            assert_eq!(
                modified.generate_column(IdentityColumnStruct((2, One))).into_column(),
                SparseVector::new(vec![(0, RB!(1, 11)), (1, RB!(253, 215)), (2, RB!(2, 43)), (3, RB!(-483, 430))], m),
            );
            assert_eq!(
                modified.generate_column(IdentityColumnStruct((3, One))).into_column(),
                SparseVector::new(vec![(1, RB!(1, 86)), (2, RB!(-1, 43)), (3, RB!(1, 86))], m),
            );
            assert_eq!(
                modified.generate_column(IdentityColumnStruct((4, One))).into_column(),
                SparseVector::new(vec![(1, RB!(1, 110)), (3, RB!(-3, 110)), (4, RB!(1, 55))], m),
            );
            // Sum of two
            assert_eq!(
                modified.generate_column(matrix_data::Column::TwoSlack([(0, RB!(1)), (1, RB!(1))], [])).into_column(),
                SparseVector::new(vec![(0, RB!(-1, 11)), (1, RB!(-363, 215)), (2, RB!(-1, 43)), (3, RB!(693, 430))], m),
            );

            // Rows
            assert_eq!(
                modified.basis_inverse_row(0),
                SparseVector::new(vec![(0, RB!(1, 11)), (1, RB!(-2, 11)), (2, RB!(1, 11))], m),
            );
            assert_eq!(
                modified.basis_inverse_row(1),
                SparseVector::new(vec![(1, RB!(-363, 215)), (2, RB!(253, 215)), (3, RB!(1, 86)), (4, RB!(1, 110))], m),
            );
            assert_eq!(
                modified.basis_inverse_row(2),
                SparseVector::new(vec![(1, RB!(-1, 43)), (2, RB!(2, 43)), (3, RB!(-1, 43))], m),
            );
            assert_eq!(
                modified.basis_inverse_row(3),
                SparseVector::new(vec![(1, RB!(693, 430)), (2, RB!(-483, 430)), (3, RB!(1, 86)), (4, RB!(-3, 110))], m),
            );
            assert_eq!(
                modified.basis_inverse_row(4),
                SparseVector::new(vec![(4, RB!(1, 55))], m),
            );
        }
    }
}
