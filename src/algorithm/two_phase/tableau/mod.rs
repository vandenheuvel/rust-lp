//! # Data structures for Simplex
//!
//! Contains the simplex tableau and logic for elementary operations which can be performed upon it.
//! The tableau is extended with supplementary data structures for efficiency.
use std::collections::HashSet;
use std::fmt::{Debug, Display, Formatter, Result as FormatResult};

use num_traits::One;
use num_traits::Zero;

use crate::algorithm::two_phase::matrix_provider::column::Column;
use crate::algorithm::two_phase::tableau::inverse_maintenance::{ColumnComputationInfo, InverseMaintener, ops as im_ops};
use crate::algorithm::two_phase::tableau::kind::Kind;
use crate::data::linear_algebra::vector::{SparseVector, Vector};
use std::cmp::max;

pub mod inverse_maintenance;
pub mod kind;

/// The most high-level data structure that is used by the Simplex algorithm: the Simplex tableau.
///
/// It holds only a reference to the (immutable) problem it solves, but owns the data structures
/// that describe the current solution basis.
#[derive(Eq, PartialEq, Debug)]
pub struct Tableau<IM, K> {
    /// Represents a matrix of size (m + 1) x (m + 1) (includes -pi, objective value, constraints).
    ///
    /// This attribute changes with a basis change.
    inverse_maintainer: IM,

    /// All columns currently in the basis.
    ///
    /// Could also be derived from `basis_indices`, but is here for faster reading and writing.
    basis_columns: HashSet<usize>,

    /// Whether this tableau has artificial variables (and is in the first phase of the two-phase
    /// algorithm) or not. See the `Kind` trait for more information.
    kind: K,
}

impl<IM, K> Tableau<IM, K>
where
    IM: InverseMaintener<F: im_ops::Column<<K::Column as Column>::F> + im_ops::Cost<K::Cost>>,
    K: Kind,
{
    /// Brings a column into the basis by updating the `self.carry` matrix and updating the
    /// data structures holding the collection of basis columns.
    pub fn bring_into_basis(
        &mut self,
        pivot_column_index: usize,
        pivot_row_index: usize,
        column: IM::ColumnComputationInfo,
        cost: IM::F
    ) {
        debug_assert!(pivot_column_index < self.kind.nr_columns());
        debug_assert!(pivot_row_index < self.nr_rows());

        let leaving_column = self.inverse_maintainer.change_basis(pivot_row_index, pivot_column_index, column, cost);
        self.update_basis_indices(pivot_column_index, leaving_column);
        self.inverse_maintainer.after_basis_change(&self.kind);
    }

    /// Update the basis index.
    ///
    /// Removes the index of the variable leaving the basis from the `basis_column_map` attribute,
    /// while inserting the entering variable index.
    ///
    /// # Arguments
    ///
    /// * `pivot_row`: Row index of the pivot, in range `0` until `self.nr_rows()`.
    /// * `pivot_column`: Column index of the pivot, in range `0` until `self.nr_columns()`. Is not
    /// yet in the basis.
    fn update_basis_indices(
        &mut self,
        pivot_column: usize,
        leaving_column: usize,
    ) {
        debug_assert!(pivot_column < self.nr_columns());
        debug_assert!(leaving_column < self.nr_columns());

        let was_there = self.basis_columns.remove(&leaving_column);
        debug_assert!(was_there);
        let was_not_there = self.basis_columns.insert(pivot_column);
        debug_assert!(was_not_there);
    }

    /// Calculates the relative cost of a column.
    ///
    /// # Arguments
    ///
    /// * `j`: Index of column to calculate the relative cost for, in range `0` through
    /// `self.nr_variables()`.
    ///
    /// # Return value
    ///
    /// The relative cost.
    ///
    /// # Note
    ///
    /// That column will typically not be a basis column. Although the method could be valid for
    /// those inputs as well, this should never be calculated, as the relative cost always equals
    /// zero in that situation.
    pub fn relative_cost(&self, j: usize) -> IM::F {
        debug_assert!(j < self.nr_columns());

        let initial = self.kind.initial_cost_value(j);
        let difference = self.inverse_maintainer.cost_difference(&self.kind.original_column(j));
        difference + initial
    }

    /// Column of original problem with respect to the current basis.
    ///
    /// Generate a column of the tableau as it would look like with the current basis by matrix
    /// multiplying the original column and the carry matrix.
    ///
    /// # Arguments
    ///
    /// * `j`: Column index of the variable, in range `0` until `self.nr_columns()`.
    ///
    /// # Return value
    ///
    /// `SparseVector<T>` of size `m`.
    pub fn generate_column(&self, j: usize) -> IM::ColumnComputationInfo {
        debug_assert!(j < self.nr_columns());

        self.inverse_maintainer.generate_column(self.kind.original_column(j))
    }

    /// Single element with respect to the current basis.
    pub fn generate_element(&self, i: usize, j: usize) -> Option<IM::F> {
        debug_assert!(i < self.nr_rows());
        debug_assert!(j < self.nr_columns());

        self.inverse_maintainer.generate_element(i, self.kind.original_column(j))
    }

    /// Whether a column is in the basis.
    ///
    /// # Return value
    ///
    /// `bool` with value true if the column is in the basis.
    ///
    /// # Note
    ///
    /// This method may not be accurate when there are artificial variables.
    pub fn is_in_basis(&self, column: usize) -> bool {
        debug_assert!(column < self.nr_columns());

        self.basis_columns.contains(&column)
    }

    /// Get the current basic feasible solution.
    ///
    /// # Return value
    ///
    /// `SparseVector<T>` of the solution.
    pub fn current_bfs(&self) -> SparseVector<IM::F, IM::F> {
        let tuples = self.inverse_maintainer.current_bfs();
        SparseVector::new(tuples, self.nr_columns())
    }

    /// Get the cost of the current solution.
    ///
    /// # Return value
    ///
    /// The current value of the objective function.
    ///
    /// # Note
    ///
    /// This function works for both artificial and non-artificial tableaus.
    pub fn objective_function_value(&self) -> IM::F {
        self.inverse_maintainer.get_objective_function_value()
    }

    /// Number of rows in the tableau.
    ///
    /// # Return value
    ///
    /// The number of rows.
    pub fn nr_rows(&self) -> usize {
        self.kind.nr_rows()
    }

    /// Number of variables in the problem.
    ///
    /// # Return value
    ///
    /// The number of variables.
    ///
    /// # Note
    ///
    /// This number might be extremely large, depending on the `AdaptedTableauProvider`.
    pub fn nr_columns(&self) -> usize {
        self.kind.nr_columns()
    }
}

impl<IM, K> Tableau<IM, K>
where
    IM: InverseMaintener<F: im_ops::FieldHR + im_ops::Column<<K::Column as Column>::F> + im_ops::Cost<K::Cost>>,
    K: Kind,
{
    /// Determine the row to pivot on.
    ///
    /// Determine the row to pivot on, given the column. This is the row with the positive but
    /// minimal ratio between the current constraint vector and the column.
    ///
    /// When there are multiple choices for the pivot row, Bland's anti cycling algorithm
    /// is used to avoid cycles.
    ///
    /// TODO(ARCHITECTURE): Reconsider the below; should this be a strategy?
    /// Because this method allows for less strategy and heuristics, it is not included in the
    /// `PivotRule` trait.
    ///
    /// # Arguments
    ///
    /// * `column`: Problem column with respect to the current basis with length `m`.
    ///
    /// # Return value
    ///
    /// Index of the row to pivot on. If not found, the problem is optimal.
    pub fn select_primal_pivot_row(&self, column: &SparseVector<IM::F, IM::F>) -> Option<usize> {
        debug_assert_eq!(Vector::len(column), self.nr_rows());

        // (chosen index, minimum ratio, corresponding leaving_column (for Bland's algorithm))
        let mut min_values: Option<(usize, IM::F, usize)> = None;
        for (row, xij) in column.iter_values() {
            if xij > &IM::F::zero() {
                let ratio = self.inverse_maintainer.get_constraint_value(*row) / xij;
                // Bland's anti cycling algorithm
                let leaving_column = self.inverse_maintainer.basis_column_index_for_row(*row);
                if let Some((min_index, min_ratio, min_leaving_column)) = &mut min_values {
                    if &ratio == min_ratio && leaving_column < *min_leaving_column {
                        *min_index = *row;
                        *min_leaving_column = leaving_column;
                    } else if &ratio < min_ratio {
                        *min_index = *row;
                        *min_ratio = ratio;
                        *min_leaving_column = leaving_column;
                    }
                } else {
                    min_values = Some((*row, ratio, leaving_column))
                }
            }
        }

        min_values.map(|(min_index, _, _)| min_index)
    }
}

/// Check whether the tableau currently has a valid basic feasible solution.
///
/// Only used for debug purposes.
pub fn debug_assert_in_basic_feasible_solution_state<IM, K>(tableau: &Tableau<IM, K>)
where
    IM: InverseMaintener<F: im_ops::Column<<K::Column as Column>::F> + im_ops::Cost<K::Cost>>,
    K: Kind,
{
    // Checking basis_columns
    // Correct number of basis columns (uniqueness is implied because it's a set)
    debug_assert_eq!(tableau.basis_columns.len(), tableau.nr_rows());

    // Checking carry matrix
    // `basis_inverse_rows` are a proper inverse by regenerating basis columns
    (0..tableau.nr_rows())
        .map(|i| (i, tableau.inverse_maintainer.basis_column_index_for_row(i)))
        .for_each(|(i, j)| {
            let e_i = SparseVector::new(vec![(i, IM::F::one())], tableau.nr_rows());
            debug_assert_eq!(
                tableau.generate_column(j).into_column(), e_i,
                "Column {} is not equal to e_{}", j, i,
            );
        });

    // `minus_pi` get to relative zero cost for basis columns
    (0..tableau.nr_rows())
        .map(|i| tableau.inverse_maintainer.basis_column_index_for_row(i))
        .for_each(|j| debug_assert_eq!(
            tableau.relative_cost(j), IM::F::zero(),
            "Relative cost of column {} is not zero", j,
        ));

    // `b` >= 0
    (0..tableau.nr_rows())
        .for_each(|i| {
            let value = &tableau.inverse_maintainer.b()[i];
            debug_assert!(
                value >= &IM::F::zero(),
                "rhs (b) is not always nonnegative: at index {} we have {} < 0", i, value,
            );
        });
}

impl<IM, K> Display for Tableau<IM, K>
where
    IM: InverseMaintener<F: im_ops::Column<<K::Column as Column>::F> + im_ops::Cost<K::Cost>>,
    K: Kind,
{
    fn fmt(&self, f: &mut Formatter) -> FormatResult {
        assert!(self.nr_rows().to_string().len() <= "cost".len());

        writeln!(f, "=== Tableau ===")?;
        let objective = (-self.inverse_maintainer.get_objective_function_value()).to_string();
        let cost = (0..self.nr_columns()).map(|j| {
            self.relative_cost(j).to_string()
        }).collect::<Vec<_>>();
        let b = self.inverse_maintainer.b()
            .data.iter().map(|v| v.to_string())
            .collect::<Vec<_>>();
        let columns = (0..self.nr_columns()).map(|j| {
            let column = self.generate_column(j);
            (0..self.nr_rows())
                .map(|i| match column.column().get(i) {
                    None => "".to_string(),
                    Some(value) => value.to_string(),
                })
                .collect::<Vec<_>>()
        }).collect::<Vec<_>>();

        let row_counter_width = "cost".len();
        let column_width = columns.iter().enumerate().map(|(j, column)| {
            [
                column.iter().map(String::len).max().unwrap(),
                j.to_string().len(),
                cost[j].len(),
            ].iter().max().copied().unwrap()
        }).collect::<Vec<_>>();

        let b_inner_width = max(b.iter().map(String::len).max().unwrap(), objective.len());

        // Column counters
        write!(f, "{0:>width$} |", "", width = row_counter_width)?;
        write!(f, " {0:^width$} |", "b", width = b_inner_width)?;
        for (j, width) in column_width.iter().enumerate() {
            write!(f, " {0:^width$}", j, width = width)?;
        }
        writeln!(f)?;

        let total_width = (row_counter_width + 1) + 1 + (1 + b_inner_width + 1) + 1 +
            column_width.iter().map(|l| 1 + l).sum::<usize>();
        // Separator
        writeln!(f, "{}", "-".repeat(total_width))?;

        // Cost row
        write!(f, "{0:>width$} |", "cost", width = row_counter_width)?;
        write!(f, " {0:^width$} |", objective, width = b_inner_width)?;
        for (j, width) in column_width.iter().enumerate() {
            write!(f, " {0:^width$}", cost[j], width = width)?;
        }
        writeln!(f)?;

        // Separator
        writeln!(f, "{}", "-".repeat(total_width))?;

        // Row counter and row data
        for i in 0..self.nr_rows() {
            write!(f, "{0:>width$} |", i, width = row_counter_width)?;
            write!(f, " {0:^width$} |", b[i], width = b_inner_width)?;
            for (j, width) in column_width.iter().enumerate() {
                write!(f, " {0:^width$}", columns[j][i], width = width)?;
            }
            writeln!(f)?;
        }
        writeln!(f)?;

        writeln!(f, "=== Basis Columns ===")?;
        let mut basis = (0..self.nr_rows())
            .map(|i| (i, self.inverse_maintainer.basis_column_index_for_row(i)))
            .collect::<Vec<_>>();
        basis.sort_by_key(|&(i, _)| i);
        writeln!(f, "{:?}", basis)?;

        writeln!(f, "=== Basis Inverse ===")?;
        self.inverse_maintainer.fmt(f)
    }
}

#[cfg(test)]
mod test {
    use std::collections::HashSet;

    use relp_num::{Rational64, RationalBig};
    use relp_num::RB;

    use crate::algorithm::two_phase::matrix_provider::matrix_data::MatrixData;
    use crate::algorithm::two_phase::strategy::pivot_rule::{FirstProfitable, PivotRule};
    use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::basis_inverse_rows::BasisInverseRows;
    use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::Carry;
    use crate::algorithm::two_phase::tableau::inverse_maintenance::ColumnComputationInfo;
    use crate::algorithm::two_phase::tableau::kind::non_artificial::NonArtificial;
    use crate::algorithm::two_phase::tableau::Tableau;
    use crate::data::linear_algebra::vector::{DenseVector, SparseVector};
    use crate::data::linear_algebra::vector::test::TestVector;
    use crate::tests::problem_2::{artificial_tableau_form, create_matrix_data_data, matrix_data_form};

    type T = Rational64;
    type S = RationalBig;

    fn tableau<'a>(
        data: &'a MatrixData<'a, T>,
    ) -> Tableau<Carry<S, BasisInverseRows<S>>, NonArtificial<'a, MatrixData<'a, T>>> {
        let carry = {
            let minus_objective = RB!(-6);
            let minus_pi = DenseVector::from_test_data(vec![1, -1, -1]);
            let b = DenseVector::from_test_data(vec![1, 2, 3]);
            let basis_indices = vec![2, 3, 4];
            let basis_inverse_rows = BasisInverseRows::new(vec![
                SparseVector::from_test_data(vec![1, 0, 0]),
                SparseVector::from_test_data(vec![-1, 1, 0]),
                SparseVector::from_test_data(vec![-1, 0, 1]),
            ]);
            Carry::<S, BasisInverseRows<S>>::new(minus_objective, minus_pi, b, basis_indices, basis_inverse_rows)
        };
        let mut basis_columns = HashSet::new();
        basis_columns.insert(2);
        basis_columns.insert(3);
        basis_columns.insert(4);

        Tableau::<_, NonArtificial<_>>::new_with_inverse_maintainer(
            data,
            carry,
            basis_columns,
        )
    }

    #[test]
    fn cost() {
        let (constraints, b, variables) = create_matrix_data_data();
        let matrix_data_form = matrix_data_form(&constraints, &b, &variables);
        let artificial_tableau = artificial_tableau_form(&matrix_data_form);
        assert_eq!(artificial_tableau.objective_function_value(), RB!(8));

        let tableau = tableau(&matrix_data_form);
        assert_eq!(tableau.objective_function_value(), RB!(6));
    }

    #[test]
    fn relative_cost() {
        let (constraints, b, variables) = create_matrix_data_data();
        let matrix_data_form = matrix_data_form(&constraints, &b, &variables);
        let artificial_tableau = artificial_tableau_form(&matrix_data_form);
        assert_eq!(artificial_tableau.relative_cost(0), RB!(0));

        assert_eq!(
            artificial_tableau.relative_cost(artificial_tableau.nr_artificial_variables() + 0),
            RB!(-10),
        );

        let tableau = tableau(&matrix_data_form);
        assert_eq!(tableau.relative_cost(0), RB!(-3));
        assert_eq!(tableau.relative_cost(1), RB!(-3));
        assert_eq!(tableau.relative_cost(2), RB!(0));
    }

    #[test]
    fn generate_column() {
        let (constraints, b, variables) = create_matrix_data_data();
        let matrix_data_form = matrix_data_form(&constraints, &b, &variables);
        let artificial_tableau = artificial_tableau_form(&matrix_data_form);

        let index_to_test = artificial_tableau.nr_artificial_variables() + 0;
        let column = artificial_tableau.generate_column(index_to_test);
        let expected = SparseVector::from_test_data(vec![3, 5, 2]);
        assert_eq!(column.column(), &expected);
        let result = artificial_tableau.relative_cost(index_to_test);
        assert_eq!(result, RB!(-10));

        let tableau = tableau(&matrix_data_form);
        let index_to_test = 0;
        let column = tableau.generate_column(index_to_test);
        let expected = SparseVector::from_test_data(vec![3, 2, -1]);
        assert_eq!(column.column(), &expected);
        let result = tableau.relative_cost(index_to_test);
        assert_eq!(result, RB!(-3));
    }

    #[test]
    fn bring_into_basis() {
        let (constraints, b, variables) = create_matrix_data_data();
        let matrix_data_form = matrix_data_form(&constraints, &b, &variables);
        let mut artificial_tableau = artificial_tableau_form(&matrix_data_form);
        let column = artificial_tableau.nr_artificial_variables() + 0;
        let column_data = artificial_tableau.generate_column(column);
        let row = artificial_tableau.select_primal_pivot_row(&column_data).unwrap();
        let cost = artificial_tableau.relative_cost(column);
        artificial_tableau.bring_into_basis(column, row, column_data, cost);

        assert!(artificial_tableau.is_in_basis(column));
        assert!(!artificial_tableau.is_in_basis(0));
        assert_eq!(artificial_tableau.objective_function_value(), RB!(14, 3));

        let mut tableau = tableau(&matrix_data_form);
        let column = 1;
        let column_data = tableau.generate_column(column);
        let row = tableau.select_primal_pivot_row(&column_data).unwrap();
        let cost = tableau.relative_cost(column);
        tableau.bring_into_basis(column, row, column_data, cost);

        assert!(tableau.is_in_basis(column));
        assert_eq!(tableau.objective_function_value(), RB!(9, 2));
    }

    fn bfs_tableau<'a>(
        data: &'a MatrixData<'a, Rational64>,
    ) -> Tableau<Carry<RationalBig, BasisInverseRows<RationalBig>>, NonArtificial<'a, MatrixData<'a, Rational64>>> {
        let m = 3;
        let carry = {
            let minus_objective = RB!(0);
            let minus_pi = DenseVector::from_test_data(vec![1, 1, 1]);
            let b = DenseVector::from_test_data(vec![1, 2, 3]);
            let basis_indices = vec![m + 2, m + 3, m + 4];
            let basis_inverse_rows = BasisInverseRows::new(vec![
                SparseVector::from_test_data(vec![1, 0, 0]),
                SparseVector::from_test_data(vec![-1, 1, 0]),
                SparseVector::from_test_data(vec![-1, 0, 1]),
            ]);
            Carry::new(minus_objective, minus_pi, b, basis_indices, basis_inverse_rows)
        };
        let mut basis_columns = HashSet::new();
        basis_columns.insert(m + 2);
        basis_columns.insert(m + 3);
        basis_columns.insert(m + 4);

        Tableau::<_, NonArtificial<_>>::new_with_inverse_maintainer(
            data,
            carry,
            basis_columns,
        )
    }

    #[test]
    fn create_tableau() {
        let (constraints, b, variables) = create_matrix_data_data();
        let matrix_data_form = matrix_data_form(&constraints, &b, &variables);
        let bfs_tableau = bfs_tableau(&matrix_data_form);
        let mut rule = <FirstProfitable as PivotRule>::new();
        assert!(rule.select_primal_pivot_column(&bfs_tableau).is_none());
    }
}
