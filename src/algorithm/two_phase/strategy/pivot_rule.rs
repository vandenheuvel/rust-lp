//! # Pivot rules
//!
//! Strategies for moving from basis to basis, whether primal or dual.
use std::ops::Range;

use num_traits::Zero;

use crate::algorithm::two_phase::matrix_provider::column::Column;
use crate::algorithm::two_phase::tableau::inverse_maintenance::{InverseMaintener, ops as im_ops};
use crate::algorithm::two_phase::tableau::kind::Kind;
use crate::algorithm::two_phase::tableau::Tableau;
use crate::data::linear_algebra::SparseTuple;

/// Deciding how to pivot.
///
/// During the Simplex method, one needs to decide how to move from basic solution to basic
/// solution. The pivot rule describes that behavior.
///
/// Once the column has been selected for a primal pivot (or the row for a dual pivot), a row
/// (column) needs to be found. This decision is currently made independent of the strategy.
pub trait PivotRule {
    /// Create a new instance.
    fn new() -> Self;
    /// Column selection rule for the primal Simplex method.
    fn select_primal_pivot_column<IM, K>(
        &mut self,
        tableau: &Tableau<IM, K>,
    ) -> Option<SparseTuple<IM::F>>
    where
        IM: InverseMaintener<F: im_ops::Column<<K::Column as Column>::F> + im_ops::Cost<K::Cost>>,
        K: Kind,
    ;
}

/// Simply pivot on the first column, which has a negative relative cost.
///
/// The behavior is the same for both the primal and dual simplex method.
pub struct FirstProfitable;
impl PivotRule for FirstProfitable {
    fn new() -> Self {
        Self
    }

    fn select_primal_pivot_column<IM, K>(
        &mut self,
        tableau: &Tableau<IM, K>,
    ) -> Option<SparseTuple<IM::F>>
    where
        IM: InverseMaintener<F: im_ops::Column<<K::Column as Column>::F> + im_ops::Cost<K::Cost>>,
        K: Kind,
    {
        // TODO(ENHANCEMENT): For artificial tableaus it's a waste to start at 0
        (0..tableau.nr_columns())
            .filter(|&column| !tableau.is_in_basis(column))
            .map(|column| (column, tableau.relative_cost(column)))
            .find(|(_, cost)| cost < &<IM::F as Zero>::zero())
    }
}

/// Small modification w.r.t. the `FirstProfitable` rule; it starts the search from the last
/// column selected.
pub struct FirstProfitableWithMemory {
    last_selected: Option<usize>,
}
impl PivotRule for FirstProfitableWithMemory {
    fn new() -> Self {
        Self { last_selected: None }
    }

    fn select_primal_pivot_column<IM, K>(
        &mut self,
        tableau: &Tableau<IM, K>,
    ) -> Option<SparseTuple<IM::F>>
    where
        IM: InverseMaintener<F: im_ops::Column<<K::Column as Column>::F> + im_ops::Cost<K::Cost>>,
        K: Kind,
    {
        let find = |to_consider: Range<usize>| to_consider
            .filter(|&column| !tableau.is_in_basis(column))
            .map(|column| (column, tableau.relative_cost(column)))
            .find(|(_, cost)| cost < &<IM::F as Zero>::zero());

        let potential = self.last_selected
            .map_or_else(
                // TODO(ENHANCEMENT): For artificial tableaus it's a waste to start at 0
                || find(0..tableau.nr_columns()),
                // TODO(ENHANCEMENT): For artificial tableaus it's a waste to start at 0
                |last| find(last..tableau.nr_columns()).or_else(|| find(0..last)),
            );

        self.last_selected = potential.as_ref().map(|&(i, _)| i);
        potential
    }
}

/// Simply pivot on the column, which has the most negative relative cost.
pub struct SteepestDescent;
impl PivotRule for SteepestDescent {
    fn new() -> Self {
        Self
    }

    fn select_primal_pivot_column<IM, K>(
        &mut self,
        tableau: &Tableau<IM, K>,
    ) -> Option<SparseTuple<IM::F>>
    where
        IM: InverseMaintener<F: im_ops::Column<<K::Column as Column>::F> + im_ops::Cost<K::Cost>>,
        K: Kind,
    {
        let mut smallest = None;
        // TODO(ENHANCEMENT): For artificial tableaus it's a waste to start at 0
        for (j, cost) in (0..tableau.nr_columns())
            .filter(|&column| !tableau.is_in_basis(column))
            .map(|column| (column, tableau.relative_cost(column)))
            .filter(|(_, cost)| cost < &<IM::F as Zero>::zero()) {
            if let Some((existing_j, existing_cost)) = smallest.as_mut() {
                if &cost < existing_cost {
                    *existing_j = j;
                    *existing_cost = cost;
                }
            } else { smallest = Some((j, cost)) }
        }

        smallest
    }
}

#[cfg(test)]
mod test {
    use crate::algorithm::two_phase::strategy::pivot_rule::{FirstProfitable, PivotRule};
    use crate::data::linear_algebra::vector::SparseVector;
    use crate::data::linear_algebra::vector::test::TestVector;
    use crate::tests::problem_2;

    #[test]
    fn find_profitable_column() {
        let (constraints, b, variables) = problem_2::create_matrix_data_data();
        let matrix_data = problem_2::matrix_data_form(&constraints, &b, &variables);
        let artificial_tableau = problem_2::artificial_tableau_form(&matrix_data);
        let mut rule = <FirstProfitable as PivotRule>::new();
        assert!(matches!(rule.select_primal_pivot_column(&artificial_tableau), Some((3, _))));

        let tableau = problem_2::tableau_form(&matrix_data);
        let mut rule = <FirstProfitable as PivotRule>::new();
        assert_eq!(rule.select_primal_pivot_column(&tableau), None);
    }

    #[test]
    fn find_pivot_row() {
        let (constraints, b, variables) = problem_2::create_matrix_data_data();
        let matrix_data = problem_2::matrix_data_form(&constraints, &b, &variables);
        let artificial_tableau = problem_2::artificial_tableau_form(&matrix_data);
        let column = SparseVector::from_test_data(vec![3, 5, 2]);
        assert_eq!(artificial_tableau.select_primal_pivot_row(&column), Some(0));

        let column = SparseVector::from_test_data(vec![2, 1, 5]);
        assert_eq!(artificial_tableau.select_primal_pivot_row(&column), Some(0));

        let tableau = problem_2::tableau_form(&matrix_data);
        let column = SparseVector::from_test_data(vec![3, 2, -1]);
        assert_eq!(tableau.select_primal_pivot_row(&column), Some(0));

        let column = SparseVector::from_test_data(vec![2, -1, 3]);
        assert_eq!(tableau.select_primal_pivot_row(&column), Some(0));
    }
}
