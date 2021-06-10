//! # Phase one: finding a basic feasible solution
use std::collections::HashSet;

use num_traits::Zero;

use crate::algorithm::two_phase::matrix_provider::column::Column;
use crate::algorithm::two_phase::matrix_provider::column::identity::IdentityColumn;
use crate::algorithm::two_phase::matrix_provider::MatrixProvider;
use crate::algorithm::two_phase::strategy::pivot_rule::{FirstProfitableWithMemory, PivotRule};
use crate::algorithm::two_phase::tableau::{debug_assert_in_basic_feasible_solution_state, Tableau};
use crate::algorithm::two_phase::tableau::inverse_maintenance::{ColumnComputationInfo, InverseMaintener, ops as im_ops};
use crate::algorithm::two_phase::tableau::kind::artificial::Artificial;
use crate::algorithm::two_phase::tableau::kind::artificial::Cost;
use crate::algorithm::two_phase::tableau::kind::artificial::fully::Fully as FullyArtificial;
use crate::algorithm::two_phase::tableau::kind::artificial::partially::Partially as PartiallyArtificial;

/// Computing a feasible solution: the first phase of the two phase method.
///
/// This can happen either using the simplex method, or some more specialized method.
///
/// TODO(ENHANCEMENT): Problems that can only have full rank
pub trait FeasibilityComputeTrait: MatrixProvider<Column: IdentityColumn> {
    /// Compute a basic feasible solution.
    ///
    /// # Returns
    ///
    /// A value representing the basic feasible solution, or an indicator that there is none.
    fn compute_bfs_giving_im<IM>(&self) -> RankedFeasibilityResult<IM>
    where
        IM: InverseMaintener<F:
            im_ops::FieldHR +
            im_ops::Column<<<Self as MatrixProvider>::Column as Column>::F> +
            im_ops::Cost<Cost> +
            im_ops::Rhs<<Self as MatrixProvider>::Rhs> +
        >,
    ;
}

/// Most generic implementation: finding a basic feasible solution using the Simplex method.
impl<MP> FeasibilityComputeTrait for MP
where
    MP: MatrixProvider<Column: IdentityColumn>
{
    default fn compute_bfs_giving_im<IM>(&self) -> RankedFeasibilityResult<IM>
    where
        IM: InverseMaintener<F:
            im_ops::FieldHR +
            im_ops::Column<<<Self as MatrixProvider>::Column as Column>::F> +
            im_ops::Cost<Cost> +
            im_ops::Rhs<<Self as MatrixProvider>::Rhs> +
        >,
    {
        // TODO(ENHANCEMENT): Consider implementing a heuristic to decide these strategies
        //  dynamically
        type PivotRule = FirstProfitableWithMemory;

        let artificial_tableau = Tableau::<_, FullyArtificial<_>>::new(self);
        primal::<_, _, MP, PivotRule>(artificial_tableau)
    }
}

/// A few basis columns are already present in the problem.
///
/// Sometimes, a few variables (like positive slack variables) are available that result in less
/// artificial variables being needed.
pub trait PartialInitialBasis: MatrixProvider {
    /// Return the indices of all positive slack variables.
    ///
    /// This is used to find a basic feasible solution faster using the two phase method.
    ///
    /// # Return value
    ///
    /// Collection of tuples with row and column index of (positive) slack coefficients.
    fn pivot_element_indices(&self) -> Vec<(usize, usize)>;

    /// How many positive slacks there are in the problem.
    ///
    /// Is equal to the length of the value returned by `positive_slack_indices`.
    fn nr_initial_elements(&self) -> usize;
}

impl<MP: PartialInitialBasis> FeasibilityComputeTrait for MP
where
    MP: MatrixProvider<Column: IdentityColumn>,
{
    default fn compute_bfs_giving_im<IM>(&self) -> RankedFeasibilityResult<IM>
    where
        IM: InverseMaintener<F:
            im_ops::FieldHR +
            im_ops::Column<<<Self as MatrixProvider>::Column as Column>::F> +
            im_ops::Cost<Cost> +
            im_ops::Rhs<<Self as MatrixProvider>::Rhs> +
        >,
    {
        // TODO(ENHANCEMENT): Consider implementing a heuristic to decide these strategies
        //  dynamically
        type PivotRule = FirstProfitableWithMemory;

        let artificial_tableau = Tableau::<_, PartiallyArtificial<_>>::new(self);
        primal::<_, _, MP, PivotRule>(artificial_tableau)
    }
}

/// A full basis is already present in the problem.
///
/// If the problem is of type Ax <= b, all (positive) slacks can be used as a basic feasible
/// solution (that is, x = 0 is feasible).
///
/// TODO(ENHANCEMENT): Is a marker trait enough, or should there be a specialized method?
/// TODO(ENHANCEMENT): What about the case where also an easy `InverseMaintainer` is available?
pub trait FullInitialBasis: PartialInitialBasis {
}

/// Reduces the artificial cost of the basic feasible solution to zero, if possible. In doing so, a
/// basic feasible solution to the standard form linear program is found.
///
/// # Arguments
///
/// * `tableau`: Artificial tableau with a valid basis. This basis will typically consist of only
/// artificial variables.
///
/// # Return value
///
/// Whether the tableau might allow a basic feasible solution without artificial variables.
pub(crate) fn primal<IM, K, MP, PR>(
    mut tableau: Tableau<IM, K>,
) -> RankedFeasibilityResult<IM>
where
    IM: InverseMaintener<F: im_ops::FieldHR + im_ops::Column<<K::Column as Column>::F> + im_ops::Cost<K::Cost>>,
    K: Artificial,
    MP: MatrixProvider,
    PR: PivotRule,
{
    let mut rule = PR::new();
    loop {
        debug_assert_in_basic_feasible_solution_state(&tableau);

        match rule.select_primal_pivot_column(&tableau) {
            Some((column_nr, cost)) => {
                let column = tableau.generate_column(column_nr);
                match tableau.select_primal_pivot_row(column.column()) {
                    Some(row_nr) => {
                        tableau.bring_into_basis(column_nr, row_nr, column, cost)
                    },
                    None => panic!("Artificial cost can not be unbounded."),
                }
            },
            None => break if tableau.objective_function_value().is_zero() {
                let rank = if tableau.has_artificial_in_basis() {
                    let rows_to_remove = remove_artificial_basis_variables(&mut tableau);
                    if rows_to_remove.is_empty() {
                        Rank::Full
                    } else {
                        Rank::Deficient(rows_to_remove)
                    }
                } else {
                    Rank::Full
                };

                let (im, nr_a, basis) = tableau.into_basis();
                RankedFeasibilityResult::Feasible {
                    rank,
                    nr_artificial_variables: nr_a,
                    inverse_maintainer: im,
                    basis,
                }
            } else {
                RankedFeasibilityResult::Infeasible
            },
        }
    }
}

/// LP's can be either feasible (allowing at least one solution) or infeasible (allowing no
/// solutions).
///
/// If the problem is feasible, it can either have full rank, or be rank deficient.
#[derive(Debug, Eq, PartialEq)]
pub enum RankedFeasibilityResult<IM> {
    /// The problem is feasible and all information necessary to construct a bfs cheaply is in this
    /// variant.
    Feasible {
        /// Whether the problem needs rows to be removed.
        rank: Rank,
        /// The amount of artificial variables that were present in the problem.
        nr_artificial_variables: usize,
        /// The inverse basis that was maintained.
        inverse_maintainer: IM,
        /// Basis indices.
        ///
        /// Note that the second element is just a set version of the elements in the first element.
        basis: HashSet<usize>,
    },
    /// The problem is not feasible.
    Infeasible,
}

/// A matrix or linear program either has full rank, or be rank deficient.
///
/// In case it is rank deficient, a sorted, deduplicated list of (row)indices should be provided,
/// that when removed, makes the matrix or linear program full rank.
#[derive(Debug, Eq, PartialEq)]
pub enum Rank {
    /// The matrix is full rank, no rows need to be removed.
    Full,
    /// The `Vec<usize>` is sorted and contains no duplicate values.
    Deficient(Vec<usize>),
}

/// Removes all artificial variables from the tableau by making a basis change "at zero level", or
/// without change of cost of the current solution.
///
/// # Arguments
///
/// * `tableau`: Tableau to change the basis for.
///
/// # Return value
///
/// A `Vec` with indices of rows that are redundant. Is sorted as a side effect of the algorithm.
///
/// Constraints containing only one variable (that is, they are actually bounds) are read as bounds
/// instead of constraints. All bounds are linearly independent among each other, and with respect
/// to all constraints. As such, they should never be among the redundant rows returned by this
/// method.
fn remove_artificial_basis_variables<IM, K>(
    tableau: &mut Tableau<IM, K>,
) -> Vec<usize>
where
    IM: InverseMaintener<F: im_ops::Column<<K::Column as Column>::F> + im_ops::Cost<K::Cost>>,
    K: Artificial,
{
    let mut artificial_variable_indices = tableau.artificial_basis_columns().into_iter().collect::<Vec<_>>();
    artificial_variable_indices.sort_unstable();
    let mut rows_to_remove = Vec::new();

    for artificial in artificial_variable_indices {
        // The problem was initialized with the artificial variable as a basis column, and it is still in the basis
        debug_assert!(tableau.is_in_basis(artificial));

        let pivot_row = tableau.pivot_row_from_artificial(artificial);
        let column_and_cost = (tableau.nr_artificial_variables()..tableau.nr_columns())
            .filter(|&j| !tableau.is_in_basis(j))
            .map(|j| (j, tableau.relative_cost(j)))
            .filter(|(_, cost)| cost.is_zero())
            .find_map(|(j, cost)| tableau.generate_element(pivot_row, j)
                .map(|_| (j, cost)));

        if let Some((pivot_column, cost)) = column_and_cost {
            let column = tableau.generate_column(pivot_column);
            tableau.bring_into_basis(pivot_column, pivot_row, column, cost);

            debug_assert!(!tableau.is_in_basis(artificial));
        } else {
            rows_to_remove.push(artificial);
        }

        debug_assert_in_basic_feasible_solution_state(&tableau);
    }

    debug_assert!(rows_to_remove.is_sorted());
    rows_to_remove
}

#[cfg(test)]
mod test {
    use relp_num::{Rational64, RationalBig};

    use crate::algorithm::two_phase::matrix_provider::matrix_data::MatrixData;
    use crate::algorithm::two_phase::phase_one::{primal, Rank, RankedFeasibilityResult};
    use crate::algorithm::two_phase::strategy::pivot_rule::FirstProfitable;
    use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::basis_inverse_rows::BasisInverseRows;
    use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::Carry;
    use crate::algorithm::two_phase::tableau::kind::artificial::partially::Partially;
    use crate::algorithm::two_phase::tableau::Tableau;
    use crate::tests::problem_2::{create_matrix_data_data, matrix_data_form};

    #[test]
    fn finding_bfs() {
        type T = Rational64;
        type S = RationalBig;

        let (constraints, b, variables) = create_matrix_data_data();
        let matrix_data_form = matrix_data_form(&constraints, &b, &variables);
        let tableau = Tableau::<Carry<S, BasisInverseRows<S>>, Partially<_>>::new(&matrix_data_form);
        assert!(matches!(
            primal::<_, _, MatrixData<T>, FirstProfitable>(tableau),
            RankedFeasibilityResult::Feasible { rank: Rank::Full, .. }
        ));
    }
}
