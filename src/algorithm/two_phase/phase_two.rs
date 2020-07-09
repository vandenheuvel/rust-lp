use crate::algorithm::OptimizationResult;
use crate::algorithm::two_phase::matrix_provider::{Column, MatrixProvider};
use crate::algorithm::two_phase::strategy::pivot_rule::PivotRule;
use crate::algorithm::two_phase::tableau::{is_in_basic_feasible_solution_state, Tableau};
use crate::algorithm::two_phase::tableau::inverse_maintenance::{ExternalOps, InternalOpsHR, InverseMaintenance};
use crate::algorithm::two_phase::tableau::kind::non_artificial::NonArtificial;

/// Reduces the cost of the basic feasible solution to the minimum.
///
/// While calling this method, a number of requirements should be satisfied:
/// - There should be a valid basis (not necessarily optimal <=> dual feasible <=> c >= 0)
/// - All constraint values need to be positive (primary feasibility)
///
/// TODO: Write debug tests for these requirements
///
/// # Return value
///
/// An `OptimizationResult` indicating whether or not the problem has a finite optimum. It cannot be
/// infeasible, as a feasible solution is needed to start using this method.
pub(crate) fn primal<'provider, IM, MP, PR>(
    tableau: &mut Tableau<IM, NonArtificial<'provider, MP>>,
) -> OptimizationResult<IM::F>
where
    IM: InverseMaintenance<F: InternalOpsHR + ExternalOps<<MP::Column as Column>::F>>,
// TODO: Is this lifetime bound necessary?
    MP: MatrixProvider + 'provider,
    PR: PivotRule,
{
    let mut rule = PR::new();
    loop {
        debug_assert!(is_in_basic_feasible_solution_state(&tableau));

        match rule.select_primal_pivot_column(tableau) {
            Some((column_index, cost)) => {
                let column = tableau.generate_column(column_index);
                match tableau.select_primal_pivot_row(&column) {
                    Some(row_index) => tableau.bring_into_basis(
                        column_index,
                        row_index,
                        &column,
                        cost,
                    ),
                    None => break OptimizationResult::Unbounded,
                }
            },
            None => break OptimizationResult::FiniteOptimum(tableau.current_bfs()),
        }
    }
}