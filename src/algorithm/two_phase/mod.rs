//! # The Simplex algorithm
//!
//! This module contains all data structures and logic specific to the simplex algorithm. The
//! algorithm is implemented as described in chapters 2 and 4 of Combinatorial Optimization, a book
//! by Christos H. Papadimitriou and Kenneth Steiglitz.
use crate::algorithm::{OptimizationResult, SolveRelaxation};
use crate::algorithm::two_phase::matrix_provider::{Column, MatrixProvider};
use crate::algorithm::two_phase::matrix_provider::filter::generic_wrapper::{IntoFilteredColumn, RemoveRows};
use crate::algorithm::two_phase::phase_one::{FeasibilityComputeTrait, FullInitialBasis, Rank, RankedFeasibilityResult};
use crate::algorithm::two_phase::strategy::pivot_rule::FirstProfitable;
use crate::algorithm::two_phase::tableau::inverse_maintenance::{ExternalOps, InternalOpsHR, InverseMaintenance};
use crate::algorithm::two_phase::tableau::kind::artificial::IdentityColumn;
use crate::algorithm::two_phase::tableau::kind::non_artificial::NonArtificial;
use crate::algorithm::two_phase::tableau::Tableau;

pub(crate) mod phase_one;
pub(crate) mod phase_two;

pub mod tableau;
pub mod matrix_provider;
pub mod strategy;

impl<MP> SolveRelaxation for MP
where
    MP: MatrixProvider<Column: IdentityColumn + IntoFilteredColumn>,
{
    default fn solve_relaxation<IM>(&self) -> OptimizationResult<IM::F>
    where
        IM: InverseMaintenance<F: InternalOpsHR + ExternalOps<<<Self as MatrixProvider>::Column as Column>::F>>,
    {
        // Default choice
        // TODO(ENHANCEMENT): Consider implementing a heuristic to decide these strategies
        //  dynamically
        type NonArtificialPR = FirstProfitable;

        match self.compute_bfs_giving_im::<IM>() {
            RankedFeasibilityResult::Feasible {
                rank,
                nr_artificial_variables,
                basis,
                inverse_maintainer,
            } => match rank {
                Rank::Deficient(rows_to_remove) if !rows_to_remove.is_empty() => {
                    let rows_removed = RemoveRows::new(self, rows_to_remove);
                    let mut non_artificial = Tableau::<_, NonArtificial<_>>::from_artificial_removing_rows(
                        inverse_maintainer,
                        nr_artificial_variables,
                        basis,
                        &rows_removed,
                    );
                    phase_two::primal::<_, _, NonArtificialPR>(&mut non_artificial)
                },
                _ => {
                    let mut non_artificial_tableau = Tableau::<_, NonArtificial<_>>::from_artificial(
                        inverse_maintainer,
                        nr_artificial_variables,
                        basis,
                        self,
                    );
                    phase_two::primal::<_, _, NonArtificialPR>(&mut non_artificial_tableau)
                },
            },
            RankedFeasibilityResult::Infeasible => OptimizationResult::Infeasible,
        }
    }
}

/// Skipping the entire first phase because a basic feasible solution could be cheaply provided.
///
/// TODO(ARCHITECTURE): If this is true, it's ot really two-phase, now is it? Should this solution
///  be provided elsewhere, or the module be renamed?
impl<MP: FullInitialBasis> SolveRelaxation for MP
where
    // TODO: The MP::Column trait bound might not be needed if no rows need to be removed for the
    //  initial basis.
    // TODO: The <MP as MatrixProvider>::Column: IdentityColumn bound is needed because of
    //  limitations of the specialization feature.
    MP: MatrixProvider<Column: IdentityColumn + IntoFilteredColumn>,
{
    fn solve_relaxation<IM>(&self) -> OptimizationResult<IM::F>
    where
        IM: InverseMaintenance<F: InternalOpsHR + ExternalOps<<<Self as MatrixProvider>::Column as Column>::F>, >,
    {
        // TODO(ENHANCEMENT): Consider implementing a heuristic to decide these strategies
        //  dynamically
        type NonArtificialPR = FirstProfitable;

        let basis_indices = self.pivot_element_indices();
        // Sorting of identity matrix columns
        let inverse_maintainer = IM::from_basis_pivots(&basis_indices, self);

        let basis_indices: Vec<usize> = basis_indices.into_iter().map(|(_row, column)| column).collect();
        let basis_columns = basis_indices.iter().copied().collect();
        let mut tableau = Tableau::<_, NonArtificial<_>>::new_with_inverse_maintainer(
            self, inverse_maintainer, basis_indices, basis_columns,
        );
        phase_two::primal::<_, _, NonArtificialPR>(&mut tableau)
    }
}

#[cfg(test)]
mod test {
    use num::FromPrimitive;

    use crate::algorithm::{OptimizationResult, SolveRelaxation};
    use crate::algorithm::two_phase::{phase_one, phase_two, Rank, RankedFeasibilityResult};
    use crate::algorithm::two_phase::matrix_provider::matrix_data::{MatrixData, Variable};
    use crate::algorithm::two_phase::strategy::pivot_rule::FirstProfitable;
    use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::Carry;
    use crate::algorithm::two_phase::tableau::kind::artificial::partially::Partially;
    use crate::algorithm::two_phase::tableau::Tableau;
    use crate::data::linear_algebra::matrix::{ColumnMajor, Order};
    use crate::data::linear_algebra::vector::{Dense as DenseVector, Sparse as SparseVector};
    use crate::data::linear_algebra::vector::test::TestVector;
    use crate::data::linear_program::elements::VariableType;
    use crate::data::number_types::rational::Rational32;
    use crate::R32;
    use crate::tests::problem_2::{create_matrix_data_data, matrix_data_form, tableau_form};

    #[test]
    fn simplex() {
        let (constraints, b) = create_matrix_data_data();
        let matrix_data_form = matrix_data_form(&constraints, &b);
        let mut tableau = tableau_form(&matrix_data_form);
        let result = phase_two::primal::<_, _, FirstProfitable>(&mut tableau);
        assert!(matches!(result, OptimizationResult::FiniteOptimum(_)));
        assert_eq!(tableau.objective_function_value(), R32!(9, 2));
    }

    #[test]
    fn finding_bfs() {
        type T = Rational32;
        // TODO(ENHANCEMENT): Decouple these two types
        type S = T;

        let (constraints, b) = create_matrix_data_data();
        let matrix_data_form = matrix_data_form(&constraints, &b);
        let tableau = Tableau::<Carry<S>, Partially<_>>::new(&matrix_data_form);
        assert!(matches!(
            phase_one::primal::<_, _, MatrixData<T>, FirstProfitable>(tableau),
            RankedFeasibilityResult::Feasible { rank: Rank::Full, .. }
        ));
    }

    #[test]
    fn solve_matrix() {
        type T = Rational32;
        // TODO(ENHANCEMENT): Decouple these two types
        type S = T;

        let (constraints, b) = create_matrix_data_data();
        let matrix_data_form = matrix_data_form(&constraints, &b);

        let result = SolveRelaxation::solve_relaxation::<Carry<S>>(&matrix_data_form);
        //  Optimal value: R64!(4.5)
        assert_eq!(result, OptimizationResult::FiniteOptimum(SparseVector::from_test_tuples(vec![
            (1, 0.5f64),
            (3, 2.5f64),
            (4, 1.5f64),
        ], 5)));
    }

    #[test]
    fn solve_relaxation_1() {
        type T = Rational32;
        // TODO(ENHANCEMENT): Decouple these two types
        type S = T;

        let constraints = ColumnMajor::from_test_data::<T, _, _>(&vec![
            vec![1, 0],
            vec![1, 1],
        ], 2);
        let b = DenseVector::from_test_data(vec![
            3f64 / 2f64,
            5f64 / 2f64,
        ]);
        let variables = vec![
            Variable {
                cost: R32!(-2),
                upper_bound: None,
                variable_type: VariableType::Integer,
            },
            Variable {
                cost: R32!(-1),
                upper_bound: None,
                variable_type: VariableType::Integer,
            },
        ];

        let data = MatrixData::new(&constraints, &b, 0, 2, 0, variables);

        let result = data.solve_relaxation::<Carry<S>>();
        assert_eq!(result, OptimizationResult::FiniteOptimum(SparseVector::from_test_tuples(vec![
            (0, 3f64 / 2f64),
            (1, 1f64),
        ], 4)));
    }
}
