//! # Algorithms
use crate::algorithm::two_phase::matrix_provider::column::Column;
use crate::algorithm::two_phase::matrix_provider::MatrixProvider;
use crate::algorithm::two_phase::tableau::inverse_maintenance::InverseMaintener;
use crate::algorithm::two_phase::tableau::inverse_maintenance::ops as im_ops;
use crate::algorithm::two_phase::tableau::kind::artificial::Cost as ArtificialCost;
use crate::data::linear_algebra::vector::SparseVector;

pub mod two_phase;
pub mod criss_cross;
pub mod primal_dual;
pub mod utilities;

/// A problem formulation of which a relaxation can be solved.
/// 
/// Implementations of this trait are specialized according to specific properties. E.g. if a 
/// specific matrix problem has a trivial basic feasible solution, this can be made clear by 
/// implementing a trait that provides this basic feasible solution, and then a specialized
/// implementation of this trait can use that basic feasible solution to bypass the search for one.
pub trait SolveRelaxation: MatrixProvider {
    /// Solve the relaxed version of this problem.
    /// 
    /// In the case of linear programming, that means that integer constraints are ignored.
    /// 
    /// # Return value
    /// 
    /// Whether the problem is feasible, and if so, a solution if the problem is bounded.
    fn solve_relaxation<IM>(&self) -> OptimizationResult<IM::F>
    where
        IM: InverseMaintener<F:
            im_ops::FieldHR +
            im_ops::Column<<Self::Column as Column>::F> +
            im_ops::Cost<ArtificialCost> +
            im_ops::Rhs<Self::Rhs> +
            im_ops::Column<Self::Rhs> +
        >,
        for<'r> IM::F: im_ops::Cost<Self::Cost<'r>>,
    ;
}

/// A linear program is either infeasible, unbounded or has a finite optimum.
///
/// This is determined as the result of an algorithm
#[allow(missing_docs)]
#[derive(Eq, PartialEq, Debug)]
pub enum OptimizationResult<F> {
    Infeasible,
    FiniteOptimum(SparseVector<F, F>),
    Unbounded,
}
