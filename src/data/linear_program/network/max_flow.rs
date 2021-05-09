//! # Maximum Flow Problem
use std::iter;
use std::ops::{Add, Mul, Range};

use num::Zero;

use crate::algorithm::two_phase::matrix_provider::column::{Column as ColumnTrait, OrderedColumn};
use crate::algorithm::two_phase::matrix_provider::column::identity::IdentityColumn;
use crate::algorithm::two_phase::matrix_provider::filter::generic_wrapper::IntoFilteredColumn;
use crate::algorithm::two_phase::matrix_provider::MatrixProvider;
use crate::algorithm::two_phase::phase_one::PartialInitialBasis;
use crate::algorithm::utilities::remove_sparse_indices;
use crate::data::linear_algebra::matrix::{ColumnMajor, Sparse as SparseMatrix};
use crate::data::linear_algebra::SparseTuple;
use crate::data::linear_algebra::traits::{SparseComparator, SparseElement};
use crate::data::linear_algebra::vector::{DenseVector, SparseVector};
use crate::data::linear_program::elements::BoundDirection;
use crate::data::linear_program::network::representation::{ArcDirection, ArcIncidenceMatrix};
use crate::data::number_types::rational::RationalBig;
use crate::data::number_types::nonzero::Nonzero;

/// Maximum flow problem.
struct Primal<F> {
    /// For each edge, two values indicating from which value the arc leaves, and where it goes to.
    arc_incidence_matrix: ArcIncidenceMatrix,
    capacity: DenseVector<F>,

    /// Source
    s: usize,
    /// Sink
    t: usize,
    /// Column indices which are arcs coming from the source.
    s_arc_range: Range<usize>,
}

enum Cost {
    Zero,
    MinusOne,
}

impl<F> Primal<F>
where
    F: SparseElement<F> + SparseComparator,
{
    pub fn new(
        adjacency_matrix: SparseMatrix<F, F, ColumnMajor>,
        s: usize,
        t: usize,
    ) -> Self {
        let nr_vertices = adjacency_matrix.nr_columns();
        debug_assert!(s < nr_vertices && t < nr_vertices);

        let arcs_before_s = adjacency_matrix.data[..s].iter().map(Vec::len).sum();
        let arcs_leaving_s = adjacency_matrix.data[s].len();
        let s_arc_range = arcs_before_s..(arcs_before_s + arcs_leaving_s);
        let (arc_incidence_matrix, capacity) = ArcIncidenceMatrix::new(adjacency_matrix, vec![s, t]);

        Self {
            arc_incidence_matrix,
            capacity,

            s,
            t,
            s_arc_range,
        }
    }

    pub fn nr_vertices(&self) -> usize {
        self.arc_incidence_matrix.nr_vertices()
    }

    pub fn nr_edges(&self) -> usize {
        self.arc_incidence_matrix.nr_edges()
    }
}

#[derive(Clone)]
struct Column {
    constraint_values: Vec<(usize, ArcDirection)>,
    slack: (usize, ArcDirection),
}

impl ColumnTrait for Column {
    type F = ArcDirection;
    type Iter<'a> = impl Iterator<Item = &'a SparseTuple<Self::F>> + Clone;

    fn iter(&self) -> Self::Iter<'_> {
        self.constraint_values.iter().chain(iter::once(&self.slack))
    }

    fn index_to_string(&self, i: usize) -> String {
        let in_constraint_values = self.constraint_values.iter().find(|&&(ii, _)| ii == i);
        if let Some((_, direction)) = in_constraint_values {
            direction.to_string()
        } else if self.slack.0 == i {
            self.slack.1.to_string()
        } else {
            "0".to_string()
        }
    }
}
impl OrderedColumn for Column {
}
impl IdentityColumn for Column {
    fn identity(i: usize, _len: usize) -> Self {
        Self {
            constraint_values: Vec::with_capacity(0),
            slack: (i, ArcDirection::Incoming)
        }
    }
}
impl IntoFilteredColumn for Column {
    type Filtered = Self;

    fn into_filtered(mut self, to_remove: &[usize]) -> Self::Filtered {
        remove_sparse_indices(&mut self.constraint_values, to_remove);
        // Slack columns are never removed.
        self
    }
}


impl<F> MatrixProvider for Primal<F>
where
    F: SparseElement<F> + Zero + Eq + Nonzero,
{
    type Column = Column;
    type Cost<'a> = Cost;
    type Rhs = F;

    fn column(&self, j: usize) -> Self::Column {
        debug_assert!(j < self.nr_columns());

        if j < self.nr_edges() {
            Column {
                constraint_values: self.arc_incidence_matrix.column(j),
                slack: (self.nr_constraints() + j, ArcDirection::Incoming),
            }
        } else {
            Column {
                constraint_values: Vec::with_capacity(0),
                slack: (self.nr_constraints() + j - self.nr_edges(), ArcDirection::Incoming)
            }
        }
    }

    fn cost_value(&self, j: usize) -> Self::Cost<'_> {
        debug_assert!(j < self.nr_columns());

        if self.s_arc_range.contains(&j) {
            Cost::MinusOne
        } else {
            Cost::Zero
        }
    }

    fn right_hand_side(&self) -> DenseVector<F> {
        let mut b = DenseVector::constant(F::zero(), self.nr_constraints());
        b.extend_with_values(self.capacity.data.clone());
        b
    }

    fn bound_row_index(&self, j: usize, bound_type: BoundDirection) -> Option<usize> {
        debug_assert!(j < self.nr_columns());

        match bound_type {
            BoundDirection::Lower => None,
            BoundDirection::Upper => if j < self.nr_edges() {
                Some(self.nr_constraints() + j)
            } else {
                None
            },
        }
    }

    fn nr_constraints(&self) -> usize {
        self.nr_vertices() - 2
    }

    fn nr_variable_bounds(&self) -> usize {
        self.nr_edges()
    }

    fn nr_columns(&self) -> usize {
        // All edges and their slacks
        self.nr_edges() + self.nr_edges()
    }

    fn reconstruct_solution<G>(&self, column_values: SparseVector<G, G>) -> SparseVector<G, G> {
        unimplemented!()
    }
}

impl<F> PartialInitialBasis for Primal<F>
where
    F: SparseElement<F> + Zero + Eq + Nonzero,
{
    fn pivot_element_indices(&self) -> Vec<(usize, usize)> {
        (0..self.nr_edges()).map(|j| (j + self.nr_constraints(), self.nr_edges() + j)).collect()
    }

    fn nr_initial_elements(&self) -> usize {
        self.nr_edges()
    }
}

impl Add<Cost> for RationalBig {
    type Output = Self;

    fn add(self, rhs: Cost) -> Self::Output {
        match rhs {
            Cost::Zero => self,
            Cost::MinusOne => {
                let (numer, denom): (num::BigInt, num::BigInt) = self.0.into();
                Self(num::BigRational::new(numer - &denom, denom))
            }
        }
    }
}

impl Mul<Cost> for &RationalBig {
    type Output = RationalBig;

    fn mul(self, rhs: Cost) -> Self::Output {
        match rhs {
            Cost::Zero => Self::Output::zero(),
            Cost::MinusOne => -self,
        }
    }
}

impl Mul<Cost> for RationalBig {
    type Output = RationalBig;

    fn mul(self, rhs: Cost) -> Self::Output {
        match rhs {
            Cost::Zero => Self::Output::zero(),
            Cost::MinusOne => -self,
        }
    }
}

#[cfg(test)]
mod test {
    use crate::algorithm::{OptimizationResult, SolveRelaxation};
    use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::basis_inverse_rows::BasisInverseRows;
    use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::Carry;
    use crate::data::linear_algebra::matrix::{ColumnMajor, Order};
    use crate::data::linear_algebra::vector::SparseVector;
    use crate::data::linear_algebra::vector::test::TestVector;
    use crate::data::linear_program::network::max_flow::Primal;
    use crate::data::number_types::rational::{Rational64, RationalBig};

    #[test]
    fn test_1() {
        type T = Rational64;
        type S = RationalBig;

        // Example from Papadimitriou's Combinatorial Optimization.
        let data = ColumnMajor::from_test_data::<T, T, _>(&vec![
            // Directed; from is top, to is on the right
            //   s  a  b  t
            vec![0, 0, 0, 0], // s
            vec![2, 0, 0, 0], // a
            vec![1, 1, 0, 0], // b
            vec![0, 1, 2, 0], // t
        ], 4);
        let problem = Primal::new(data, 0, 3);

        SolveRelaxation::solve_relaxation::<Carry<S, BasisInverseRows<S>>>(&problem);

        debug_assert_eq!(
            problem.solve_relaxation::<Carry<S, BasisInverseRows<S>>>(),
            OptimizationResult::FiniteOptimum(SparseVector::from_test_data(
                vec![2, 1, 1, 1, 2, 0, 0, 0, 0, 0]
            )),
        );
    }
}
