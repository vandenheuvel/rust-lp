//! # Carry matrix
//!
//! A carry matrix as described in Papadimitriou's book Combinatorial Optimization, but with
//! the ability to use a different basis inverse representation.
use std::cmp::Ordering;
use std::collections::HashSet;
use std::fmt::{Display, Formatter, Debug};
use std::fmt;
use std::ops::Neg;

use relp_num::NonZero;
use relp_num::One;

use crate::algorithm::two_phase::matrix_provider::column::{Column, ColumnNumber, OrderedColumn, SparseColumn};
use crate::algorithm::two_phase::matrix_provider::column::identity::IdentityColumnStruct;
use crate::algorithm::two_phase::matrix_provider::filter::Filtered;
use crate::algorithm::two_phase::matrix_provider::MatrixProvider;
use crate::algorithm::two_phase::tableau::inverse_maintenance::{ColumnComputationInfo, InverseMaintener, ops};
use crate::algorithm::two_phase::tableau::kind::Kind;
use crate::algorithm::utilities::remove_indices;
use crate::data::linear_algebra::SparseTuple;
use crate::data::linear_algebra::traits::Element;
use crate::data::linear_algebra::vector::{DenseVector, SparseVector, Vector};

pub mod basis_inverse_rows;
pub mod lower_upper;

/// The carry matrix represents a basis inverse.
///
/// The carry matrix looks like:
///
///   obj  ||     -pi
/// ----------------------
///    |   ||
///    |   ||    basis
///    b   ||   inverse
///    |   ||     B^-1
///    |   ||
///
/// The dimensions of the matrix are (m + 1) * (m + 1), where m is the number of rows in the
/// problem. Every basis change, this basis inverse matrix changes. As such, the `b` (and
/// `minus_pi`) in this struct don't correspond to the `b` (and `minus_pi`) of the original problem.
///
/// The `b` and `-pi` are stored densely as they are computed with often and are relatively dense.
#[derive(Eq, PartialEq, Clone, Debug)]
pub struct Carry<F, BI> {
    /// Negative of the objective function value.
    ///
    /// Is non-negative (objective value is non-positive) for artificial tableaus.
    minus_objective: F,
    /// Used to compute updated versions of relative cost of columns.
    minus_pi: DenseVector<F>,
    /// Latest version of the constraints.
    pub(super) b: DenseVector<F>,

    /// Maps the rows to the column containing its pivot.
    ///
    /// The rows are indexed 0 through `self.nr_rows()`, while the columns are indexed 0 through
    /// `self.nr_columns()`.
    ///
    /// This attribute changes with a basis change.
    basis_indices: Vec<usize>,

    /// Represents the inverse of the basis columns.
    basis_inverse: BI,
}

/// Facilitating computations with the basis inverse.
pub trait BasisInverse: Display {
    /// Results any computations are yielded in this type.
    type F;
    /// Information gathered while computing a column.
    ///
    /// Introduced to avoid recomputation of the spike during a basis change.
    type ColumnComputationInfo: ColumnComputationInfo<Self::F>;

    /// Create a representation of the identity matrix.
    ///
    /// # Arguments
    ///
    /// * `m`: The number of rows in the problem. The size of the `Carry` matrix that wraps the
    /// implementor is of size `m + 1` by `m + 1`.
    fn identity(m: usize) -> Self;

    /// Invert an ordered collection of columns.
    ///
    /// Note that the implementor can choose to internally permute the columns in order to improve
    /// sparsity.
    fn invert<C: Column + OrderedColumn>(columns: Vec<C>) -> Self
    where
        Self::F: ops::Column<C::F>,
    ;

    /// Update the basis by replacing a basis column.
    ///
    /// # Arguments
    ///
    /// * `pivot_row_index`: The row index of the pivot in the new column `column`. The column that
    /// currently has a `1` at this index (and zeros otherwise) will be removed from the basis.
    /// * `column`: To determine the pivot row index, the column was already explicitly computed
    /// (probably by the same implementor that now accepts this method call). It is now provided
    /// again for insertion into the basis.
    fn change_basis(
        &mut self,
        pivot_row_index: usize,
        column: Self::ColumnComputationInfo,
    );

    /// Compute a column w.r.t. the current basis by computing `B^-1 c`.
    ///
    /// # Arguments
    ///
    /// * `original_column`: Column produced by a `MatrixProvider` instance. Currently needs to
    /// yield elements in order.
    ///
    /// # Return value
    ///
    /// A `SparseVector<T>` of length `m`.
    /// TODO(ENHANCEMENT): Drop the `OrderedColumn` trait bound once it is possible to specialize on
    ///  it, some implementations don't need it.
    fn generate_column<C: Column + OrderedColumn>(
        &self,
        original_column: C,
    ) -> Self::ColumnComputationInfo
    where
        Self::F: ops::Column<C::F>,
    ;

    /// Generate a single element in the tableau with respect to the current basis.
    ///
    /// # Arguments
    ///
    /// * `i`: Row index
    /// * `original_column`: Column with respect to the original basis.
    /// TODO(ENHANCEMENT): Drop the `OrderedColumn` trait bound once it is possible to specialize on
    ///  it.
    fn generate_element<C: Column + OrderedColumn>(
        &self,
        i: usize,
        original_column: C,
    ) -> Option<Self::F>
    where
        Self::F: ops::Column<C::F>,
    ;

    /// Whether this basis inverse should be recomputed.
    ///
    /// A decision rule for when the basis representation has degenerated enough. Meaning, it is
    /// cheaper to recompute than it is to continue updating the existing basis, probably due to
    /// fill-in.
    fn should_refactor(&self) -> bool;

    /// Iterate over a row of the basis inverse matrix.
    fn basis_inverse_row(&self, row: usize) -> SparseVector<Self::F, Self::F>;

    /// Size of the basis who's inverse is represented.
    fn m(&self) -> usize;
}

/// If a basis inverse can more cheaply remove rows from the problem than it would be to recompute
/// complete, this trait can be implemented.
///
/// This operation might be done when converting a basis inverse used to compute a basic feasible
/// solution using artificial variables, that turned out to have redundant rows in them.
pub trait RemoveBasisPart: BasisInverse {
    /// Modify such that it represents the inverse of a modification of the original matrix. That
    /// modification is the removal of some "dimensions" or indices.
    fn remove_basis_part(&mut self, indices: &[usize]);
}

impl<F, BI> Carry<F, BI>
where
    F: ops::Field + ops::FieldHR,
    BI: BasisInverse<F=F>,
{
    /// Create a `Carry` for a tableau with a known basis inverse.
    ///
    /// # Arguments
    ///
    /// * `b`: Constraint values of the original problem with respect to the provided basis.
    /// * `basis_inverse_rows`: Rows of the matrix B^-1, where B is the matrix of columns of the
    /// basis.
    ///
    /// # Return value
    ///
    /// `Carry` with the provided values and a `minus_pi` zero vector.
    pub fn new(
        minus_objective: F,
        minus_pi: DenseVector<F>,
        b: DenseVector<F>,
        basis_indices: Vec<usize>,
        basis_inverse: BI,
    ) -> Self {
        let m = basis_inverse.m();
        debug_assert_eq!(minus_pi.len(), m);
        debug_assert_eq!(b.len(), m);
        debug_assert_eq!(basis_indices.len(), m);

        Carry {
            minus_objective,
            minus_pi,
            b,
            basis_indices,
            basis_inverse,
        }
    }

    /// Create the `minus_pi` field from an existing basis.
    ///
    /// # Arguments
    ///
    /// * `basis_inverse_rows`: A basis inverse that represents a basic feasible solution.
    /// * `provider`: Matrix provider.
    /// * `basis`: Indices of the basis elements.
    fn create_minus_pi_from_artificial<'a, MP: MatrixProvider>(
        basis_inverse: &BI,
        provider: &'a MP,
        basis: &[usize],
    ) -> DenseVector<F>
    where
        F: ops::Column<<MP::Column as Column>::F> + ops::Cost<MP::Cost<'a>>,
    {
        let m = basis_inverse.m();
        debug_assert_eq!(provider.nr_rows(), m);
        debug_assert_eq!(basis.len(), m);

        let b_inverse_columns = (0..m)
            .map(|i| IdentityColumnStruct((i, One)))
            .map(|column| basis_inverse.generate_column(column))
            .map(BI::ColumnComputationInfo::into_column)
            .map(SparseVector::into_iter);
        let mut b_inverse_rows = vec![Vec::new(); m];
        for (j, column) in b_inverse_columns.enumerate() {
            for (i, v) in column {
                b_inverse_rows[i].push((j, v));
            }
        }

        let mut pi = vec![F::zero(); m];
        for (i, inverse_row) in b_inverse_rows.into_iter().enumerate() {
            for (j, value) in inverse_row {
                pi[j] += value * provider.cost_value(basis[i]);
            }
        }

        let data = pi.into_iter().map(Neg::neg).collect::<Vec<_>>();
        let len = data.len();
        DenseVector::new(data, len)
    }

    /// Create the -pi value from an existing basis.
    ///
    /// # Arguments
    ///
    /// * `provider`: Matrix provider.
    /// * `basis`: Basis indices (elements are already shifted, no compensation for the artificial
    /// variables is needed).
    /// * `b`: Constraint values with respect to this basis.
    fn create_minus_obj_from_artificial<'a, MP: MatrixProvider>(
        provider: &'a MP,
        basis: &[usize],
        b: &DenseVector<F>,
    ) -> F
    where
        F: ops::Column<<MP::Column as Column>::F> + ops::Cost<MP::Cost<'a>>,
    {
        let mut objective = F::zero();
        for row in 0..provider.nr_rows() {
            objective += &b[row] * provider.cost_value(basis[row]);
        }
        -objective
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
    fn update_b(
        &mut self,
        pivot_row_index: usize,
        column: &SparseVector<F, F>,
    ) {
        debug_assert!(pivot_row_index < self.m());
        debug_assert_eq!(Vector::len(column), self.m());

        let pivot_value = column.get(pivot_row_index)
            .expect("Pivot value can't be zero.");

        // First normalize the pivot index
        self.b[pivot_row_index] /= pivot_value;

        // Then add multiples of the resulting value to the other values
        // TODO(ARCHITECTURE): Is there a nicer way to go about this?
        let (b_left, b_right) = self.b.data.split_at_mut(pivot_row_index);
        let (b_middle, b_right) = b_right.split_first_mut().unwrap();

        for (edit_row_index, column_value) in column.iter_values() {
            match edit_row_index.cmp(&pivot_row_index) {
                Ordering::Less => {
                    b_left[*edit_row_index] -= column_value * &*b_middle;
                },
                Ordering::Equal => {},
                Ordering::Greater => {
                    b_right[*edit_row_index - (pivot_row_index + 1)] -= column_value * &*b_middle;
                }
            }
        }
    }

    /// Update `self.minus_pi` and the objective function value by performing a row reduction
    /// operation.
    ///
    /// # Arguments
    ///
    /// * `pivot_row_index`: Index of the pivot row.
    /// * `column_value`: Relative cost value for the pivot column.
    ///
    /// # Note
    ///
    /// This method requires a normalized pivot element.
    fn update_minus_pi_and_obj(&mut self, pivot_row_index: usize, relative_cost: F) {
        let basis_inverse_row = self.basis_inverse.basis_inverse_row(pivot_row_index);
        for (column_index, value) in basis_inverse_row.iter_values() {
            self.minus_pi[*column_index] -= &relative_cost * value;
        }

        self.minus_objective -= relative_cost * &self.b[pivot_row_index];
    }

    /// A property of the dimensions of this matrix.
    ///
    /// A `Carry` is always square matrix of size `m + 1` times `m + 1`. Here, `m` is the
    /// length of the constraint column vector and the width of the `minus_pi` row vector. Also, the
    /// basis inverse submatrix B^-1 has dimension `m` times `m`.
    fn m(&self) -> usize {
        self.b.len()
        //           == self.minus_pi.len()
        //           == self.basis_inverse_rows.len()
        //           == self.basis_inverse_rows[0].len()
        //           == ...
        //           == self.basis_inverse_rows[self.m() - 1].len()
    }
}

impl<F, BI> InverseMaintener for Carry<F, BI>
where
    F: ops::Field + ops::FieldHR,
    BI: BasisInverse<F=F>,
{
    type F = F;
    type ColumnComputationInfo = BI::ColumnComputationInfo;

    fn create_for_fully_artificial<Rhs: Element>(
        b: DenseVector<Rhs>
    ) -> Self
    where
        Self::F: ops::Rhs<Rhs>,
    {
        let m = b.len();

        let mut b_sum = Self::F::zero();
        for v in b.iter_values() {
            b_sum += v;
        }

        Self {
            minus_objective: -b_sum,
            minus_pi: DenseVector::constant(-Self::F::one(), m),
            b: DenseVector::new(b.data.into_iter().map(|v| v.into()).collect(), m),
            // Identity matrix
            basis_indices: (0..m).collect(),
            basis_inverse: BI::identity(m),
        }
    }

    fn create_for_partially_artificial<Rhs: Element>(
        artificial_rows: &[usize],
        free_basis_values: &[(usize, usize)],
        b: DenseVector<Rhs>,
        basis_indices: Vec<usize>,
    ) -> Self
    where
        Self::F: ops::Rhs<Rhs>,
    {
        let m = b.len();
        debug_assert_eq!(artificial_rows.len() + free_basis_values.len(), m);  // Correct sizes
        let merged = artificial_rows.iter().copied()
            .chain(free_basis_values.iter().map(|&(i, _)| i)).collect::<HashSet<_>>();
        debug_assert!(merged.iter().all(|&i| i < m));  // Correct range
        debug_assert_eq!(merged.len(), m);  // Uniqueness

        // Initial value of zero is the value that the objective value has when a feasible solution
        // is reached.
        let mut objective = Self::F::zero();
        for &index in artificial_rows {
            // One because we minimize a simple sum of non-negative artificial variables in the
            // basis.
            objective += &b[index];
        }

        // Only the artificial columns "had a cost to them" before they were added to the basis.
        // Putting elements in the basis also influences the minus_pi field.
        let mut counter = 0;
        let minus_pi_values = (0..m).map(|row| {
            if counter < artificial_rows.len() && artificial_rows[counter] == row {
                counter += 1;
                -Self::F::one()
            } else {
                Self::F::zero()
            }
        }).collect();

        Self {
            minus_objective: -objective,
            minus_pi: DenseVector::new(minus_pi_values, m),
            b: DenseVector::new(b.data.into_iter().map(|v| v.into()).collect(), m),
            // Identity matrix
            basis_indices,
            basis_inverse: BI::identity(m),
        }
    }

    fn from_basis<'a, MP: MatrixProvider>(basis: &[usize], provider: &'a MP) -> Self
    where
        Self::F:
            ops::Column<<MP::Column as Column>::F> +
            ops::Rhs<MP::Rhs> +
            ops::Cost<MP::Cost<'a>> +
            ops::Column<MP::Rhs> +
        ,
        MP::Rhs: 'static,
    {
        let columns = basis.iter().map(|&j| provider.column(j)).collect::<Vec<_>>();
        let basis_inverse = BI::invert(columns);

        let b_data = provider.right_hand_side().data
            .into_iter().enumerate()
            .filter(|(_, v)| v.is_not_zero())
            .collect::<Vec<_>>();
        let b_column = SparseColumn { inner: b_data, };
        let mut b_values = vec![F::zero(); provider.nr_rows()];
        for (i, v) in basis_inverse.generate_column(b_column)
            .into_column().into_iter() {
            b_values[i] = v;
        }
        let b = DenseVector::new(b_values, provider.nr_rows());

        let minus_objective = Carry::<_, BI>::create_minus_obj_from_artificial(provider, basis, &b);
        let minus_pi = Carry::create_minus_pi_from_artificial(&basis_inverse, provider, basis);

        Self {
            minus_objective,
            minus_pi,
            b,
            basis_indices: Vec::from(basis),
            basis_inverse,
        }
    }

    fn from_basis_pivots<'a, MP: MatrixProvider>(
        basis_columns: &[(usize, usize)],
        provider: &'a MP,
    ) -> Self
    where
        Self::F:
            ops::Column<<MP::Column as Column>::F> +
            ops::Rhs<MP::Rhs> +
            ops::Cost<MP::Cost<'a>> +
            ops::Column<MP::Rhs> +
        ,
        MP::Rhs: 'static + ColumnNumber,
    {
        let mut elements = Vec::from(basis_columns);
        elements.sort_by_key(|&(row, _)| row);
        let columns = elements.into_iter().map(|(_, column)| column).collect::<Vec<_>>();
        Self::from_basis(&columns, provider)
    }

    fn from_artificial<'provider, MP: MatrixProvider>(
        mut artificial: Self,
        provider: &'provider MP,
        nr_artificial: usize,
    ) -> Self
    where
        F: ops::Column<<MP::Column as Column>::F> + ops::Cost<MP::Cost<'provider>>,
    {
        debug_assert_eq!(artificial.m(), provider.nr_rows());

        for index in &mut artificial.basis_indices {
            *index -= nr_artificial;
        }

        let minus_pi = Carry::create_minus_pi_from_artificial(
            &artificial.basis_inverse,
            provider,
            &artificial.basis_indices,
        );
        let minus_objective = Carry::<_, BI>::create_minus_obj_from_artificial(
            provider,
            &artificial.basis_indices,
            &artificial.b,
        );

        Self::new(minus_objective, minus_pi, artificial.b, artificial.basis_indices, artificial.basis_inverse)
    }

    default fn from_artificial_remove_rows<'provider, MP: Filtered>(
        mut artificial: Self,
        rows_removed: &'provider MP,
        nr_artificial: usize,
    ) -> Self
    where
        Self::F: ops::Column<<<MP as MatrixProvider>::Column as Column>::F> + ops::Cost<MP::Cost<'provider>>,
    {
        debug_assert_eq!(artificial.basis_indices.len(), rows_removed.nr_rows() + rows_removed.filtered_rows().len());

        remove_indices(&mut artificial.basis_indices, rows_removed.filtered_rows());
        for basis_column in &mut artificial.basis_indices {
            *basis_column -= nr_artificial;
        }

        let basis_columns = artificial.basis_indices.iter()
            .map(|&j| rows_removed.column(j))
            .collect();
        let basis_inverse = BI::invert(basis_columns);

        let minus_pi = Carry::create_minus_pi_from_artificial(
            &basis_inverse,
            rows_removed,
            &artificial.basis_indices,
        );

        artificial.b.remove_indices(rows_removed.filtered_rows());

        let minus_obj = Carry::<_, BI>::create_minus_obj_from_artificial(
            rows_removed,
            &artificial.basis_indices,
            &artificial.b,
        );

        Self::new(minus_obj, minus_pi, artificial.b, artificial.basis_indices, basis_inverse)
    }

    fn change_basis(
        &mut self,
        pivot_row_index: usize,
        pivot_column_index: usize,
        column: Self::ColumnComputationInfo,
        relative_cost: Self::F,
    ) -> usize {
        debug_assert!(pivot_row_index < self.m());
        debug_assert_eq!(column.column().len(), self.m());

        // The order of these calls matters: the first of the two normalizes the pivot row
        self.update_b(pivot_row_index, column.column());

        self.basis_inverse.change_basis(pivot_row_index, column);
        self.update_minus_pi_and_obj(pivot_row_index, relative_cost);

        // Update the indices
        let leaving_column = self.basis_indices[pivot_row_index];
        self.basis_indices[pivot_row_index] = pivot_column_index;

        leaving_column
    }

    fn cost_difference<G, C: Column<F=G> + OrderedColumn>(&self, original_column: &C) -> Self::F
    where
        Self::F: ops::Column<G>,
        G: Display + Debug,
    {
        self.minus_pi.sparse_inner_product::<Self::F, _, _>(original_column.iter())
    }

    fn generate_column<G, C: Column<F=G> + OrderedColumn>(
        &self,
        original_column: C,
    ) -> Self::ColumnComputationInfo
    where
        Self::F: ops::Column<G>,
    {
        self.basis_inverse.generate_column(original_column)
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

        self.basis_inverse.generate_element(i, original_column)
    }

    fn after_basis_change<K: Kind>(&mut self, kind: &K)
    where
        Self::F: ops::Column<<<K as Kind>::Column as Column>::F>,
    {
        // TODO(ENHANCEMENT): Refactoring after the basis change means that the last change was made
        //  is discarded.
        if self.basis_inverse.should_refactor() {
            let columns = self.basis_indices.iter()
                .map(|&j| kind.original_column(j))
                .collect::<Vec<_>>();
            self.basis_inverse = BI::invert(columns);
        }
    }

    fn current_bfs(&self) -> Vec<SparseTuple<Self::F>> {
        let mut tuples = self.b.iter_values()
            .enumerate()
            .map(|(i, v)| (self.basis_column_index_for_row(i), v))
            .filter(|(_, v)| !v.is_zero())
            .map(|(i, v)| (i, v.clone()))
            .collect::<Vec<_>>();
        tuples.sort_by_key(|&(i, _)| i);
        tuples
    }

    fn basis_column_index_for_row(&self, row: usize) -> usize {
        self.basis_indices[row]
    }

    fn b(&self) -> DenseVector<Self::F> {
        // TODO(ARCHITECTURE): Avoid this clone, perhaps also alter trait
        self.b.clone()
    }

    fn get_objective_function_value(&self) -> Self::F {
        -self.minus_objective.clone()
    }

    fn get_constraint_value(&self, i: usize) -> &Self::F {
        &self.b[i]
    }
}

impl<F, BI> InverseMaintener for Carry<F, BI>
where
    F: ops::Field + ops::FieldHR,
    BI: BasisInverse<F=F> + RemoveBasisPart,
{
    fn from_artificial_remove_rows<'provider, MP: Filtered>(
        mut artificial: Self,
        rows_removed: &'provider MP,
        nr_artificial: usize,
    ) -> Self
    where
        Self::F: ops::Column<<<MP as MatrixProvider>::Column as Column>::F> + ops::Cost<MP::Cost<'provider>>,
    {
        debug_assert_eq!(artificial.basis_indices.len(), rows_removed.nr_rows() + rows_removed.filtered_rows().len());

        remove_indices(&mut artificial.basis_indices, rows_removed.filtered_rows());
        debug_assert_eq!(artificial.basis_indices.len(), rows_removed.nr_rows());
        for basis_column in &mut artificial.basis_indices {
            *basis_column -= nr_artificial;
        }

        artificial.basis_inverse.remove_basis_part(rows_removed.filtered_rows());

        let minus_pi = Carry::create_minus_pi_from_artificial(
            &artificial.basis_inverse,
            rows_removed,
            &artificial.basis_indices,
        );

        artificial.b.remove_indices(rows_removed.filtered_rows());

        let minus_objective = Carry::<_, BI>::create_minus_obj_from_artificial(
            rows_removed,
            &artificial.basis_indices,
            &artificial.b,
        );

        Self::new(
            minus_objective,
            minus_pi,
            artificial.b,
            artificial.basis_indices,
            artificial.basis_inverse,
        )
    }
}

impl<F, BI> Display for Carry<F, BI>
where
    F: ops::Field + ops::FieldHR,
    BI: BasisInverse<F=F>,
{
    fn fmt(&self, f: &mut Formatter) -> fmt::Result {
        // writeln!(f, "Carry:\n============")?;
        // writeln!(f, "Objective function value: {}", self.get_objective_function_value())?;
        writeln!(f, "Column ordering:")?;
        writeln!(f, "{:?}", self.basis_indices)?;
        writeln!(f, "Minus PI:")?;
        <DenseVector<F> as Display>::fmt(&self.minus_pi, f)?;
        // writeln!(f, "b:")?;
        // <DenseVector<F> as Display>::fmt(&self.b, f)?;
        writeln!(f, "B^-1:")?;

        self.basis_inverse.fmt(f)?;
        writeln!(f)
    }
}
