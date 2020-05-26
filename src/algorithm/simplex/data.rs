//! # Data structures for Simplex
//!
//! Contains the simplex tableau and logic for elementary operations which can be perfomed upon it.
//! The tableau is extended with supplementary data structures for efficiency.
use std::collections::HashSet;
use std::fmt::{Display, Formatter, Result as FormatResult};
use std::iter::repeat;
use std::marker::PhantomData;

use itertools::repeat_n;

use crate::algorithm::simplex::matrix_provider::{MatrixProvider, VariableFeasibilityLogic, BoundType};
use crate::algorithm::utilities::remove_indices;
use crate::data::linear_algebra::matrix::{ColumnMajorOrdering, MatrixOrder, SparseMatrix};
use crate::data::linear_algebra::vector::{DenseVector, SparseVector, Vector};
use crate::data::number_types::traits::{Field, OrderedField};

/// The most high-level data structure that is used by the Simplex algorithm: the Simplex tableau.
///
/// It holds only a reference to the (immutable) problem it solves, but owns the data structures
/// that describe the current solution basis.
#[derive(Eq, PartialEq, Debug)]
pub struct Tableau<'a, F: Field, TT: TableauType, MP: MatrixProvider<F>> {
    /// Supplies data about the problem.
    ///
    /// This data doesn't change throughout the lifetime of this `Tableau`, and it is independent of
    /// the current basis as described by the `carry` and `basis_columns` attributes.
    provider: &'a MP,

    /// Matrix of size (m + 1) x (m + 1).
    ///
    /// This attribute changes with a basis change.
    carry: CarryMatrix<F>,

    /// Maps the rows to the column containing its pivot.
    ///
    /// The rows are indexed 0 through self.nr_rows(), while the columns are indexed 0 through
    /// self.nr_columns().
    ///
    /// This attribute changes with a basis change.
    basis_indices: Vec<usize>,
    /// All columns currently in the basis.
    ///
    /// Could also be derived from `basis_indices`, but is here for faster reading and writing.
    basis_columns: HashSet<usize>,

    phantom_data: PhantomData<TT>,
}

impl<'a, F: Field, TT: TableauType, MP: MatrixProvider<F>> Tableau<'a, F, TT, MP> {
    /// Creates a Simplex tableau with a specific basis.
    ///
    /// Currently only used for testing.
    ///
    /// # Arguments
    ///
    /// * `provider`: Provides the original problem for which the other arguments describe a basis.
    /// * `carry`: `CarryMatrix` with the basis transformation. Corresponds to `basis_indices`.
    /// * `basis_indices`: Maps each row to a column, describing a basis. Corresponds to `carry`.
    ///
    /// # Return value
    ///
    /// `Tableau` with for the provided problem with the provided basis.
    pub(crate) fn new_with_basis(
        provider: &'a MP,
        carry: CarryMatrix<F>,
        basis_indices: Vec<usize>,
        basis_columns: HashSet<usize>,
    ) -> Self {
        Tableau {
            provider,
            carry,
            basis_indices,
            basis_columns,
            phantom_data: PhantomData,
        }
    }

    /// Brings a column into the basis by updating the `self.carry` matrix and updating the
    /// data structures holding the collection of basis columns.
    pub fn bring_into_basis(
        &mut self,
        pivot_column_index: usize,
        pivot_row_index: usize,
        column: &SparseVector<F>,
        cost: F
    ) {
        debug_assert!(pivot_column_index < self.nr_columns());
        debug_assert!(pivot_row_index < self.nr_rows());

        self.carry.row_reduce_for_basis_change(pivot_row_index, column, cost);
        self.update_basis_indices(pivot_row_index, pivot_column_index);
    }

    /// Update the basis index.
    ///
    /// Removes the index of the variable leaving the basis from the `basis_column_map` attribute,
    /// while inserting the entering variable index.
    ///
    /// # Arguments
    ///
    /// * `pivot_row`: Row index of the pivot, in range 0 until self.nr_rows().
    /// * `pivot_column`: Column index of the pivot, in range 0 until self.nr_columns(). Is not yet
    /// in the basis.
    fn update_basis_indices(&mut self, pivot_row: usize, pivot_column: usize) {
        debug_assert!(pivot_row < self.nr_rows());
        debug_assert!(pivot_column < self.nr_columns());

        let leaving_column = self.basis_indices[pivot_row];
        self.basis_columns.remove(&leaving_column);
        self.basis_indices[pivot_row] = pivot_column;
        self.basis_columns.insert(pivot_column);
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
    pub fn relative_cost(&self, j: usize) -> F {
        debug_assert!(j < self.nr_columns());

        self.initial_cost_value(j) + self.carry.cost_difference(TT::original_column(&self, j))
    }

    /// Coefficient of variable `j` in the objective function.
    ///
    /// # Arguments
    ///
    /// * `j`: Column index of the variable, in range 0 until self.nr_columns().
    ///
    /// # Return value
    ///
    /// The cost of variable `j`.
    fn initial_cost_value(&self, j: usize) -> F {
        TT::initial_cost_value(&self, j)
    }

    /// Column of original problem with respect to the current basis.
    ///
    /// Generate a column of the tableau as it would look like with the current basis by matrix
    /// multiplying the original column and the carry matrix.
    ///
    /// # Arguments
    ///
    /// * `j`: Column index of the variable, in range 0 until self.nr_columns().
    ///
    /// # Return value
    ///
    /// `SparseVector<T>` of size `m`.
    pub fn generate_column(&self, j: usize) -> SparseVector<F> {
        debug_assert!(j < self.nr_columns());

        self.carry.generate_column(TT::original_column(&self, j))
    }

    /// Row of the original problem with respect to the current basis.
    pub fn generate_row(&self, i: usize) -> SparseVector<F> {
        debug_assert!(i < self.nr_rows());

        self.carry.generate_row(i, TT::original_matrix(&self))
    }

    /// Single element with respect to the current basis.
    pub fn generate_element(&self, i: usize, j: usize) -> F {
        debug_assert!(i < self.nr_rows());
        debug_assert!(j < self.nr_columns());

        self.carry.generate_element(i, TT::original_column(&self, j))
    }

    /// The constraint values relative to the current basis.
    pub fn constraint_values(&self) -> &DenseVector<F> {
        &self.carry.b
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
    pub fn is_in_basis(&self, column: &usize) -> bool {
        debug_assert!(*column < self.nr_columns());

        self.basis_columns.contains(column)
    }

    /// Get the current basic feasible solution.
    ///
    /// # Return value
    ///
    /// `SparseVector<T>` of the solution.
    pub fn current_bfs(&self) -> SparseVector<F> {
        let mut tuples = self.carry.b.iter_values()
            .enumerate()
            .map(|(i, &v)| (self.basis_indices[i], v))
            .filter(|&(_, v)| v != F::additive_identity())
            .collect::<Vec<_>>();
        tuples.sort_by_key(|&(i, _)| i);
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
    /// This function works for both artificial and non-artificial tableau's.
    pub fn objective_function_value(&self) -> F {
        self.carry.get_objective_function_value()
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
        TT::nr_columns(&self)
    }

    /// Number of rows in the tableau.
    ///
    /// # Return value
    ///
    /// The number of rows.
    pub fn nr_rows(&self) -> usize {
        self.provider.nr_rows()
    }
}

impl<'a, OF: OrderedField, TT: TableauType, MP: MatrixProvider<OF>> Tableau<'a, OF, TT, MP> {
    /// Determine the row to pivot on.
    ///
    /// Determine the row to pivot on, given the column. This is the row with the positive but
    /// minimal ratio between the current constraint vector and the column.
    ///
    /// When there are multiple choices for the pivot row, Bland's anti cycling algorithm
    /// is used to avoid cycles.
    ///
    /// TODO: Reconsider the below
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
    pub fn select_primal_pivot_row(&self, column: &SparseVector<OF>) -> Option<usize> {
        debug_assert_eq!(column.len(), self.nr_rows());

        // (chosen index, minimum ratio, corresponding leaving_column (for Bland's algorithm))
        let mut min_values: Option<(usize, OF, usize)> = None;
        for &(row, xij) in column.iter_values() {
            if xij > OF::additive_identity() {
                let ratio = self.carry.get_constraint_value(row) / xij;
                // Bland's anti cycling algorithm
                let leaving_column = self.basis_indices[row];
                if let Some((min_index, min_ratio, min_leaving_column)) = &mut min_values {
                    if ratio == *min_ratio && leaving_column < *min_leaving_column {
                        *min_index = row;
                        *min_leaving_column = leaving_column;
                    } else if ratio < *min_ratio {
                        *min_index = row;
                        *min_ratio = ratio;
                        *min_leaving_column = leaving_column;
                    }
                } else {
                    min_values = Some((row, ratio, leaving_column))
                }
            }
        }

        min_values.map(|(min_index, _, _)| min_index)
    }

    /// Select the pivot column for a dual-simplex step.
    ///
    /// There is far less room for smart strategies and heuristics in this selection than for the
    /// dual pivot row selection. As such, this function is not included in the `PivotRule` trait.
    ///
    /// # Arguments
    ///
    /// * `row`: The chosen pivot row, on which the pivot element is to be selected.
    ///
    /// # Return value
    ///
    /// If a suitable value is found, the column index is returned. If not, a `None` is returned. In
    /// this case, the problem is (primal) infeasible.
    pub fn select_dual_pivot_column(&self, row: usize) -> Option<usize> {
        let row_data = self.generate_row(row);
        let row_ratio = row_data
            .iter_values()
            .filter(|&&(_, xij)| xij < OF::additive_identity())
            .map(|&(column, xij)| (column, self.relative_cost(column) / xij));

        // (chosen column index, maximum ratio)
        let mut max_values: Option<(usize, OF)> = None;
        for (column, ratio) in row_ratio {
            if let Some((chosen_column_index, max_ratio)) = &mut max_values {
                if ratio > *max_ratio {
                    *chosen_column_index = column;
                    *max_ratio = ratio;
                }
            } else {
                max_values = Some((column, ratio));
            }
        }

        max_values.map(|(column, _)| column)
    }

    /// Check whether the tableau currently has a valid basic feasible solution.
    ///
    /// Only used for debug_purposes.
    pub fn debug_assert_is_in_basic_feasible_solution_state(&self) -> bool {
        // Checking basis_columns
        // Correct number of basis columns (uniqueness is implied because it's a set)
        debug_assert_eq!(self.basis_columns.len(), self.provider.nr_rows());

        // Checking basis_indices
        // Correct number of basis columns
        debug_assert_eq!(self.basis_indices.len(), self.provider.nr_rows());
        // Uniqueness of the basis columns
        debug_assert_eq!(self.basis_indices.iter().collect::<HashSet<_>>().len(), self.provider.nr_rows());
        // Same columns as in `basis_indices`
        let as_set = self.basis_indices.iter().map(|&v| v).collect::<HashSet<_>>();
        debug_assert_eq!(as_set.difference(&self.basis_columns).collect::<Vec<_>>().len(), 0);

        // Checking carry matrix
        // `basis_inverse_rows` are a proper inverse by regenerating basis columns
        for (i, &basis_column) in self.basis_indices.iter().enumerate() {
            let generated = self.generate_column(basis_column);
            debug_assert_eq!(generated, SparseVector::standard_basis_vector(i, self.provider.nr_rows()));
        }
        // `minus_pi` get to relative zero cost for basis columns
        for &basis_column in self.basis_indices.iter() {
            debug_assert_eq!(self.relative_cost(basis_column), OF::additive_identity());
        }
        // `b` >= 0
        let b = self.constraint_values();
        for row in 0..self.nr_rows() {
            debug_assert!(b[row] >= OF::additive_identity());
        }

        true
    }

    /// A simple getter for the internal provider.
    ///
    /// # Return value
    ///
    /// A borrow to the internal provider.
    pub fn provider(&self) -> &MP {
        self.provider
    }
}

/// The tableau type provides two different ways for the `Tableau` to function, depending on whether
/// any virtual artificial variables should be included in the problem.
pub trait TableauType: Sized + Eq + PartialEq {
    /// Coefficient of variable `j` in the objective function.
    ///
    /// # Arguments
    ///
    /// * `j`: Column index of the variable, in range 0 until self.nr_columns().
    ///
    /// # Return value
    ///
    /// The cost of variable `j`.11
    fn initial_cost_value<F: Field, MP: MatrixProvider<F>>(
        tableau: &Tableau<F, Self, MP>,
        j: usize,
    ) -> F;

    /// Get the column from the original problem.
    ///
    /// Depending on whether the tableau is artificial or not, this requires either an artificial
    /// basis column, or a column from the original problem.
    fn original_column<F: Field, MP: MatrixProvider<F>>(
        tableau: &Tableau<F, Self, MP>,
        j: usize,
    ) -> SparseVector<F>;

    /// Tableau as a matrix with respect to the original identity basis.
    fn original_matrix<F: Field, MP: MatrixProvider<F>>(
        tableau: &Tableau<F, Self, MP>,
    ) -> SparseMatrix<F, ColumnMajorOrdering>;

    /// Number of variables in the problem.
    ///
    /// # Return value
    ///
    /// The number of variables.
    ///
    /// # Note
    ///
    /// This number might be extremely large, depending on the `AdaptedTableauProvider`.
    fn nr_columns<F: Field, MP: MatrixProvider<F>>(
        tableau: &Tableau<F, Self, MP>,
    ) -> usize;
}

/// The `TableauType` in case the `Tableau` contains artificial variables.
#[derive(Eq, PartialEq, Debug)]
pub struct Artificial;
impl TableauType for Artificial {
    /// Coefficient of variable `j` in the objective function.
    ///
    /// # Arguments
    ///
    /// * `j`: Column index of the variable, in range 0 until self.nr_columns().
    ///
    /// # Return value
    ///
    /// The cost of variable `j`.
    fn initial_cost_value<F: Field, MP: MatrixProvider<F>>(
        tableau: &Tableau<F, Self, MP>,
        j: usize,
    ) -> F {
        debug_assert!(j < tableau.nr_columns());

        if j < tableau.nr_artificial_variables() {
            F::multiplicative_identity()
        } else {
            F::additive_identity()
        }
    }

    /// Retrieve an original column.
    ///
    /// # Arguments
    ///
    /// * `tableau`: Tableau to retrieve the column from.
    /// * `j`: Column index.
    ///
    /// # Return value
    ///
    /// The generated column, relative to the basis represented in the `Tableau`.
    fn original_column<F: Field, MP: MatrixProvider<F>>(
        tableau: &Tableau<F, Self, MP>,
        j: usize,
    ) -> SparseVector<F> {
        debug_assert!(j < tableau.nr_columns());

        if j < tableau.nr_artificial_variables() {
            SparseVector::standard_basis_vector(j, tableau.nr_rows())
        } else {
            tableau.provider.column(j - tableau.nr_artificial_variables())
        }
    }

    fn original_matrix<F: Field, MP: MatrixProvider<F>>(
        tableau: &Tableau<F, Self, MP>,
    ) -> SparseMatrix<F, ColumnMajorOrdering> {
        ColumnMajorOrdering::identity(tableau.nr_rows()).concatenate_horizontally(tableau.provider.as_sparse_matrix())
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
    fn nr_columns<F: Field, MP: MatrixProvider<F>>(
        tableau: &Tableau<F, Self, MP>,
    ) -> usize {
        tableau.nr_artificial_variables() + tableau.provider.nr_columns()
    }
}

impl<'a, F: Field, MP: MatrixProvider<F>> Tableau<'a, F, Artificial, MP> {
    /// Create a `Tableau` augmented with artificial variables.
    ///
    /// The tableau is then in a basic feasible solution having only the artificial variables in the
    /// basis.
    ///
    /// # Arguments
    ///
    /// * `provider`: Provides the problem to find a basic feasible solution for.
    ///
    /// # Return value
    ///
    /// The tableau.
    pub(crate) fn new(provider: &'a MP) -> Self {
        let m = provider.nr_rows();
        let carry = CarryMatrix::create_for_artificial(provider.constraint_values());

        Tableau {
            provider,
            carry,
            basis_indices: (0..m).collect(),
            basis_columns: (0..m).collect(),
            phantom_data: PhantomData,
        }
    }

    /// Number of artificial variables in this tableau.
    pub fn nr_artificial_variables(&self) -> usize {
        self.provider.nr_rows()
    }

    /// Whether there are any artificial variables in the basis.
    pub fn has_artificial_in_basis(&self) -> bool {
        self.basis_columns.iter().any(|&c| c < self.nr_artificial_variables())
    }

    /// Get the indices of the artificial variables that are still in the basis.
    pub fn artificial_basis_columns(&self) -> HashSet<usize> {
        self.basis_columns
            .iter()
            .filter(|&&v| v < self.nr_artificial_variables())
            .map(|&v| v)
            .collect()
    }
}

/// The `TableauType` in case the `Tableau` does not contain any artificial variables.
///
/// This `Tableau` variant should only be constructed with a known feasible basis.
#[derive(Eq, PartialEq, Debug)]
pub struct NonArtificial;
impl TableauType for NonArtificial {
    /// Coefficient of variable `j` in the objective function.
    ///
    /// # Arguments
    ///
    /// * `tableau`: Tableau to retrieve the cost value from.
    /// * `j`: Column index of the variable, in range 0 until self.nr_columns().
    ///
    /// # Return value
    ///
    /// The cost of variable `j`.
    fn initial_cost_value<F: Field, MP: MatrixProvider<F>>(
        tableau: &Tableau<F, NonArtificial, MP>,
        j: usize,
    ) -> F {
        tableau.provider.cost_value(j)
    }

    /// Retrieve an original column.
    ///
    /// # Arguments
    ///
    /// * `tableau`: Tableau to retrieve the column from.
    /// * `j`: Column index.
    ///
    /// # Return value
    ///
    /// The generated column, relative to the basis represented in the `Tableau`.
    fn original_column<F: Field, MP: MatrixProvider<F>>(
        tableau: &Tableau<F, Self, MP>,
        j: usize,
    ) -> SparseVector<F> {
        debug_assert!(j < tableau.nr_columns());

        tableau.provider.column(j)
    }

    fn original_matrix<F: Field, MP: MatrixProvider<F>>(
        tableau: &Tableau<F, Self, MP>,
    ) -> SparseMatrix<F, ColumnMajorOrdering> {
        tableau.provider.as_sparse_matrix()
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
    fn nr_columns<F: Field, MP: MatrixProvider<F>>(
        tableau: &Tableau<F, NonArtificial, MP>,
    ) -> usize {
        tableau.provider.nr_columns()
    }
}

impl<'a, F: Field, MP: MatrixProvider<F>> Tableau<'a, F, NonArtificial, MP> {
    /// Create a `Tableau` from an artificial tableau.
    ///
    /// # Arguments
    ///
    /// * `artificial_tableau`: `Tableau` instance created with artificial variables.
    ///
    /// # Return value
    ///
    /// `Tableau` with the same basis, but non-artificial cost row.
    pub fn from_artificial(
        artificial_tableau: Tableau<'a, F, Artificial, MP>,
    ) -> Self {
        debug_assert!(artificial_tableau.artificial_basis_columns().is_empty());

        let nr_artificial = artificial_tableau.nr_artificial_variables();
        let basis_indices = artificial_tableau.basis_indices.into_iter()
            .map(|column| column - nr_artificial)
            .collect();
        let basis_columns = artificial_tableau.basis_columns
            .into_iter()
            .map(|column| column - nr_artificial)
            .collect();
        let carry = CarryMatrix::from_artificial(
            artificial_tableau.carry,
            artificial_tableau.provider,
            &basis_indices,
        );

        Tableau {
            provider: artificial_tableau.provider,
            carry,
            basis_indices,
            basis_columns,
            phantom_data: PhantomData,
        }
    }

    /// Create a `Tableau` from an artificial tableau while removing some rows.
    ///
    /// # Arguments
    ///
    /// * `artificial_tableau`: `Tableau` instance created with artificial variables.
    /// * `rows_removed`: `RemoveRows` instance containing all the **sorted** rows that should be
    /// removed from e.g. the basis inverse matrix.
    ///
    /// # Return value
    ///
    /// `Tableau` with the same basis, but non-artificial cost row.
    pub fn from_artificial_removing_rows<'b: 'a>(
        artificial_tableau: Tableau<'a, F, Artificial, MP>,
        rows_removed: &'b RemoveRows<'a, F, MP>,
    ) -> Tableau<'a, F, NonArtificial, RemoveRows<'a, F, MP>> {
        debug_assert!(
            artificial_tableau.basis_indices.iter()
                .all(|&v| v >= artificial_tableau.nr_artificial_variables() || rows_removed.rows_to_skip.contains(&v))
        );

        let nr_artificial = artificial_tableau.nr_artificial_variables();

        let mut basis_indices = artificial_tableau.basis_indices;
        remove_indices(&mut basis_indices, &rows_removed.rows_to_skip);
        basis_indices.iter_mut().for_each(|index| *index -= nr_artificial);

        // Remove same row and column from carry matrix
        let carry = CarryMatrix::from_artificial_remove_rows(
            artificial_tableau.carry,
            rows_removed,
            &basis_indices,
        );

        let mut basis_columns = artificial_tableau.basis_columns;
        for index in &rows_removed.rows_to_skip {
            let was_there = basis_columns.remove(index);
            debug_assert!(was_there);
        }
        let basis_columns = basis_columns.into_iter().map(|v| v - nr_artificial).collect();


        Tableau {
            provider: rows_removed,
            carry,
            basis_indices,
            basis_columns,
            phantom_data: PhantomData,
        }
    }
}

impl<'a, F, TT, MP> Display for Tableau<'a, F, TT, MP>
    where F: Field, TT: TableauType, MP: MatrixProvider<F> {
    fn fmt(&self, f: &mut Formatter) -> FormatResult {
        writeln!(f, "Tableau:")?;
        writeln!(f, "=== Current Tableau ===")?;
        let column_width = 10;
        let counter_width = 5;
        // Column counter
        write!(f, "{0:width$}", "", width = counter_width)?;
        for column_index in 0..self.nr_columns() {
            write!(f, "{0:>width$}", column_index, width = column_width)?;
        }
        writeln!(f, "")?;

        // Separator
        writeln!(f, "{}", repeat_n("-", counter_width + self.nr_columns() * column_width)
            .collect::<String>())?;

        // Cost row
        write!(f, "{0:width$}", "cost", width = counter_width)?;
        for column_index in 0..self.nr_columns() {
            let number = format!("{}", self.relative_cost(column_index));
            write!(f, "{0:>width$.5}", number, width = column_width)?;
        }
        writeln!(f, "")?;

        // Separator
        writeln!(f, "{}", repeat("-")
            .take(counter_width + self.nr_columns() * column_width)
            .collect::<String>())?;

        // Row counter and row data
        for row_index in 0..self.nr_rows() {
            write!(f, "{0: <width$}", row_index, width = counter_width)?;
            for column_index in 0..self.nr_columns() {
                let number = format!("{}", self.generate_column(column_index)[row_index]);
                write!(f, "{0:>width$.5}", number, width = column_width)?;
            }
            writeln!(f, "")?;
        }
        writeln!(f, "")?;

        writeln!(f, "=== Basis Columns ===")?;
        let mut basis = self.basis_indices.iter()
            .enumerate()
            .map(|(i, &j)| (i, j))
            .collect::<Vec<_>>();
        basis.sort_by_key(|&(i, _)| i);
        writeln!(f, "{:?}", basis)?;

        writeln!(f, "=== Basis Inverse ===")?;
        <CarryMatrix<F> as Display>::fmt(&self.carry, f)

        // TODO
//        writeln!(f, "=== Data Provider ===")?;
//        self.provider.fmt(f)
    }
}

/// The carry matrix looks like:
///
/// obj | minus_pi
/// --------------
///  |  |
///  |  |  basis
///  b  | inverse
///  |  |   B^-1
///  |  |
///
/// The b (and minus_pi) in this struct change continuously; they don't correspond to the b (and
/// minus_pi) of the original problem.
///
/// Used in simplex method, inside the loop, in the methods.
///
/// TODO: Write better docstring
#[derive(Eq, PartialEq, Clone, Debug)]
pub(crate) struct CarryMatrix<T> {
    minus_objective: T,
    minus_pi: DenseVector<T>,
    b: DenseVector<T>,
    basis_inverse_rows: Vec<SparseVector<T>>,
    number_of_basis_changes: u64,
}

impl<F: Field> CarryMatrix<F> {
    /// Create a `CarryMatrix` for a tableau with artificial variables.
    ///
    /// # Arguments
    ///
    /// * `b`: Constraint values of the original problem, i.e. the original `b` with respect to the
    /// unit basis.
    ///
    /// # Return value
    ///
    /// `CarryMatrix` with a `minus_pi` equal to -1's and the standard basis.
    fn create_for_artificial(b: DenseVector<F>) -> Self {
        let mut objective = F::additive_identity();
        for &value in b.iter_values() {
            objective += value;
        }

        let m = b.len();
        CarryMatrix {
            minus_objective: -objective,
            minus_pi: DenseVector::constant(-F::multiplicative_identity(), m),
            b,
            basis_inverse_rows: (0..m).map(|i| SparseVector::standard_basis_vector(i, m)).collect(),
            number_of_basis_changes: 0,
        }
    }

    /// Create a `CarryMatrix` for a tableau with a known basis inverse.
    ///
    /// # Arguments
    ///
    /// * `b`: Constraint values of the original problem with respect to the provided basis.
    /// * `basis_inverse_rows`: Rows of the matrix B^-1, where B is the matrix of columns of the
    /// basis.
    ///
    /// # Return value
    ///
    /// `CarryMatrix` with the provided values and a `minus_pi` zero vector.
    pub fn create(
        minus_objective: F,
        minus_pi: DenseVector<F>,
        b: DenseVector<F>,
        basis_inverse_rows: Vec<SparseVector<F>>,
    ) -> Self {
        let m = minus_pi.len();
        debug_assert_eq!(minus_pi.len(), m);
        debug_assert_eq!(b.len(), m);
        debug_assert!(basis_inverse_rows.iter().all(|r| r.len() == m));

        CarryMatrix {
            minus_objective,
            minus_pi,
            b,
            basis_inverse_rows,
            number_of_basis_changes: 0,
        }
    }

    /// Create a `CarryMatrix` from a carry matrix with an artificial cost row.
    ///
    /// # Arguments
    ///
    /// * `artificial_carry_matrix`: A `CarryMatrix` which has an artificial cost row.
    /// * `provider`: Provider of the problem, used for original cost in `minus_pi` restoration.
    /// * `basis`: Vector with, ordered by pivot row, the basis indices of the new basis.
    ///
    /// # Return value
    ///
    /// `CarryMatrix` with a restored `minus_pi` vector and objective value.
    pub fn from_artificial<MP: MatrixProvider<F>>(
        artificial_carry_matrix: Self,
        provider: &MP,
        basis: &Vec<usize>,
    ) -> Self {
        let minus_pi = CarryMatrix::create_minus_pi_from_artificial(
            &artificial_carry_matrix.basis_inverse_rows,
            provider,
            basis,
        );
        let minus_obj = CarryMatrix::create_minus_obj_from_artificial(
            provider,
            basis,
            &artificial_carry_matrix.b,
        );

        Self::create(
            minus_obj,
             minus_pi,
             artificial_carry_matrix.b,
             artificial_carry_matrix.basis_inverse_rows,
        )
    }

    /// Create a `CarryMatrix` from a carry matrix with an artificial cost row while removing
    /// several row and column indices.
    ///
    /// # Arguments
    ///
    /// * `artificial_carry_matrix`: A `CarryMatrix` which has an artificial cost row.
    /// * `provider`: Provider of the reduced problem.
    /// * `basis`: Vector with, ordered by pivot row, the basis indices of the new basis.
    /// * `to_remove`: All rows and columns to remove from the carry matrix.
    ///
    /// # Return value
    ///
    /// `CarryMatrix` of reduced dimension with a restored `minus_pi` vector and objective value.
    pub fn from_artificial_remove_rows<MP: MatrixProvider<F>>(
        artificial_carry_matrix: Self,
        remove_rows_provider: &RemoveRows<F, MP>,
        basis: &Vec<usize>,
    ) -> Self {
        debug_assert_eq!(basis.len(), remove_rows_provider.nr_rows());

        // Remove the rows
        let mut basis_inverse_rows = artificial_carry_matrix.basis_inverse_rows;
        remove_indices(&mut basis_inverse_rows, &remove_rows_provider.rows_to_skip);
        // Remove the columns
        for element in basis_inverse_rows.iter_mut() {
            element.remove_indices(&remove_rows_provider.rows_to_skip);
        }

        let minus_pi = CarryMatrix::create_minus_pi_from_artificial(
            &basis_inverse_rows,
            remove_rows_provider,
            basis,
        );

        let mut b = artificial_carry_matrix.b;
        b.remove_indices(&remove_rows_provider.rows_to_skip);

        let minus_obj = CarryMatrix::create_minus_obj_from_artificial(
            remove_rows_provider,
            basis,
            &b,
        );

        Self::create(
            minus_obj,
            minus_pi,
            b,
            basis_inverse_rows,
        )
    }

    fn create_minus_pi_from_artificial(
        basis_inverse_rows: &Vec<SparseVector<F>>,
        provider: &impl MatrixProvider<F>,
        basis: &Vec<usize>,
    ) -> DenseVector<F> {
        let m = basis_inverse_rows.len();
        debug_assert_eq!(provider.nr_rows(), m);
        debug_assert_eq!(basis.len(), m);

        let mut pi = repeat_n(F::additive_identity(), m).collect::<Vec<_>>();
        for row in 0..m {
            for &(column, value) in basis_inverse_rows[row].iter_values() {
                pi[column] += provider.cost_value(basis[row]) * value;
            }
        }

        let data = pi.into_iter().map(|v| -v).collect::<Vec<_>>();
        let len = data.len();
        DenseVector::new(data, len)
    }

    fn create_minus_obj_from_artificial<MP: MatrixProvider<F>>(
        provider: &MP,
        basis: &Vec<usize>,
        b: &DenseVector<F>,
    ) -> F {
        let mut objective = F::additive_identity();
        for row in 0..provider.nr_rows() {
            objective += provider.cost_value(basis[row]) * b[row];
        }
        -objective
    }

    /// Update the basis by doing row reduction operations.
    ///
    /// Supply a column, that should become part of the basis, to this function. Then this
    /// `CarryMatrix` gets updated such that if one were to matrix-multiply the concatenation column
    /// `[cost, c for c in column]` with this `CarryMatrix`, an (m + 1)-dimensional unitvector would
    /// be the result.
    ///
    /// # Arguments
    ///
    /// * `pivot_row_index`: Index of the pivot row.
    /// * `column`: Column relative to the current basis to be entered into that basis.
    /// * `cost`: Relative cost of that column. The objective function value will change by this
    /// amount.
    fn row_reduce_for_basis_change(
        &mut self,
        pivot_row_index: usize,
        column: &SparseVector<F>,
        cost: F,
    ) {
        debug_assert!(pivot_row_index < self.m());
        debug_assert_eq!(column.len(), self.m());

        // The order of these calls matters: the first of the two normalizes the pivot row

        self.normalize_pivot_row(pivot_row_index, column);
        self.row_reduce_update_basis_inverse_and_b(pivot_row_index, column);
        self.row_reduce_update_minus_pi_and_obj(pivot_row_index, cost);
        self.number_of_basis_changes += 1;
    }

    /// Normalize the pivot row.
    ///
    /// That is, the pivot value will be set to `1`.
    ///
    /// # Arguments
    ///
    /// * `pivot_row_index`: Index of the pivot row.
    /// * `column`: Column relative to the current basis to be entered into that basis.
    fn normalize_pivot_row(&mut self, pivot_row_index: usize, column: &SparseVector<F>) {
        let pivot_value = column[pivot_row_index];
        debug_assert_ne!(pivot_value, F::additive_identity());
        self.basis_inverse_rows[pivot_row_index].element_wise_divide(pivot_value);
        self.b[pivot_row_index] /= pivot_value;
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
    fn row_reduce_update_basis_inverse_and_b(
        &mut self,
        pivot_row_index: usize,
        column: &SparseVector<F>,
    ) {
        debug_assert!(pivot_row_index < self.m());
        debug_assert_eq!(column.len(), self.m());

        // Do the other non-minus_pi rows
        for &(edit_row_index, column_value) in column.iter_values() {
            if edit_row_index != pivot_row_index {
                // Update the constraint
                let pivot_row_b = self.b[pivot_row_index];
                self.b.shift_value(edit_row_index,-column_value * pivot_row_b);

                // Update the basis inverse
                self.basis_inverse_rows[edit_row_index] = self.basis_inverse_rows[edit_row_index]
                    .add_multiple_of_row(
                        -column_value,
                        &self.basis_inverse_rows[pivot_row_index],
                    );
            }
        }
    }

    /// Multiplies the submatrix consisting of `minus_pi` and B^-1 by a original_column.
    ///
    /// # Arguments
    ///
    /// * `original_column`: A `SparseVector<T>` of length `m`.
    ///
    /// # Return value
    ///
    /// A `SparseVector<T>` of length `m`.
    fn generate_column(&self, original_column: SparseVector<F>) -> SparseVector<F> {
        debug_assert_eq!(original_column.len(), self.m());

        let tuples = (0..self.m())
            .map(|i| self.basis_inverse_rows[i].inner_product(&original_column))
            .enumerate()
            .filter(|&(_, v)| v != F::additive_identity())
            .collect();

        SparseVector::new(tuples, self.m())
    }

    fn generate_row(&self, i: usize, original_matrix: SparseMatrix<F, ColumnMajorOrdering>) -> SparseVector<F> {
        debug_assert!(i < self.m());
        debug_assert_eq!(original_matrix.nr_rows(), self.m());

        self.basis_inverse_rows[i].multiply(original_matrix)
    }

    fn generate_element(&self, i: usize, original_column: SparseVector<F>) -> F {
        debug_assert!(i < self.m());
        debug_assert_eq!(original_column.len(), self.m());

        self.basis_inverse_rows[i].inner_product(&original_column)
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
    fn row_reduce_update_minus_pi_and_obj(&mut self, pivot_row_index: usize, column_value: F) {
        self.minus_objective -= column_value * self.b[pivot_row_index];

        for &(column_index, value) in self.basis_inverse_rows[pivot_row_index].iter_values() {
            self.minus_pi.shift_value(column_index, - column_value * value);
        }
    }

    /// Calculates the cost difference `c_j`.
    ///
    /// This cost difference is the inner product of `minus_pi` and the column.
    fn cost_difference(&self, column: SparseVector<F>) -> F {
        debug_assert_eq!(column.len(), self.m());

        let mut total = F::additive_identity();
        for (j, value) in column.values_into_iter() {
            total += self.minus_pi[j] * value;
        }
        total
    }

    /// Get the `i`th constraint value.
    ///
    /// # Arguments
    ///
    /// * `i`: Index of the constraint to retrieve value from, in range `0` until `m`.
    ///
    /// # Return value
    ///
    /// The constraint value.
    fn get_constraint_value(&self, i: usize) -> F {
        self.b[i]
    }

    /// Get the objective function value for the current basis.
    ///
    /// # Return value
    ///
    /// The objective value.
    fn get_objective_function_value(&self) -> F {
        -self.minus_objective
    }

    /// A property of the dimensions of this matrix.
    ///
    /// A `CarryMatrix` is always square matrix of size `m + 1` times `m + 1`. Here, `m` is the
    /// length of the constraint column vector and the width of the `minus_pi` row vector. Also, the
    /// basis inverse submatrix B^-1 has dimension `m` times `m`.
    fn m(&self) -> usize {
        self.b.len() // == self.minus_pi.len()
                     // == self.basis_inverse_rows.len()
                     // == ...
                     // == ...
                     // == ...
                     // == self.basis_inverse_rows[self.m() - 1].len()
    }
}

impl<F: Field> Display for CarryMatrix<F> {
    fn fmt(&self, f: &mut Formatter) -> FormatResult {
        writeln!(f, "CarryMatrix:\n============")?;
        writeln!(f, "Objective function value: {}", -self.minus_objective)?;
        writeln!(f, "Minus PI:")?;
        <DenseVector<F> as Display>::fmt(&self.minus_pi, f)?;
        writeln!(f, "b:")?;
        <DenseVector<F> as Display>::fmt(&self.b, f)?;
        writeln!(f, "B^-1:")?;
        for row in 0..self.m() {
            for column in 0..self.m() {
                write!(f, "\t{}", self.basis_inverse_rows[row][column])?;
            }
            writeln!(f, "")?;
        }
        writeln!(f, "")
    }
}

/// Wraps a `MatrixProvider` deleting some of it's constraint rows, variable bounds should not be
/// deleted.
///
/// Used for deleting duplicate constraints after finding primal feasibility.
pub struct RemoveRows<'a, F: Field, MP: MatrixProvider<F>> {
    provider: &'a MP,
    /// List of rows that this method removes.
    ///
    /// Sorted at all times.
    pub rows_to_skip: Vec<usize>,
    // TODO: Consider implementing a cache

    phantom_data: PhantomData<F>,
}

impl<'a, F: Field, MP: MatrixProvider<F>> RemoveRows<'a, F, MP> {
    /// Create a new `RemoveRows` instance.
    ///
    /// # Arguments
    ///
    /// * `provider`: Reference to an instance implementing the `MatrixProvider` trait. Rows from
    /// this provider will be removed.
    /// * `rows_to_skip`: A **sorted** list of rows that are skipped.
    ///
    /// # Return value
    ///
    /// A new `RemoveRows` instance.
    pub fn new(provider: &'a MP, rows_to_skip: Vec<usize>) -> Self {
        debug_assert!(rows_to_skip.is_sorted());
        debug_assert_eq!(rows_to_skip.iter().collect::<HashSet<_>>().len(), rows_to_skip.len());

        RemoveRows { provider, rows_to_skip, phantom_data: PhantomData, }
    }

    /// Get the index of the same row in the original `MatrixProvider`.
    ///
    /// # Arguments
    ///
    /// * `i`: Index of row in the version of the problem from which rows were removed (this
    /// struct).
    ///
    /// # Return value
    ///
    /// Index of row in the original problem.
    pub fn get_underlying_row_index(&self, i: usize) -> usize {
        debug_assert!(i < self.provider.nr_rows() - self.rows_to_skip.len());

        Self::get_underlying_index(&self.rows_to_skip, i)
    }

    /// Delete a row
    ///
    /// # Arguments
    ///
    /// * `i`: Index of row in the version of the problem from which rows were removed (this
    /// struct)
    pub fn delete_row(&mut self, i: usize) {
        debug_assert!(i < self.provider.nr_rows() - self.rows_to_skip.len());

        let in_original_problem = self.get_underlying_row_index(i);
        debug_assert!(self.rows_to_skip.contains(&in_original_problem));
        let insertion_index = match self.rows_to_skip.binary_search(&in_original_problem) {
            Ok(_) => panic!("Deleting a row that already was deleted!"),
            Err(nr) => nr,
        };
        self.rows_to_skip.insert(insertion_index, in_original_problem);
    }

    pub fn nr_constraints_deleted(&self) -> usize {
        self.rows_to_skip.len()
    }
}

impl<'a, F: Field, MP: 'a + MatrixProvider<F>> RemoveRows<'a, F, MP> {
    /// Method abstracting over the row and column getter methods.
    ///
    /// # Arguments
    ///
    /// * `i`: Index in the reduced version of the problem.
    ///
    /// # Return value
    ///
    /// Index in the original problem.
    fn get_underlying_index(skip_indices_array: &Vec<usize>, i: usize) -> usize {
        if skip_indices_array.len() == 0 {
            // If no indices have been deleted, it's just the original value
            i
        } else if skip_indices_array.len() == 1 {
            // If one index has been deleted, see if that was before of after the value tested
            if i < skip_indices_array[0] {
                i
            } else {
                i + skip_indices_array.len()
            }
        } else {
            if i < skip_indices_array[0] {
                i
            } else {
                // Binary search with invariants:
                //   1. skip_indices_array[lower_bound] - lower_bound <= i
                //   2. skip_indices_array[upper_bound] - upper_bound > i
                let (mut lower_bound, mut upper_bound) = (0, skip_indices_array.len());
                while upper_bound - lower_bound != 1 {
                    let middle = (lower_bound + upper_bound) / 2;
                    if skip_indices_array[middle] - middle <= i {
                        lower_bound = middle
                    } else {
                        upper_bound = middle
                    }
                }

                i + upper_bound
            }
        }
    }

    /// Method abstracting over the row and column deletion methods.
    ///
    /// # Arguments
    ///
    /// * `i`: Index in the reduced version of the problem to be deleted from the original problem.
    fn delete_index(skip_indices_array: &mut Vec<usize>, i: usize) {
        let in_original_problem = Self::get_underlying_index(skip_indices_array, i);
        debug_assert!(skip_indices_array.contains(&in_original_problem));

        let insertion_index = match skip_indices_array.binary_search(&in_original_problem) {
            Ok(_) => panic!("Deleting an index that already was deleted."),
            Err(nr) => nr,
        };
        skip_indices_array.insert(insertion_index, in_original_problem);
    }
}


impl<'a, F, MP> MatrixProvider<F> for RemoveRows<'a, F, MP>
where F: Field, MP: MatrixProvider<F> {
    fn column(&self, j: usize) -> SparseVector<F> {
        let mut column = self.provider.column(j);
        column.remove_indices(&self.rows_to_skip);
        column
    }

    fn cost_value(&self, j: usize) -> F {
        self.provider.cost_value(j)
    }

    fn constraint_values(&self) -> DenseVector<F> {
        let mut all = self.provider.constraint_values();
        all.remove_indices(&self.rows_to_skip);
        all
    }

    fn bound_row_index(&self, j: usize, bound_type: BoundType) -> Option<usize> {
        self.provider.bound_row_index(j, bound_type).map(|nr| nr - self.nr_constraints_deleted())
    }

    fn bounds(&self, j: usize) -> (F, Option<F>) {
        self.provider.bounds(j)
    }

    /// This implementation assumes that only constraint rows are removed from the `MatrixProvider`.
    fn nr_constraints(&self) -> usize {
        self.provider.nr_constraints() - self.nr_constraints_deleted()
    }

    /// This implementation assumes that only constraint rows are removed from the `MatrixProvider`.
    fn nr_bounds(&self) -> usize {
        self.provider.nr_bounds()
    }

    fn nr_rows(&self) -> usize {
        self.provider.nr_rows() - self.nr_constraints_deleted()
    }

    fn nr_columns(&self) -> usize {
        self.provider.nr_columns()
    }
}

impl<'a, F, MP> VariableFeasibilityLogic<F> for RemoveRows<'a, F, MP>
where F: Field, MP: MatrixProvider<F> + VariableFeasibilityLogic<F> {
    fn is_feasible(&self, j: usize, value: F) -> bool {
        self.provider.is_feasible(j, value)
    }

    fn closest_feasible(&self, j: usize, value: F) -> (Option<F>, Option<F>) {
        self.provider.closest_feasible(j, value)
    }
}


#[cfg(test)]
mod test {
    use std::collections::HashSet;

    use num::FromPrimitive;
    use num::rational::Ratio;

    use crate::algorithm::simplex::data::{Artificial, CarryMatrix, NonArtificial, RemoveRows, Tableau};
    use crate::algorithm::simplex::matrix_provider::matrix_data::MatrixData;
    use crate::algorithm::simplex::strategy::pivot_rule::{FirstProfitable, PivotRule};
    use crate::data::linear_algebra::vector::{DenseVector, SparseVector};
    use crate::data::linear_algebra::vector::test::TestVector;
    use crate::data::number_types::traits::RealField;
    use crate::RF;
    use crate::tests::problem_2::{artificial_tableau_form, matrix_data_form};

    fn tableau<RF: RealField>(data: &MatrixData<RF>) -> Tableau<RF, NonArtificial, MatrixData<RF>> {
        let carry = {
            let minus_objective = RF!(-6);
            let minus_pi = DenseVector::from_test_data(vec![1f64, -1f64, -1f64]);
            let b = DenseVector::from_test_data(vec![1f64, 2f64, 3f64]);
            let basis_inverse_rows = vec![
                SparseVector::from_test_data(vec![1f64, 0f64, 0f64]),
                SparseVector::from_test_data(vec![-1f64, 1f64, 0f64]),
                SparseVector::from_test_data(vec![-1f64, 0f64, 1f64]),
            ];
            CarryMatrix::<RF>::create(minus_objective, minus_pi, b, basis_inverse_rows)
        };
        let basis_indices = vec![2, 3, 4];
        let mut basis_columns = HashSet::new();
        basis_columns.insert(2);
        basis_columns.insert(3);
        basis_columns.insert(4);

        Tableau::new_with_basis(
            data,
            carry,
            basis_indices,
            basis_columns,
        )
    }

    #[test]
    fn test_cost() {
        let data = matrix_data_form::<Ratio<i32>>();
        let artificial_tableau = artificial_tableau_form(&data);
        assert_eq!(artificial_tableau.objective_function_value(), Ratio::from_i32(8).unwrap());

        let tableau = tableau(&data);
        assert_eq!(tableau.objective_function_value(), Ratio::from_i32(6).unwrap());
    }

    #[test]
    fn test_relative_cost() {
        let data = matrix_data_form::<Ratio<i32>>();
        let artificial_tableau = artificial_tableau_form(&data);
        assert_eq!(artificial_tableau.relative_cost(0), Ratio::from_i32(0).unwrap());

        assert_eq!(
            artificial_tableau.relative_cost(artificial_tableau.nr_artificial_variables() + 0),
            Ratio::from_i32(-10).unwrap(),
        );

        let tableau = tableau(&data);
        assert_eq!(tableau.relative_cost(0), Ratio::from_i32(-3).unwrap());
        assert_eq!(tableau.relative_cost(1), Ratio::from_i32(-3).unwrap());
        assert_eq!(tableau.relative_cost(2), Ratio::from_i32(0).unwrap());
    }

    #[test]
    fn test_generate_column() {
        let data = matrix_data_form::<Ratio<i32>>();
        let artificial_tableau = artificial_tableau_form(&data);

        let index_to_test = artificial_tableau.nr_artificial_variables() + 0;
        let column = artificial_tableau.generate_column(index_to_test);
        let expected = SparseVector::from_test_data(vec![3f64, 5f64, 2f64]);
        assert_eq!(column, expected);
        let result = artificial_tableau.relative_cost(index_to_test);
        assert_eq!(result, Ratio::<i32>::from_i32(-10).unwrap());

        let tableau = tableau(&data);
        let index_to_test = 0;
        let column = tableau.generate_column(index_to_test);
        let expected = SparseVector::from_test_data(vec![3f64, 2f64, -1f64]);
        assert_eq!(column, expected);
        let result = tableau.relative_cost(index_to_test);
        assert_eq!(result, Ratio::<i32>::from_i32(-3).unwrap());
    }

    #[test]
    fn test_bring_into_basis() {
        let data = matrix_data_form::<Ratio<i64>>();
        let mut artificial_tableau = artificial_tableau_form(&data);
        let column = artificial_tableau.nr_artificial_variables() + 0;
        let column_data = artificial_tableau.generate_column(column);
        let row = artificial_tableau.select_primal_pivot_row(&column_data).unwrap();
        let cost = artificial_tableau.relative_cost(column);
        artificial_tableau.bring_into_basis(column, row, &column_data, cost);

        assert!(artificial_tableau.is_in_basis(&column));
        assert!(!artificial_tableau.is_in_basis(&0));
        assert_eq!(artificial_tableau.objective_function_value(), Ratio::<i64>::new(14, 3));

        let mut tableau = tableau(&data);
        let column = 1;
        let column_data = tableau.generate_column(column);
        let row = tableau.select_primal_pivot_row(&column_data).unwrap();
        let cost = tableau.relative_cost(column);
        tableau.bring_into_basis(column, row, &column_data, cost);

        assert!(tableau.is_in_basis(&column));
        assert_eq!(tableau.objective_function_value(), Ratio::<i64>::new(9, 2));
    }

    fn bfs_tableau<F: RealField>(data: &MatrixData<F>) -> Tableau<F, Artificial, MatrixData<F>> {
        let m = 3;
        let carry = {
            let minus_objective = F::from_i32(0).unwrap();
            let minus_pi = DenseVector::from_test_data(vec![1f64, 1f64, 1f64]);
            let b = DenseVector::from_test_data(vec![1f64, 2f64, 3f64]);
            let basis_inverse_rows = vec![
                SparseVector::from_test_data(vec![1f64, 0f64, 0f64]),
                SparseVector::from_test_data(vec![-1f64, 1f64, 0f64]),
                SparseVector::from_test_data(vec![-1f64, 0f64, 1f64]),
            ];
            CarryMatrix::<F>::create(minus_objective, minus_pi, b, basis_inverse_rows)
        };
        let basis_indices = vec![m + 2, m + 3, m + 4];
        let mut basis_columns = HashSet::new();
        basis_columns.insert(m + 2);
        basis_columns.insert(m + 3);
        basis_columns.insert(m + 4);

        Tableau::new_with_basis(data, carry, basis_indices, basis_columns)
    }

    #[test]
    fn test_create_tableau() {
        let data = matrix_data_form::<Ratio<i32>>();
        let bfs_tableau = bfs_tableau(&data);
        let mut rule = <FirstProfitable as PivotRule<Artificial>>::new();
        assert!(rule.select_primal_pivot_column(&bfs_tableau).is_none());
    }

    #[test]
    fn test_get_underlying_index() {
        for (deleted, size) in vec![
            (vec![2, 5, 7, 9, 12, 15, 16, 19, 20, 21], 25),
            (vec![2], 5),
            (vec![2, 3], 6),
        ] {
            let left_from_original = (0..size).filter(|i| !deleted.contains(i)).collect::<Vec<_>>();
            for (in_reduced, in_original) in (0..left_from_original.len()).zip(left_from_original.into_iter()) {
                assert_eq!(RemoveRows::<Ratio<i32>, MatrixData<Ratio<i32>>>::get_underlying_index(&deleted, in_reduced), in_original);
            }
        }
    }
}
