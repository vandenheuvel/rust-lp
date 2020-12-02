//! # Basis initially with a few slacks
//!
//! If a model can provide a few initial slack variables, those are used such that a basic feasible
//! solution is found quicker. Less variables need to be driven out of the basis.
use std::collections::HashSet;
use std::marker::PhantomData;

use crate::algorithm::two_phase::matrix_provider::MatrixProvider;
use crate::algorithm::two_phase::PartialInitialBasis;
use crate::algorithm::two_phase::tableau::inverse_maintenance::InverseMaintenance;
use crate::algorithm::two_phase::tableau::kind::artificial::{Artificial, IdentityColumn};
use crate::algorithm::two_phase::tableau::kind::Kind;
use crate::algorithm::two_phase::tableau::Tableau;
use crate::data::linear_algebra::traits::SparseElementZero;
use crate::data::number_types::traits::{Field, FieldRef};

/// The `TableauType` in case the `Tableau` contains artificial variables.
#[derive(Eq, PartialEq, Debug)]
pub struct Partially<'a, F, FZ, MP> {
    /// For the `i`th artificial variable was originally a simple basis vector with the coefficient
    /// at index `column_to_row[i]`.
    column_to_row: Vec<usize>,

    /// Values that can be referred to when unsized constants need to be returned.
    ///
    /// TODO(ARCHITECTURE): Replace with values that are Copy, or an enum?
    /// TODO(ARCHITECTURE): Rename the (single) method that uses these to shift the the relevant
    ///  value to be able to remove these fields.
    ONE: F,
    ZERO: F,
    /// Supplies data about the problem.
    ///
    /// This data doesn't change throughout the lifetime of this `Tableau`, and it is independent of
    /// the current basis as described by the `carry` and `basis_columns` attributes.
    provider: &'a MP,
    phantom_zero: PhantomData<FZ>,
}
impl<'a, F: 'static, FZ, MP> Partially<'a, F, FZ, MP>
where
    F: Field,
    FZ: SparseElementZero<F>,
    MP: MatrixProvider<F, FZ>,
{
    fn nr_artificial_variables(&self) -> usize {
        self.column_to_row.len()
    }
}
impl<'a, F: 'static, FZ, MP> Kind<F, FZ> for Partially<'a, F, FZ, MP>
where
    F: Field,
    FZ: SparseElementZero<F>,
    MP: MatrixProvider<F, FZ, Column: IdentityColumn<F>>,
{
    type Column = MP::Column;

    /// Coefficient of variable `j` in the objective function.
    ///
    /// # Arguments
    ///
    /// * `j`: Column index of the variable, in range 0 until self.nr_columns().
    ///
    /// # Return value
    ///
    /// The cost of variable `j`.
    fn initial_cost_value(&self, j: usize) -> &F {
        debug_assert!(j < self.nr_artificial_variables() + self.provider.nr_columns());

        if j < self.nr_artificial_variables() {
            &self.ONE
        } else {
            &self.ZERO
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
    fn original_column(&self, j: usize) -> Self::Column {
        debug_assert!(j < self.nr_columns());

        if j < self.nr_artificial_variables() {
            <Self::Column as IdentityColumn<F>>::identity(self.column_to_row[j], self.nr_rows())
        } else {
            self.provider.column(j - self.nr_artificial_variables())
        }
    }

    fn nr_rows(&self) -> usize {
        self.provider.nr_rows()
    }

    fn nr_columns(&self) -> usize {
        self.nr_artificial_variables() + self.provider.nr_columns()
    }
}
impl<'provider, F: 'static, FZ, MP> Artificial<F, FZ> for Partially<'provider, F, FZ, MP>
where
    F: Field + 'provider,
    for<'r> &'r F: FieldRef<F>,
    FZ: SparseElementZero<F>,
    MP: MatrixProvider<F, FZ, Column: IdentityColumn<F>>,
{
    fn nr_artificial_variables(&self) -> usize {
        self.column_to_row.len()
    }

    fn pivot_row_from_artificial(&self, artificial_index: usize) -> usize {
        debug_assert!(artificial_index < self.nr_artificial_variables());

        self.column_to_row[artificial_index]
    }
}

impl<'provider, F: 'static, FZ, IM, MP> Tableau<F, FZ, IM, Partially<'provider, F, FZ, MP>>
where
    F: Field + 'provider,
    for<'r> &'r F: FieldRef<F>,
    FZ: SparseElementZero<F>,
    IM: InverseMaintenance<F, FZ>,
    MP: MatrixProvider<F, FZ>,
{
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
    pub(crate) fn new(provider: &'provider MP) -> Self
    where
        MP: PartialInitialBasis
    {
        let m = provider.nr_rows();

        // (row index, column index) coordinate tuples of suitable pivots in a slack column.
        let real = provider.pivot_element_indices();
        debug_assert!(real.is_sorted_by_key(|&(row, _column)| row));
        // Amount of slack variables that can be used for the initial artificial basis.
        let nr_real = real.len();
        // Amount of artificial variables that are still needed. Could be zero.
        let nr_artificial = m - nr_real;

        // Collect the rows that are not in the `real` value. Simple complement of the `real` set.
        let mut artificial = Vec::with_capacity(nr_artificial);
        let mut i = 0;
        for ith_artificial in 0..nr_artificial {
            while i < nr_real && ith_artificial + i == real[i].0 {
                i += 1;
            }
            artificial.push(ith_artificial + i);
            debug_assert!(!real.iter().any(|&(row, _)| row == ith_artificial + i));
        }
        debug_assert!(artificial.iter().all(|&row| row < m));

        // Create the map from the row to column indices where the initial basis elements are.
        //
        // Neither the artificial variables nor the real initial basis variables are necessarily in
        // a uninterrupted range, which is why this block is needed.
        debug_assert!(artificial.is_sorted() && real.is_sorted_by_key(|&(row, _column)| row));
        let mut artificial_counter = 0;
        let basis_indices = (0..m).map(|i| {
            let artificial_column = artificial_counter;
            let real_column = || nr_artificial + real[i - artificial_counter].1;

            let can_take_from_artificial = artificial_counter < nr_artificial;
            let can_take_from_real = i - artificial_counter < nr_real;

            match (can_take_from_artificial, can_take_from_real) {
                (true, true) => {
                    // Which index row is lower?
                    if artificial[artificial_column] < real[i - artificial_counter].0 {
                        artificial_counter += 1;
                        artificial_column
                    } else {
                        real_column()
                    }
                },
                (true, false) => {
                    artificial_counter += 1;
                    artificial_column
                },
                (false, true) => {
                    real_column()
                },
                (false, false) => {
                    unreachable!("If both `ac == nr_a` and `i - ac == nr_r`, loop should have ended")
                },
            }
        }).collect::<Vec<_>>();
        let basis_columns = basis_indices.iter().copied().collect();

        let inverse_maintainer = IM::create_for_partially_artificial(
            &artificial,
            &real,
            provider.constraint_values(),
        );

        Tableau {
            inverse_maintainer,
            basis_indices,
            basis_columns,

            kind: Partially {
                column_to_row: artificial,

                ONE: F::one(),
                ZERO: F::zero(),

                provider,
                phantom_zero: PhantomData,
            },

            phantom: PhantomData,
        }
    }

    /// Creates a Simplex tableau with a specific basis.
    ///
    /// Currently only used for testing.
    ///
    /// # Arguments
    ///
    /// * `provider`: Provides the original problem for which the other arguments describe a basis.
    /// * `carry`: `Carry` with the basis transformation. Corresponds to `basis_indices`.
    /// * `basis_indices`: Maps each row to a column, describing a basis. Corresponds to `carry`.
    ///
    /// # Return value
    ///
    /// `Tableau` with for the provided problem with the provided basis.
    pub(crate) fn new_with_basis(
        provider: &'provider MP,
        inverse_maintainer: IM,
        basis_indices: Vec<usize>,
        basis_columns: HashSet<usize>,
        column_to_row_artificials: Vec<usize>,
    ) -> Self {
        debug_assert!(column_to_row_artificials.iter().all(|i| basis_indices.contains(&i)));

        Tableau {
            inverse_maintainer,
            basis_indices,
            basis_columns,

            kind: Partially {
                column_to_row: column_to_row_artificials,

                ONE: F::one(),
                ZERO: F::zero(),

                provider,

                phantom_zero: PhantomData
            },

            phantom: PhantomData,
        }
    }
}
