use std::collections::HashSet;
use std::marker::PhantomData;

use crate::algorithm::two_phase::matrix_provider::MatrixProvider;
use crate::algorithm::two_phase::matrix_provider::filter::generic_wrapper::{IntoFilteredColumn, RemoveRows};
use crate::algorithm::two_phase::tableau::inverse_maintenance::InverseMaintenance;
use crate::algorithm::two_phase::tableau::kind::Kind;
use crate::algorithm::two_phase::tableau::Tableau;
use crate::algorithm::utilities::remove_indices;
use crate::data::linear_algebra::traits::SparseElementZero;
use crate::data::number_types::traits::{Field, FieldRef};

/// The `TableauType` in case the `Tableau` does not contain any artificial variables.
///
/// This `Tableau` variant should only be constructed with a known feasible basis.
#[derive(Eq, PartialEq, Debug)]
pub struct NonArtificial<'a, F: 'static, FZ, MP> {
    /// Supplies data about the problem.
    ///
    /// This data doesn't change throughout the lifetime of this `Tableau`, and it is independent of
    /// the current basis as described by the `carry` and `basis_columns` attributes.
    provider: &'a MP,

    phantom: PhantomData<(F, FZ)>,
}
impl<'a, F: 'static, FZ, MP> Kind<F, FZ> for NonArtificial<'a, F, FZ, MP>
where
    F: Field,
    FZ: SparseElementZero<F>,
    MP: MatrixProvider<F, FZ>,
{
    type Column = MP::Column;

    /// Coefficient of variable `j` in the objective function.
    ///
    /// # Arguments
    ///
    /// * `tableau`: Tableau to retrieve the cost value from.
    /// * `j`: Column index of the variable, in range 0 until `self.nr_columns()`.
    ///
    /// # Return value
    ///
    /// The cost of variable `j`.
    fn initial_cost_value(&self, j: usize) -> &F {
        debug_assert!(j < self.provider.nr_columns());

        self.provider.cost_value(j)
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
        debug_assert!(j < self.provider.nr_columns());

        self.provider.column(j)
    }

    fn nr_rows(&self) -> usize {
        self.provider.nr_rows()
    }

    fn nr_columns(&self) -> usize {
        self.provider.nr_columns()
    }
}

impl<'provider, F: 'static, FZ, IM, MP> Tableau<F, FZ, IM, NonArtificial<'provider, F, FZ, MP>>
where
    F: Field,
    for<'r> &'r F: FieldRef<F>,
    FZ: SparseElementZero<F>,
    IM: InverseMaintenance<F, FZ>,
    MP: MatrixProvider<F, FZ>,
{
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
    pub(crate) fn new_with_inverse_maintainer(
        provider: &'provider MP,
        inverse_maintainer: IM,
        basis_indices: Vec<usize>,
        basis_columns: HashSet<usize>,
    ) -> Self {
        Tableau {
            inverse_maintainer,
            basis_indices,
            basis_columns,

            kind: NonArtificial {
                provider,

                phantom: PhantomData,
            },

            phantom: PhantomData,
        }
    }

    /// Creates a Simplex tableau with a specific basis.
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
        // TODO(OPTIMIZATION): Order doesn't matter, document and make this a Vec
        basis: &HashSet<usize>,
    ) -> Self {
        let arbitrary_order = basis.iter().copied().collect::<Vec<_>>();

        Tableau {
            inverse_maintainer: IM::from_basis(&arbitrary_order, provider),
            basis_indices: arbitrary_order,
            basis_columns: basis.clone(),

            kind: NonArtificial {
                provider,

                phantom: PhantomData,
            },

            phantom: PhantomData,
        }
    }

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
        inverse_maintainer: IM,
        nr_artificial: usize,
        basis: (Vec<usize>, HashSet<usize>),
        provider: &'provider MP,
    ) -> Self {
        let (mut basis_indices, basis_columns) = basis;

        // Shift the basis column indices back
        basis_indices.iter_mut().for_each(|column| *column -= nr_artificial);

        Tableau {
            inverse_maintainer: IM::from_artificial(
                inverse_maintainer,
                provider,
                // TODO(CORRECTNESS): Should these be shifted?
                &basis_indices,
            ),
            basis_indices,
            basis_columns: basis_columns.into_iter()
                .map(|column| column - nr_artificial)
                .collect(),

            kind: NonArtificial {
                provider,

                phantom: PhantomData,
            },

            phantom: PhantomData,
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
    pub fn from_artificial_removing_rows<'b: 'provider>(
        inverse_maintainer: IM,
        nr_artificial: usize,
        basis: (Vec<usize>, HashSet<usize>),
        rows_removed: &'b RemoveRows<'provider, F, FZ, MP>,
    ) -> Tableau<F, FZ, IM, NonArtificial<'b, F, FZ, RemoveRows<'provider, F, FZ, MP>>>
    where
        MP::Column: IntoFilteredColumn<F>,
    {
        debug_assert!(basis.0.iter().all(|&v| v >= nr_artificial || rows_removed.rows_to_skip.contains(&v)));

        let (mut basis_indices, mut basis_columns) = basis;
        for &row in &rows_removed.rows_to_skip {
            let was_there = basis_columns.remove(&basis_indices[row]);
            debug_assert!(was_there);
        }
        let basis_columns = basis_columns.into_iter().map(|j| j - nr_artificial).collect();

        remove_indices(&mut basis_indices, &rows_removed.rows_to_skip);
        basis_indices.iter_mut().for_each(|index| *index -= nr_artificial);

        // Remove same row and column from carry matrix
        let inverse_maintainer = IM::from_artificial_remove_rows(
            inverse_maintainer,
            rows_removed,
            &basis_indices,
        );

        Tableau {
            inverse_maintainer,
            basis_indices,
            basis_columns,

            kind: NonArtificial {
                provider: rows_removed,

                phantom: PhantomData,
            },

            phantom: PhantomData,
        }
    }
}