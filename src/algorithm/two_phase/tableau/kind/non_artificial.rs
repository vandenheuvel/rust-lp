//! # Non-Artificial Tableau
//!
//! Contains a tableau `Kind` type that acts as a simple passthrough for a matrix provider and
//! `Tableau` logic that is only relevant in the second phase.
use std::collections::HashSet;

use crate::algorithm::two_phase::matrix_provider::column::{Column, ColumnNumber};
use crate::algorithm::two_phase::matrix_provider::filter::Filtered;
use crate::algorithm::two_phase::matrix_provider::MatrixProvider;
use crate::algorithm::two_phase::tableau::inverse_maintenance::{InverseMaintener, ops as im_ops};
use crate::algorithm::two_phase::tableau::kind::Kind;
use crate::algorithm::two_phase::tableau::Tableau;

/// The `TableauType` in case the `Tableau` does not contain any artificial variables.
///
/// This `Tableau` variant should only be constructed with a known feasible basis.
#[derive(Eq, PartialEq, Debug)]
pub struct NonArtificial<'a, MP> {
    /// Supplies data about the problem.
    ///
    /// This data doesn't change throughout the lifetime of this `Tableau`, and it is independent of
    /// the current basis as described by the `carry` and `basis_columns` attributes.
    provider: &'a MP,
}
impl<'provider, MP> Kind for NonArtificial<'provider, MP>
where
    MP: MatrixProvider,
{
    type Column = MP::Column;
    type Cost = MP::Cost<'provider>;

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
    fn initial_cost_value(&self, j: usize) -> Self::Cost {
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

impl<'provider, IM, MP> Tableau<IM, NonArtificial<'provider, MP>>
where
    IM: InverseMaintener<F:
        im_ops::FieldHR +
        im_ops::Column<<MP::Column as Column>::F> +
        im_ops::Cost<MP::Cost<'provider>> +
    >,
    MP: MatrixProvider,
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
    pub fn new_with_inverse_maintainer(
        provider: &'provider MP,
        inverse_maintainer: IM,
        basis_columns: HashSet<usize>,
    ) -> Self {
        Tableau {
            inverse_maintainer,
            basis_columns,

            kind: NonArtificial {
                provider,
            },
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
    ) -> Self
    where
        IM::F: im_ops::Column<MP::Rhs>,
        MP::Rhs: 'static,
    {
        let arbitrary_order = basis.iter().copied().collect::<Vec<_>>();

        Tableau {
            inverse_maintainer: IM::from_basis(&arbitrary_order, provider),
            basis_columns: basis.clone(),

            kind: NonArtificial {
                provider,
            },
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
        basis_indices: HashSet<usize>,
        provider: &'provider MP,
    ) -> Self {
        Tableau {
            inverse_maintainer: IM::from_artificial(
                inverse_maintainer,
                provider,
                nr_artificial,
            ),
            // Shift the basis column indices back
            basis_columns: basis_indices.into_iter()
                .map(|column| column - nr_artificial)
                .collect(),

            kind: NonArtificial {
                provider,
            },
        }
    }
}

impl<'provider, IM, MP> Tableau<IM, NonArtificial<'provider, MP>>
where
    IM: InverseMaintener<F: im_ops::FieldHR + im_ops::Column<<MP::Column as Column>::F> + im_ops::Cost<MP::Cost<'provider>>>,
    MP: Filtered,
{
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
    pub fn from_artificial_removing_rows(
        inverse_maintainer: IM,
        nr_artificial: usize,
        mut basis: HashSet<usize>,
        provider: &'provider MP,
    ) -> Self {
        debug_assert!(basis.iter().all(|&v| v >= nr_artificial || provider.filtered_rows().contains(&v)));

        for &row in provider.filtered_rows() {
            let was_there = basis.remove(&inverse_maintainer.basis_column_index_for_row(row));
            debug_assert!(was_there);
        }
        let basis_columns = basis.into_iter().map(|j| j - nr_artificial).collect();

        // Remove same row and column from carry matrix
        let inverse_maintainer = IM::from_artificial_remove_rows(
            inverse_maintainer,
            provider,
            nr_artificial,
        );

        Tableau {
            inverse_maintainer,
            basis_columns,

            kind: NonArtificial {
                provider,
            },
        }
    }
}
