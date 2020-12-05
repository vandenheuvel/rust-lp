//! # Artificial variables in the tableau
//!
//! Representing artificial variables in the tableau "virtually". That is, the values aren't
//! actually stored but instead logic provides the algorithm with the input it needs to drive the
//! artificial variables out of the basis.
use std::collections::HashSet;

use crate::algorithm::two_phase::matrix_provider::Column;
use crate::algorithm::two_phase::tableau::inverse_maintenance::{ExternalOps, InverseMaintenance};
use crate::algorithm::two_phase::tableau::kind::Kind;
use crate::algorithm::two_phase::tableau::Tableau;

pub mod fully;
pub mod partially;

/// tableaus with artificial variables.
///
/// There are currently two implementations; either all variables are artificial, or not necessarily
/// all variables are. See the two submodules.
pub trait Artificial: Kind<Column: IdentityColumn> {
    /// How many artificial variables are in the tableau.
    ///
    /// This number varies, because slack variables might have been recognized as practical
    /// candidates for basic feasible solutions by the `MatrixProvider` (the
    /// `positive_slack_indices` method).
    ///
    /// # Return value
    ///
    /// This number can be zero (for non artificial tableaus, represented by the `NonArtificial`
    /// struct), or any number through the number of rows (`self.nr_rows`).
    fn nr_artificial_variables(&self) -> usize;

    /// At which row is the pivot from a specific artificial variable located?
    ///
    /// # Arguments
    ///
    /// * `artificial_index`: Index of artificial variable.
    ///
    /// # Returns
    ///
    /// Row index where the pivot is located.
    fn pivot_row_from_artificial(&self, artificial_index: usize) -> usize;
}

/// Identity columns are needed for artificial matrices.
///
/// When a matrix provider is to be used in the first phase, it should be possible to represent
/// identity columns.
pub trait IdentityColumn: Column {
    /// Create an identity column, placing a "1" at a certain index and "0"'s otherwise.
    ///
    /// # Arguments
    ///
    /// * `i`: Index at which the "1" should be placed.
    /// * `len`: Length of the column. Might not be used in an actual implementation.
    fn identity(i: usize, len: usize) -> Self;
}

/// Functionality needed only, and for all, artificial tableaus.
///
/// Most of these functions get called in the artificial simplex method, or the method that removes
/// artificial variables from the problem at zero level.
impl<'provider, IM, A> Tableau<IM, A>
where
    IM: InverseMaintenance<F: ExternalOps<<A::Column as Column>::F>>,
    A: Artificial,
{
    /// Whether there are any artificial variables in the basis.
    pub fn has_artificial_in_basis(&self) -> bool {
        self.basis_columns.iter().any(|&c| c < self.nr_artificial_variables())
    }

    /// Get the indices of the artificial variables that are still in the basis.
    pub fn artificial_basis_columns(&self) -> HashSet<usize> {
        self.basis_columns
            .iter()
            .filter(|&&v| v < self.nr_artificial_variables())
            .copied()
            .collect()
    }

    /// At which row is the pivot from a specific artificial variable located?
    ///
    /// # Arguments
    ///
    /// * `artificial_index`: Index of artificial variable.
    ///
    /// # Returns
    ///
    /// Row index where the pivot is located.
    pub fn pivot_row_from_artificial(&self, artificial_index: usize) -> usize {
        debug_assert!(artificial_index < self.nr_artificial_variables());
        // Only used to remove variables from basis
        debug_assert!(self.is_in_basis(&artificial_index));

        self.kind.pivot_row_from_artificial(artificial_index)
    }

    /// Number of artificial variables in this tableau.
    pub fn nr_artificial_variables(&self) -> usize {
        self.kind.nr_artificial_variables()
    }

    /// Extract information necessary to construct a `NonArtificial` tableau.
    ///
    /// # Returns
    ///
    /// The inverse maintainer (typically contains the inverse of a basis, expensive to redo), the
    /// number of artificial variables that were in the basis (all assumed to have been located at
    /// the lowest indices), and a tuple of which the first element maps rows of the problem (before
    /// any rows were removed, if necessary) to the columns holding the pivots from the current
    /// basis, as well as a set copy of those columns.
    pub fn into_basis(self) -> (IM, usize, (Vec<usize>, HashSet<usize>)) {
        let nr_artificial = self.nr_artificial_variables();
        (self.inverse_maintainer, nr_artificial, (self.basis_indices, self.basis_columns))
    }
}
