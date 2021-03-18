//! # Variable logic
//!
//! If this project is ever extended to a branch and bound framework, we can generalize variables
//! as the trait in this module specifies.
use crate::algorithm::two_phase::matrix_provider::column::Column;
use crate::algorithm::two_phase::matrix_provider::MatrixProvider;

/// Logic for testing whether variables are feasible.
///
/// Defined as a separate trait from `MatrixProvider`. Matrices are defined over fields, and so
/// the `MatrixProvider` is. Some of the logic of variable feasibility is more part of linear
/// programming algorithms specifically, which are only defined over ordered fields. This logic is
/// thus separated into a different trait, which depends on the other trait.
pub trait FeasibilityLogic: MatrixProvider {
    /// Whether a variable is integer.
    ///
    /// # Arguments
    ///
    /// * `j`: Variable index.
    ///
    /// # Return value
    ///
    /// `true` if the variable is integer, `false` otherwise.
    fn is_feasible(&self, j: usize, value: <Self::Column<'_> as Column>::F) -> bool;

    /// Closest feasible variable to the left and right.
    ///
    /// # Arguments
    ///
    /// * `j`: Variable index.
    ///
    /// # Return value
    ///
    /// Two `Option`s, one for the closest feasible value to the left, one for the closest feasible
    /// value to the right. Note that these values might be equal, if there is only one feasible
    /// value.
    fn closest_feasible(
        &self,
        j: usize,
        value: <Self::Column<'_> as Column>::F,
    ) -> (Option<<Self::Column<'_> as Column>::F>, Option<<Self::Column<'_> as Column>::F>);
}
