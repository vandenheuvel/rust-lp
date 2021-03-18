//! # Basis initially completely from artificial variables
//!
//! When no initial pivots can be derived, this module is used. It is slightly quicker than the
//! partially artificial tableau kind.
use crate::algorithm::two_phase::matrix_provider::column::Column;
use crate::algorithm::two_phase::matrix_provider::column::identity::IdentityColumn;
use crate::algorithm::two_phase::matrix_provider::MatrixProvider;
use crate::algorithm::two_phase::tableau::inverse_maintenance::{InverseMaintener, ops as im_ops};
use crate::algorithm::two_phase::tableau::kind::artificial::{Artificial, Cost};
use crate::algorithm::two_phase::tableau::kind::Kind;
use crate::algorithm::two_phase::tableau::Tableau;

/// All variables are artificial.
pub struct Fully<'a, MP: MatrixProvider> {
    /// Supplies data about the problem.
    ///
    /// This data doesn't change throughout the lifetime of this `Tableau`, and it is independent of
    /// the current basis as described by the `carry` and `basis_columns` attributes.
    provider: &'a MP,
}

impl<'provider, MP> Kind for Fully<'provider, MP>
where
    MP: MatrixProvider<Column: Column + IdentityColumn>,
{
    type Column = MP::Column<'provider>;
    type Cost = Cost;

    fn initial_cost_value(&self, j: usize) -> Self::Cost {
        if j < self.nr_artificial_variables() {
            Cost::One
        } else {
            Cost::Zero
        }
    }

    fn original_column(&self, j: usize) -> Self::Column {
        if j < self.nr_rows() {
            // TODO(ENHANCEMENT): Would it be possible to specialize the code where this identity
            //  column is used?
            <Self::Column as IdentityColumn>::identity(j, self.nr_rows())
        } else {
            self.provider.column(j - self.nr_rows())
        }
    }

    fn nr_rows(&self) -> usize {
        self.provider.nr_rows()
    }

    fn nr_columns(&self) -> usize {
        self.nr_rows() + self.provider.nr_columns()
    }
}

impl<'provider, MP> Artificial for Fully<'provider, MP>
where
    MP: MatrixProvider<Column: IdentityColumn>,
{
    fn nr_artificial_variables(&self) -> usize {
        self.nr_rows()
    }

    fn pivot_row_from_artificial(&self, artificial_index: usize) -> usize {
        debug_assert!(artificial_index < self.nr_artificial_variables());

        artificial_index
    }
}

impl<'provider, IM, MP> Tableau<IM, Fully<'provider, MP>>
where
    IM: InverseMaintener<F: im_ops::Rhs<MP::Rhs>>,
    MP: MatrixProvider,
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
    pub(crate) fn new(provider: &'provider MP) -> Self {
        let m = provider.nr_rows();

        Tableau {
            inverse_maintainer: IM::create_for_fully_artificial(provider.right_hand_side()),
            basis_columns: (0..m).collect(),

            kind: Fully { provider },
        }
    }
}
