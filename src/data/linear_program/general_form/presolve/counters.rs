//! # Counters
//!
//! A set of counters that makes searching during presolving unnecessary.
use relp_num::{Field, OrderedField, OrderedFieldRef};

use crate::data::linear_algebra::matrix::{RowMajor, Sparse};
use crate::data::linear_algebra::SparseTuple;
use crate::data::linear_algebra::traits::SparseElement;
use crate::data::linear_program::general_form::GeneralForm;
use crate::data::linear_program::general_form::presolve::NonZeroSign;

/// Avoiding searching during presolving.
pub(super) struct Counters<'a, F: Field> {
    /// Amount of meaningful elements still in the column or row.
    /// The elements should at least be considered when the counter drops below 2. This also
    /// depends on whether the variable appears in the cost function.
    pub variable: Vec<usize>,
    /// The elements should at least be reconsidered when the counter drops below 2.
    pub constraint: Vec<usize>,
    /// Number of bounds that are missing before the activity bound can be computed.
    ///
    /// If only one bound is missing, a variable bound can be computed. If none are missing, the
    /// entire variable bound can be computed.
    pub activity: Vec<(usize, usize)>,

    /// Row major representation of the constraint matrix (a copy of `generalform.constraints`).
    rows: Sparse<&'a F, F, RowMajor>,
    general_form: &'a GeneralForm<F>,
}

impl<'a, OF> Counters<'a, OF>
where
    OF: OrderedField + SparseElement<OF>,
    for<'r> &'r OF: OrderedFieldRef<OF>,
{
    /// Create a new instance.
    ///
    /// Create a row major representation of the problem, which is expensive. Doing this once allows
    /// quick iteration over rows, which is helpful for interacting with the constraints.
    ///
    /// # Arguments
    ///
    /// * `general_form`: Problem being presolved.
    ///
    /// # Return value
    ///
    /// A new instance.
    pub fn new(general_form: &'a GeneralForm<OF>) -> Self {
        let rows: Sparse<_, _, _> = Sparse::from_column_major_ordered_matrix_although_this_is_expensive(
            &general_form.constraints,
        );

        Self {
            constraint: (0..general_form.nr_active_constraints())
                .map(|i| rows.data[i].len())
                .collect(),
            variable: (0..general_form.nr_active_variables())
                .map(|j| general_form.constraints.data[j].len())
                .collect(),
            activity: rows.iter_rows().map(|row| {
                row.iter()
                    .map(|&(j, coefficient)| {
                        let (lower, upper) =  (&general_form.variables[j].lower_bound, &general_form.variables[j].upper_bound);
                        match coefficient.signum() {
                            NonZeroSign::Positive => (lower, upper),
                            NonZeroSign::Negative => (upper, lower),
                        }
                    })
                    .fold((0, 0), |(lower_total, upper_total), (lower, upper)| {
                        let is_missing = |option: &Option<_>| match option.as_ref() {
                            Some(_) => 0,
                            None => 1,
                        };
                        (lower_total + is_missing(lower), upper_total + is_missing(upper))
                    })
            }).collect(),

            rows,
            general_form,
        }
    }

    /// Iterate over the constraints of a column who have not (yet) been eliminated.
    ///
    /// # Arguments
    ///
    /// * `column`: Column to iter over.
    /// * `constraints`: The row major representation of the constraints of the general form for
    /// fast iteration.
    ///
    /// # Return value
    ///
    /// An iterator of (row index, reference to value) tuples.
    pub(crate) fn iter_active_column(
        &self,
        variable: usize,
    ) -> impl Iterator<Item = SparseTuple<&OF>> {
        self.general_form.constraints.iter_column(variable)
            .map(|&(i, ref v)| (i, v))
            .filter(move |&(i, _)| self.is_constraint_still_active(i))
    }

    /// Iterate over the columns of a constraint that have not yet been eliminated.
    ///
    /// During presolving, for each column, a count is being kept of the number of active (belonging
    /// to a constraint that has not yet been removed) column. When that count is equal to zero, the
    /// coefficient in the original matrix of active variable coefficients is neglected.
    ///
    /// # Arguments
    ///
    /// * `constraint`: Constraint to iter over.
    ///
    /// # Return value
    ///
    /// A collection of (column index, coefficient value) tuples.
    pub(crate) fn iter_active_row(&self, constraint: usize) -> impl Iterator<Item = SparseTuple<&OF>> {
        debug_assert!(self.is_constraint_still_active(constraint));

        self.rows.iter_row(constraint)
            .copied()
            .filter(move |&(j, _)| self.is_variable_still_active(j))
    }

    /// The constraint counter indicates whether the constraint still has any variables left.
    ///
    /// Note that this can be zero, even though the constraint has not yet been eliminated, because
    /// it still needs to be checked by `presolve_empty_constraint`. It can, in that case, however
    /// be ignored during the application of presolving rules.
    pub(crate) fn is_constraint_still_active(&self, constraint: usize) -> bool {
        self.constraint[constraint] > 0
    }

    /// The variable counter indicates whether the variable still appears in any constraint.
    pub(crate) fn is_variable_still_active(&self, variable: usize) -> bool {
        self.variable[variable] > 0
    }
}
