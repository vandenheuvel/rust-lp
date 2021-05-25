//! # Subsitute a fixed variable
//!
//! Adapt the right hand side for the remaining constraint and potentially, update the problem's
//! fixed cost.
use relp_num::{OrderedField, OrderedFieldRef};

use crate::data::linear_algebra::traits::SparseElement;
use crate::data::linear_program::elements::LinearProgramType;
use crate::data::linear_program::general_form::presolve::Index;
use crate::data::linear_program::general_form::RemovedVariable;

impl<'a, OF> Index<'a, OF>
where
    OF: OrderedField + SparseElement<OF>,
    for<'r> &'r OF: OrderedFieldRef<OF>,
{
    /// Substitute a variable with a known value in the constraints in which it appears.
    ///
    /// # Arguments
    ///
    /// * `variable`: Index of column under consideration.
    pub(in crate::data::linear_program::general_form::presolve) fn presolve_fixed_variable(
        &mut self,
        variable: usize,
    ) -> Result<(), LinearProgramType<OF>> {
        debug_assert!(!self.updates.removed_variables.iter().any(|&(j, _)| j == variable));
        debug_assert!(self.updates.is_variable_fixed(variable).is_some());

        let value = self.updates.is_variable_fixed(variable).unwrap().clone();

        for (constraint, coefficient) in self.counters.iter_active_column(variable) {
            self.updates.change_b(constraint, -coefficient * &value);
        }
        self.updates.fixed_cost += &self.general_form.variables[variable].cost * &value;

        // Can't combine these loops because of the reference in `coefficient` clashing with the
        // counter update.
        let rows_to_substitute = self.counters.iter_active_column(variable)
            .map(|(i, _)| i)
            .collect::<Vec<_>>();
        for constraint in rows_to_substitute {
            self.counters.variable[variable] -= 1;
            self.counters.constraint[constraint] -= 1;
            self.queue_constraint_by_counter(constraint)?;
        }

        self.remove_variable(variable, RemovedVariable::Solved(value));
        Ok(())
    }
}
