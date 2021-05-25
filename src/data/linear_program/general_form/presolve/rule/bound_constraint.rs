//! # Remove a constraint that is a bound on a variable.
//!
//! Triggered when there is only a single variable in a constraint.
use relp_num::{OrderedField, OrderedFieldRef};
use relp_num::NonZeroSign;

use crate::data::linear_algebra::traits::SparseElement;
use crate::data::linear_program::elements::{BoundDirection, LinearProgramType, RangedConstraintRelation};
use crate::data::linear_program::general_form::presolve::Index;
use crate::data::linear_program::general_form::presolve::updates::BoundChange;

impl<'a, OF> Index<'a, OF>
where
    OF: OrderedField + SparseElement<OF>,
    for<'r> &'r OF: OrderedFieldRef<OF>,
{
    /// Remove a constraint that is a bound on a variable.
    ///
    /// This bound is either a bound on a slack variable, or it has no slack variables in it.
    ///
    /// # Arguments
    ///
    /// * `constraint`: Index of a row with only a bound.
    ///
    /// # Return value
    ///
    /// `Result::Err` if the linear program was determined infeasible.
    pub (in crate::data::linear_program::general_form::presolve) fn presolve_bound_constraint(
        &mut self,
        constraint: usize,
    ) -> Result<(), LinearProgramType<OF>> {
        use RangedConstraintRelation::{Equal, Range, Less, Greater};
        use NonZeroSign::{Positive, Negative};
        debug_assert_eq!(self.counters.constraint[constraint], 1);
        debug_assert_eq!(self.counters.iter_active_row(constraint).count(), 1);

        let (variable, coefficient) = self.counters.iter_active_row(constraint).next().unwrap();
        debug_assert_ne!(self.counters.variable[variable], 0);

        let bound_value = self.updates.b(constraint) / coefficient;
        let mut changes = Vec::with_capacity(2);
        match (self.updates.constraint_type(constraint), coefficient.signum()) {
            (Greater, Positive) | (Less, Negative) => {
                changes.push((BoundDirection::Lower, bound_value));
            },
            (Less, Positive) | (Greater, Negative) => {
                changes.push((BoundDirection::Upper, bound_value));
            },
            (Equal, _) => {
                changes.push((BoundDirection::Lower, bound_value.clone()));
                changes.push((BoundDirection::Upper, bound_value));
            },
            (Range(range), _) => {
                let bound1 = (self.updates.b(constraint) - range) / coefficient;
                let bound2 = bound_value;

                match coefficient.signum() {
                    Positive => {
                        changes.push((BoundDirection::Lower, bound1));
                        changes.push((BoundDirection::Upper, bound2));
                    }
                    Negative => {
                        changes.push((BoundDirection::Lower, bound2));
                        changes.push((BoundDirection::Upper, bound1));
                    }
                }
            },
        }

        self.counters.variable[variable] -= 1;
        self.counters.constraint[constraint] -= 1;
        self.remove_constraint(constraint);

        for (direction, bound_value) in changes {
            match self.updates.update_bound(variable, direction, bound_value) {
                BoundChange::None => {}
                BoundChange::NewBound => self.after_bound_change(variable, direction, None),
                BoundChange::BoundShift(shift) => self.after_bound_change(variable, direction, Some(shift)),
            }
        }

        match self.updates.variable_feasible_value(variable) {
            None => Err(LinearProgramType::Infeasible),
            Some(_) => self.queue_variable_by_counter(variable),
        }
    }
}
