//! # Remove a variable that is a slack.
//!
//! Triggered when there is only a single constraint in a variable and it does not appear in the
//! objective function.
use crate::data::linear_algebra::traits::SparseElement;
use crate::data::linear_program::elements::{BoundDirection, LinearProgramType};
use crate::data::linear_program::elements::RangedConstraintRelation;
use crate::data::linear_program::general_form::presolve::Index;
use crate::data::linear_program::general_form::RemovedVariable;
use crate::data::number_types::traits::{OrderedField, OrderedFieldRef};
use crate::data::number_types::nonzero::NonzeroSign;
use crate::data::number_types::nonzero::sign::Signed;

impl<'a, OF> Index<'a, OF>
where
    OF: OrderedField + SparseElement<OF>,
    for<'r> &'r OF: OrderedFieldRef<OF>,
{
    /// Try to remove a variable that appears in exactly one constraint.
    ///
    /// This method attempts to remove slack variable that can be viewed as a slack variable. This
    /// is a variable that does not appear in the cost function and only in a single constraint.
    ///
    /// TODO(CORRECTNESS): A slack constraint might also become irrelevant, remove it in that case.
    ///
    /// # Cause and effect
    ///
    /// One of three actions will be taken by this function:
    /// * Adapt the right-hand side b and remove the column
    /// * Remove the column and the row
    /// * Do nothing
    ///
    /// A combination of three factors determines what needs to happen:
    /// * Constraint direction (either <=, ==, >= or =r=)
    /// * Bound on the "slack" (either <=, >=, <=>= or none)
    /// * Sign of the "slack"'s constant (+ or -)
    ///
    /// # Specification
    ///
    /// We start with a constraint `<a, x> + cs (<= == >= =r=) b`
    ///
    /// When c > 0:
    ///     ||      s >= l      |      s <= u      |         l <= s <= u        | free |
    /// --------------------------------------------------------------------------------
    /// <=  || <a, x> <= b - cl |                  |      <a, x> <= b - cl      |      |
    /// ==  || <a, x> <= b - cl | <a, x> >= b - cu |   <a, x> =c(u-l)= b - cl   |      |
    /// >=  ||                  | <a, x> >= b - cu |      <a, x> >= b - cu      |      |
    /// =r= || <a, x> <= b - cl | <a, x> >= b - cu | <a, x> =r + c(u-l)= b - cl |      |
    ///
    /// When c < 0:
    ///     ||      s >= l      |      s <= u      |         l <= s <= u        | free |
    /// --------------------------------------------------------------------------------
    /// <=  ||                  | <a, x> <= b - cu |      <a, x> <= b - cu      |      |
    /// ==  || <a, x> >= b - cl | <a, x> <= b - cu |   <a, x> =c(l-u)= b - cu   |      |
    /// >=  || <a, x> >= b - cl |                  |      <a, x> >= b - cl      |      |
    /// =r= || <a, x> >= b - cl | <a, x> <= b - cu | <a, x> =r + c(l-u)= b - cu |      |
    ///
    /// Note that nothing changes when the constraint is of type == and s has both a lower and
    /// upper bound (regardless of the sign of c). If a bound is adapted, the column is removed
    /// afterward. If no constraint remains, both the column and the row are removed.
    ///
    /// # Arguments
    ///
    /// * `variable_index`: Index of variable that should be removed if it is a slack with suitable
    /// bounds.
    pub (in crate::data::linear_program::general_form::presolve) fn presolve_slack(
        &mut self,
        variable: usize,
    ) -> Result<(), LinearProgramType<OF>> {
        use RangedConstraintRelation::{Less, Equal, Greater, Range};
        use NonzeroSign::{Positive, Negative};
        debug_assert_eq!(self.counters.variable[variable], 1);
        debug_assert_eq!(self.counters.iter_active_column(variable).count(), 1);
        debug_assert_eq!(self.general_form.variables[variable].cost, OF::zero());
        debug_assert!(self.updates.is_variable_fixed(variable).is_none());

        let (constraint, coefficient) = self.counters.iter_active_column(variable)
            // Only coefficient in the problem
            .next().unwrap();
        // TODO(ARCHITECTURE): This clone is necessary because the borrow can't be split enough.
        let coefficient = coefficient.clone();
        let constraint_type = self.updates.constraint_type(constraint);

        let lower = self.updates.variable_bound(variable, BoundDirection::Lower);
        let upper = self.updates.variable_bound(variable, BoundDirection::Upper);
        // TODO(ARCHITECTURE): This clone is necessary because the borrow can't be split enough.
        let bounds = (lower.cloned(), upper.cloned());
        // We see whether these are `None` here already such that `bounds` can be moved in the match
        let bounds_is_none = (lower.is_none(), upper.is_none());
        let coefficient_sign = coefficient.signum();

        // Patterns are ordered left to right, top to bottom and then by coefficient (positive to
        // negative)
        let (new_constraint_type, bound) = match (constraint_type, bounds, coefficient_sign) {
            // Remove both row and column
            (Greater, (Some(_), None), Positive) | (Less, (None, Some(_)), Positive) |
            (Less, (Some(_), None), Negative) | (Greater, (None, Some(_)), Negative) |
            (_, (None, None), _) => {
                // Remove also the row, 2 * (1 + 1 + 3) = 10 cases

                // The order of the calls below is a bit critical because of the strict assumptions
                // they make about the state of counters, updates etc.

                let solution = self.compute_removed_variable_solution(constraint, variable, &coefficient);

                let variables_to_scan = self.counters.iter_active_row(constraint)
                    .map(|(j, _)| j).collect::<Vec<_>>();
                for other_variable in variables_to_scan {
                    self.counters.constraint[constraint] -= 1;
                    self.counters.variable[other_variable] -= 1;
                    if other_variable != variable {
                        self.queue_variable_by_counter(other_variable)?;
                    }
                }

                self.remove_variable(variable, solution);
                self.remove_constraint(constraint);

                return Ok(());
            },
            // All cases below remove the column only
            (Equal, (Some(lower), Some(upper)), Positive) => {
                // <a, x> =c(u-l)= b - cl
                // `upper` and `lower` are not equal, because then the variable would have been
                // substituted already as a fixed variable.
                (RangedConstraintRelation::Range(&coefficient * (upper - &lower)), lower)
            },
            (Equal, (Some(lower), Some(upper)), Negative) => {
                // <a, x> =c(l-u)= b - cu
                // `upper` and `lower` are not equal, because then the variable would have been
                // substituted already as a fixed variable.
                (RangedConstraintRelation::Range(&coefficient * (lower - &upper)), upper)
            }
            (Less | Equal | Range(_), (Some(bound), None), Positive) | (Less, (Some(bound), Some(_)), Positive) => {
                // <a, x> <= b - cl
                (RangedConstraintRelation::Less, bound)
            },
            (Equal | Greater | Range(_), (None, Some(bound)), Positive) | (Greater, (Some(_), Some(bound)), Positive) => {
                // <a, x> >= b - cu
                (RangedConstraintRelation::Greater, bound)
            },
            (Equal | Greater | Range(_), (Some(bound), None), Negative) | (Greater, (Some(bound), Some(_)), Negative) => {
                // <a, x> >= b - cl
                (RangedConstraintRelation::Greater, bound)
            },
            (Less | Equal | Range(_), (None, Some(bound)), Negative) | (Less, (Some(_), Some(bound)), Negative) => {
                // <a, x> <= b - cu
                (RangedConstraintRelation::Less, bound)
            },
            (Range(range), (Some(lower), Some(upper)), Positive) => {
                // <a, x> =r + c(u-l)= b - cl
                (RangedConstraintRelation::Range(range + &coefficient * (upper - &lower)), lower)
            },
            (Range(range), (Some(lower), Some(upper)), Negative) => {
                // <a, x> =r + c(l-u)= b - cu
                (RangedConstraintRelation::Range(range + &coefficient * (lower - &upper)), upper)
            },
        };

        let change = -&coefficient * &bound;

        let removed = match constraint_type {
            Equal | Range(_) => {
                self.compute_removed_variable_solution(constraint, variable, &coefficient)
            }
            Greater | Less => RemovedVariable::Solved(bound),
        };
        self.counters.variable[variable] -= 1;
        self.remove_variable(variable, removed);

        self.update_activity_queues_if_needed(constraint, bounds_is_none, coefficient_sign);
        self.counters.constraint[constraint] -= 1;
        self.queue_constraint_by_counter(constraint).map(|_| ())?;

        self.updates.change_b(constraint, change);
        self.updates.constraints.insert(constraint, new_constraint_type);

        Ok(())
    }

    fn update_activity_queues_if_needed(
        &mut self,
        constraint: usize,
        bounds: (bool, bool),
        coefficient_sign: NonzeroSign,
    ) {
        use NonzeroSign::{Positive, Negative};
        if matches!((bounds, coefficient_sign), ((true, _), Positive) | ((_, true), Negative)) {
            self.counters.activity[constraint].0 -= 1;
            if self.counters.activity[constraint].0 <= 1 {
                self.queues.activity.insert(constraint, BoundDirection::Lower);
            }
        }
        if matches!((bounds, coefficient_sign), ((_, true), Positive) | ((true, _), Negative)) {
            self.counters.activity[constraint].1 -= 1;
            if self.counters.activity[constraint].1 <= 1 {
                self.queues.activity.insert(constraint, BoundDirection::Upper);
            }
        }
    }

    fn compute_removed_variable_solution(
        &self,
        constraint: usize,
        variable: usize,
        coefficient: &OF,
    ) -> RemovedVariable<OF> {
        let constant = self.updates.b(constraint) / coefficient;
        let coefficients = self.counters.iter_active_row(constraint)
            .filter(|&(j, _)| j != variable)
            .map(|(j, other_coefficient)| {
                (self.general_form.from_active_to_original[j], other_coefficient / coefficient)
            })
            .collect();
        RemovedVariable::FunctionOfOthers { constant, coefficients, }
    }
}
