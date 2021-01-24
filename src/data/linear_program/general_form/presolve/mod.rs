//! # Presolving linear programs
//!
//! A `GeneralForm` can be presolved by building an index, repeatedly applying reduction rules and
//! applying the changes proposed. This module contains data structures and logic for presolving.
use std::iter::Iterator;

use crate::data::linear_algebra::traits::SparseElement;
use crate::data::linear_program::elements::{BoundDirection, LinearProgramType, RangedConstraintRelation};
use crate::data::linear_program::general_form::{GeneralForm, RemovedVariable};
use crate::data::linear_program::general_form::presolve::counters::Counters;
use crate::data::linear_program::general_form::presolve::queues::Queues;
use crate::data::linear_program::general_form::presolve::updates::Updates;
use crate::data::number_types::traits::{Field, OrderedField, OrderedFieldRef};

mod rule;
mod queues;
pub(super) mod updates;
mod counters;
pub(super) mod scale;
pub use scale::Scalable as Prescalable;
use crate::data::number_types::nonzero::NonzeroSign;
use crate::data::number_types::nonzero::sign::Signed;

/// Container data structure to keep track of presolve status.
///
/// Queues are used to determine which rules still have to be applied to which constraint or
/// variable indices. Proposed changes that don't change the solution to the final problem are
/// collected in the Updates field. A few indices, counters, are kept to speed up the process: they
/// indicate when which rules can be applied.
pub(super) struct Index<'a, F: Field> {
    /// Which rule should still be applied to which constraint or variable.
    queues: Queues,
    /// Set of changes that doesn't change the solution to the problem.
    pub updates: Updates<'a, F>,
    /// Indices to speed up solving.
    counters: Counters<'a, F>,

    /// We maintain the computed activity bounds.
    activity_bounds: Vec<(Option<F>, Option<F>)>,

    /// Column major representation of the constraint matrix, borrowed from problem being presolved.
    general_form: &'a GeneralForm<F>,
}

/// Whether a change that was made was "meaningful".
///
/// It can happen that the presolve algorithm gets caught in a never ending loop where it repeatedly
/// tightens bounds, but never reaches the solution, slowly converging. This return value indicates
/// to the caller whether anything "meaningful" happened, such as the creation of a new bound or the
/// removal of a constraint, or not, such as when a bound is tightened.
#[derive(Eq, PartialEq, Debug)]
pub(super) enum Change {
    Meaningful,
    NotMeaningful,
    None,
}

impl<'a, OF> Index<'a, OF>
where
    OF: OrderedField + SparseElement<OF>,
    for<'r> &'r OF: OrderedFieldRef<OF>,
{
    /// Create a new instance.
    ///
    /// This operation is expensive, it creates a row major copy of the constraint matrix. The
    /// creation of some of the indices might also be expensive.
    /// TODO(OPTIMIZATION): Can this copy be created only when needed?
    /// TODO(OPTIMIZATION): Are these counters and queues created in a single pass? (compiler)
    ///
    /// Note that all indices of columns and variables in the attributes of this struct are
    /// relative to the active variables in the problem in its state before the presolving started.
    /// These indices get updated (e.g. due to the removal of constraints from the general form's
    /// data structures) only after this struct is dropped.
    ///
    /// # Arguments
    ///
    /// * `general_form`: Problem description which is being presolved.
    pub(super) fn new(general_form: &'a GeneralForm<OF>) -> Result<Self, LinearProgramType<OF>> {
        let counters = Counters::new(general_form);
        let updates = Updates::new(general_form, &counters)?;
        let queues = Queues::new(general_form, &counters);

        Ok(Self {
            queues,
            updates,
            counters,

            activity_bounds: vec![(None, None); general_form.nr_constraints()],

            general_form: &general_form,
        })
    }

    /// Apply a single presolve rule.
    ///
    /// The order of the rules is an estimate of the ratio between how likely a rule yields a useful
    /// result, and how expensive it is to apply the rule.
    ///
    /// The following rules are applied:
    ///
    /// 1. Substitute fixed variable
    /// 2. Remove a bound constraint
    /// 3. Eliminate slack variable
    /// 4. Domain propagation (constraint activation bounds)
    ///
    /// Generally, a constraint can only directly influence a variable, not another constraint. The
    /// same holds true for variables. As such, we can view the constraints and variables as a
    /// biparte graph.
    ///
    /// TODO(ENHANCEMENT): What is the best order to apply the rules in?
    /// TODO(ENHANCEMENT): Which element should be removed from the queue first?
    ///
    /// # Arguments
    ///
    /// * `&mut self`: A `PresolveIndex` that lives across repeated calls to this function. It is
    /// used to store which constraints and variables need to be checked for which rules. After the
    /// application of each rule, constraints and variables might be added to queues in this struct.
    ///
    /// # Return value
    ///
    /// `Ok` indicating whether a "meaningful change" was made. This could for example be a row or
    /// column being removed, but also a constraint type being adapted. This is used to avoid
    /// never-ending "improvements" of bounds that don't make the problem easier (but numerically
    /// harder, unless the problem is rescaled). If the program is determined to be infeasible, an
    /// `Err` type.
    pub(super) fn presolve_step(&mut self) -> Result<Change, LinearProgramType<OF>> {
        // Rules that are guaranteed to make the problem smaller
        if let Some(variable) = self.queues.substitution.pop() {
            return self.presolve_fixed_variable(variable)
                // always removes a variable
                .map(|()| Change::Meaningful);
        }
        while let Some(constraint) = self.queues.bound.pop() {
            if self.counters.is_constraint_still_active(constraint) {
                return self.presolve_bound_constraint(constraint)
                    // always removes a constraint
                    .map(|()| Change::Meaningful);
            }
        }
        while let Some(variable) = self.queues.slack.pop() {
            if self.counters.is_variable_still_active(variable) {
                return self.presolve_slack(variable)
                    // always removes a variable
                    .map(|()| Change::Meaningful);
            }
        }

        // Rules not guaranteed to make the problem smaller
        while let Some((constraint, direction)) = self.queues.activity.pop() {
            if self.counters.is_constraint_still_active(constraint) {
                // Rules with higher priority are applied when a row or column counter equals 1.
                debug_assert!(self.counters.constraint[constraint] > 1, "was {}", self.counters.constraint[constraint]);
                return self.presolve_domain_propagation(constraint, direction);
            }
        }

        // TODO(ENHANCEMENT): Duplicate rows and columns.

        // No rule was applied
        Ok(Change::NotMeaningful)
    }

    /// Performs actions that should be performed after a new variable bound is found.
    ///
    /// These include:
    ///
    /// * If a variable is now fixed, it should be substituted in the problem.
    /// * Otherwise, some rows might need to have the existing activity bound recomputed
    ///
    /// # Arguments
    ///
    /// * `variable`: Variable who's bound was changed.
    /// * `direction`: Whether the lower- or upper variable bound was updated.
    /// * `change`: Whether another bound was previously known, and the difference between the
    /// former and current bound.
    fn after_bound_change(
        &mut self,
        variable: usize,
        direction: BoundDirection,
        change: Option<OF>,
    ) {
        debug_assert!(self.updates.removed_variables.iter().all(|&(j, _)| variable != j));
        debug_assert!(match direction {
            BoundDirection::Lower => change.as_ref().map_or(true, |v| v > &OF::zero()),
            BoundDirection::Upper => change.as_ref().map_or(true, |v| v < &OF::zero()),
        });

        if self.updates.is_variable_fixed(variable).is_some() && self.counters.is_variable_still_active(variable) {
            self.queues.substitution.push(variable);
        }

        match change {
            Some(difference) => self.update_activity_bounds(variable, direction, difference),
            None => self.update_activity_counters(variable, direction),
        }
    }

    /// Update an activity bound if it exists.
    ///
    /// Activity bounds don't have to be recomputed entirely after a single variable bound change.
    /// You need to know by how much the bound was changed, which gets lost. So it is recomputed
    /// now, even though the result might not be directly used, and other recomputations might be
    /// triggered first (in the future, perhaps due to the counter for the number of recomputations
    /// reaching it's limit (precision)).
    ///
    /// # Arguments
    ///
    /// * `variable`: Variable who's bound was changed.
    /// * `direction`: Whether a variable upper or lower bound was changed.
    /// * `by_how_much`: Size of the change.
    fn update_activity_bounds(
        &mut self,
        variable: usize,
        direction: BoundDirection,
        by_how_much: OF,
    ) {
        debug_assert!(match direction {
            BoundDirection::Lower => by_how_much > OF::zero(),
            BoundDirection::Upper => by_how_much < OF::zero(),
        });

        let rows_to_check = self.counters.iter_active_column(variable).collect::<Vec<_>>();
        for (row, coefficient) in rows_to_check {
            if !self.counters.is_constraint_still_active(row) {
                continue;
            }

            let bound_to_edit = direction * coefficient.signum();
            if let Some(ref mut bound) = match bound_to_edit {
                BoundDirection::Lower => &mut self.activity_bounds[row].0,
                BoundDirection::Upper => &mut self.activity_bounds[row].1,
            } {
                *bound += &by_how_much * coefficient;
                debug_assert!(match bound_to_edit {
                    BoundDirection::Lower => self.counters.activity[row].0,
                    BoundDirection::Upper => self.counters.activity[row].1,
                } <= 1);
                self.queues.activity.insert(row, bound_to_edit);
            }
        }
    }

    /// Update the activity counters after a new bound was found.
    ///
    /// # Arguments
    ///
    /// * `variable`: Variable who's bound was added.
    /// * `direction`: Whether a variable upper or lower bound was added.
    fn update_activity_counters(
        &mut self,
        variable: usize,
        direction: BoundDirection,
    ) {
        let constraints_to_check = self.counters.iter_active_column(variable)
            // TODO(ARCHITECTURE): Avoid This clone
            .map(|(i, v)| (i, v.clone())).collect::<Vec<_>>();
        for (constraint, coefficient) in constraints_to_check {
            let activity_direction = direction * coefficient.signum();
            let counter = match activity_direction {
                BoundDirection::Lower => &mut self.counters.activity[constraint].0,
                BoundDirection::Upper => &mut self.counters.activity[constraint].1,
            };
            *counter -= 1;
            if *counter <= 1 {
                self.queues.activity.insert(constraint, activity_direction);
            }
        }
    }

    /// Mark a constraint as removed.
    ///
    /// There should be more than one element in this row, otherwise, the column that can be removed
    /// should be known from the rule application and removed directly, such that an iteration can
    /// be avoided.
    ///
    /// # Arguments
    ///
    /// * `constraint`: Constraint to iter over.
    fn remove_constraint_values(
        &mut self,
        constraint: usize,
    ) -> Result<(), LinearProgramType<OF>> {
        // TODO(ARCHITECTURE): Avoid this clone and collect.
        let variables_to_scan = self.counters.iter_active_row(constraint)
            .map(|(j, _)| j)
            .collect::<Vec<_>>();
        for variable in variables_to_scan {
            self.counters.constraint[constraint] -= 1;
            self.counters.variable[variable] -= 1;
            self.queue_variable_by_counter(variable)?;
        }

        debug_assert_eq!(self.counters.constraint[constraint], 0);
        Ok(())
    }

    /// When the variable counter drops low, this has implications for rules that should be tested.
    ///
    /// Note that if a variable has no coefficients, we solve it directly in this function.
    fn queue_variable_by_counter(&mut self, variable: usize) -> Result<(), LinearProgramType<OF>> {
        match self.counters.variable[variable] {
            0 => {
                // If a variable is unfeasible before presolving, should have been detected during
                // loading as a trivial infeasibility. If unfeasible later, should be detected at
                // bound change.
                debug_assert!(self.updates.variable_feasible_value(variable).is_some());

                let value = if self.general_form.variables[variable].cost.is_zero() {
                    RemovedVariable::Solved(self.updates.variable_feasible_value(variable).unwrap())
                } else {
                    self.updates.optimize_column_independently(variable)?
                };
                self.remove_variable(variable, value);
            },
            1 => if self.general_form.variables[variable].cost.is_zero() {
                self.queues.slack.push(variable);
            },
            _ => (),
        }

        Ok(())
    }

    /// If the number of elements in a constraint drops below a certain number, something might need
    /// to happen.
    ///
    /// We try to remove elements from the constraint matrix the entire time; the number of elements
    /// in each row is counted such that we don't need to search for the rows where we might apply
    /// some of the most powerful rules, the ones that work with certainty and make the problem
    /// smaller.
    ///
    /// Empty rows are handled directly in this function.
    ///
    /// # Arguments
    ///
    /// * `constraint`: The constraint whos counter might have becomes sufficiently low to activate
    /// a rule.
    ///
    /// # Returns
    ///
    /// Whether a meaningful change was made in the function (such as when the constraint is empty).
    fn queue_constraint_by_counter(&mut self, constraint: usize) -> Result<Change, LinearProgramType<OF>> {
        match self.counters.constraint[constraint] {
            0 => {
                let right_hand_side = self.updates.b(constraint);
                let constraint_type = self.updates.constraint_type(constraint);
                if is_empty_constraint_feasible(right_hand_side, constraint_type) {
                    self.remove_constraint(constraint);
                    Ok(Change::Meaningful)
                } else {
                    Err(LinearProgramType::Infeasible)
                }
            }
            1 => {
                self.queues.bound.push(constraint);
                Ok(Change::None)
            },
            _ => Ok(Change::None),
        }
    }

    /// Mark a constraint as removed.
    ///
    /// Contains a debug check to see whether all variables have been removed. While this makes the
    /// order of statements that execute a removing operation tricky, it's a simple check that
    /// catches bugs.
    fn remove_constraint(&mut self, constraint: usize) {
        debug_assert_eq!(self.counters.constraint[constraint], 0);

        self.updates.constraints_marked_removed.push(constraint);
    }

    /// Mark a variable as removed.
    ///
    /// Contains a debug check to see whether all variables have been removed. While this makes the
    /// order of statements that execute a removing operation tricky, it's a simple check that
    /// catches bugs.
    fn remove_variable(&mut self, variable: usize, solution: RemovedVariable<OF>) {
        debug_assert_eq!(self.counters.variable[variable], 0);

        self.updates.removed_variables.push((variable, solution));
    }

    /// Whether there are no more rules in the queue to be applied.
    pub fn are_queues_empty(&self) -> bool {
        self.queues.are_empty()
    }
}

/// Whether an empty constraint indicates infeasibility.
///
/// This method will be called when a constraints has no coefficients left. In the case, the
/// constraint should still be satisfied, or the problem is infeasible.
///
/// # Arguments
///
/// * `right_hand_side`: Bound value.
/// * `constraint_type`: Relation that is supposed to hold between right hand side and the inner
/// product of the relevant matrix row with the variable (<a, x>), which in this case is zero.
fn is_empty_constraint_feasible<OF>(
    right_hand_side: &OF,
    constraint_type: &RangedConstraintRelation<OF>,
) -> bool
where
    OF: OrderedField,
    for<'r> &'r OF: OrderedFieldRef<OF>,
{
    match constraint_type {
        RangedConstraintRelation::Equal => {
            right_hand_side == &OF::zero()
        },
        RangedConstraintRelation::Range(range) => {
            right_hand_side >= &OF::zero() && right_hand_side - range <= OF::zero()
        },
        RangedConstraintRelation::Less => {
            right_hand_side >= &OF::zero()
        },
        RangedConstraintRelation::Greater => {
            right_hand_side <= &OF::zero()
        },
    }
}

#[cfg(test)]
mod test;
