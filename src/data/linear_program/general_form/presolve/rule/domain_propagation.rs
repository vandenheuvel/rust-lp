//! # Domain propagation
//!
//! Use variable bounds to see whether constraints can be dropped or simplified. Use bounds on
//! constraints to derive variable bounds.
use std::cmp::Ordering;

use crate::data::linear_algebra::traits::SparseElement;
use crate::data::linear_program::elements::{BoundDirection, InequalityRelation, LinearProgramType};
use crate::data::linear_program::elements::{RangedConstraintRelation};
use crate::data::linear_program::general_form::presolve::{Change, Index};
use crate::data::linear_program::general_form::presolve::updates::BoundChange;
use crate::data::number_types::traits::{OrderedField, OrderedFieldRef};

impl<'a, OF> Index<'a, OF>
where
    OF: OrderedField + SparseElement<OF>,
    for<'r> &'r OF: OrderedFieldRef<OF>,
{
    /// Attempt to tighten bounds using activity bounds.
    ///
    /// As described in Achterberg (2007), algorithm 7.1.
    ///
    /// TODO(ENHANCEMENT): Variable bounds might be derived without a constraint eventually being
    ///  removed. Should these variable bounds be tracked, and perhaps not emitted in the end?
    ///
    /// # Arguments
    ///
    /// * `constraint`: Index of the constraint under consideration.
    /// * `direction`: Whether the activity lower- or upperbound should be used.
    ///
    /// # Return value
    ///
    /// `Result::Err` if the problem was determined to be infeasible.
    pub (in crate::data::linear_program::general_form::presolve) fn presolve_domain_propagation(
        &mut self,
        constraint: usize,
        direction: BoundDirection,
    ) -> Result<Change, LinearProgramType<OF>> {
        let counter = match direction {
            BoundDirection::Lower => self.counters.activity[constraint].0,
            BoundDirection::Upper => self.counters.activity[constraint].1,
        };
        let nr_bounds_missing_in_row = self.counters.iter_active_row(constraint)
            .filter(|&(j, c)| {
                let bound_direction = direction * c.signum();
                self.updates.variable_bound(j, bound_direction).is_none()
            })
            .count();
        debug_assert_eq!(nr_bounds_missing_in_row, counter);

        match counter {
            0 => self.for_entire_constraint(constraint, direction),
            1 => Ok(self.create_variable_bound(constraint, direction)),
            _ => unreachable!("Constraint should not be added to the activity queue before the counter is 0 or 1."),
        }
    }

    /// Apply domain propagation through activation bounds for a constraint where all relevant
    /// variable bounds are known.
    ///
    /// Attempts to both prove infeasibility and removal of the entire constraint. Afterwards, it
    /// attempts to add any missing "opposite" bounds.
    ///
    /// # Arguments
    ///
    /// * `constraint`: Index of the constraint under consideration.
    /// * `direction`: Whether the activity lower- or upperbound should be used.
    fn for_entire_constraint(
        &mut self,
        constraint: usize,
        direction: BoundDirection, // Activity bound that can be calculated
    ) -> Result<Change, LinearProgramType<OF>> {
        debug_assert_eq!(
            match direction {
                BoundDirection::Lower => self.counters.activity[constraint].0,
                BoundDirection::Upper => self.counters.activity[constraint].1,
            },
            0,
        );
        // In this variable we collect the most change made by the function, upgrading when needed.
        let mut most_meaningful_change = Change::None;

        // TODO(ARCHITECTURE): Try to avoid this clone
        let activity_bound = self.compute_activity_bound_if_needed(constraint, direction).clone();

        let (remove_constraint, apply_variable_part) = self.constraint_part(
            constraint,
            &activity_bound,
            direction,
            &mut most_meaningful_change,
        )?;

        if apply_variable_part {
            if let Some(right_hand_side) = self.can_variable_rule_be_applied(constraint, direction) {
                self.variable_part(
                    constraint,
                    right_hand_side,
                    &activity_bound,
                    direction,
                    &mut most_meaningful_change,
                );
            }
        }

        if remove_constraint {
            self.remove_constraint_values(constraint)?;
            self.remove_constraint(constraint);
        }

        Ok(most_meaningful_change)
    }

    /// If no value is known for one of the activity bounds, compute them.
    ///
    /// This is a helper method for `presolve_constraint_by_domain_propagation`.
    ///
    /// # Arguments
    ///
    /// * `constraint`: Index of the constrain being checked.
    ///
    /// # Return value
    ///
    /// A tuple with a (lower, upper) activity bound. Can be `None` if not all of the relevant
    /// lower or upper bounds are known.
    fn compute_activity_bound_if_needed(
        &mut self,
        constraint: usize,
        direction: BoundDirection,
    ) -> &OF {
        debug_assert_eq!(match direction {
            BoundDirection::Lower => self.counters.activity[constraint].0,
            BoundDirection::Upper => self.counters.activity[constraint].1,
        }, 0);

        // Relevant activity counter.
        let bound = match direction {
            BoundDirection::Lower => &mut self.activity_bounds[constraint].0,
            BoundDirection::Upper => &mut self.activity_bounds[constraint].1,
        };

        // If there is no bound, compute it.
        if bound.is_none() {
            let updates = &self.updates;
            *bound = Some(
                self.counters.iter_active_row(constraint)
                    .map(|(variable, coefficient)| {
                        let bound_direction = direction * coefficient.signum();
                        let bound = updates.variable_bound(variable, bound_direction).unwrap();
                        coefficient * bound
                    })
                    .sum()
            );
        }

        // We unwrap because we just computed it.
        bound.as_ref().unwrap()
    }

    fn constraint_part(
        &mut self,
        constraint: usize,
        bound: &OF,
        direction: BoundDirection,
        made_change: &mut Change,
    ) -> Result<(bool, bool), LinearProgramType<OF>> {
        let constraint_update = self.constraint_update(constraint, &bound, direction)?;
        if let Some(change) = constraint_update {
            let (remove_constraint, apply_variable_part) = match change {
                ConstraintUpdate::Remove => (true, true),
                ConstraintUpdate::Replace(new_inequality, right_hand_side_shift) => {
                    self.updates.constraints.insert(constraint, new_inequality.into());
                    self.updates.change_b(constraint, right_hand_side_shift);

                    (false, true)
                },
                ConstraintUpdate::SetVariablesToBound => {
                    let mut activity_counters_to_update = Vec::new();
                    for (variable, coefficient) in self.counters.iter_active_row(constraint) {
                        let variable_direction = direction * coefficient.signum();
                        let bound = self.updates.variable_bound(variable, variable_direction)
                            // This value exists because the activity bound we're working with is based on it.
                            .unwrap()
                            // TODO(ARCHITECTURE): Try to avoid this clone.
                            .clone();

                        // Move the bound over to the "permanent" bound changes
                        if let Some(bound_value) = self.updates.activity_bounds.remove(&(variable, variable_direction)) {
                            self.updates.bounds.insert((variable, variable_direction), bound_value);
                        }
                        // Update the other bound to be the same
                        let change = self.updates.update_bound(variable, !variable_direction, bound.clone());
                        if change == BoundChange::NewBound {
                            activity_counters_to_update.push((variable, !variable_direction));
                        }
                        debug_assert!(self.updates.is_variable_fixed(variable).is_some());

                        self.queues.substitution.push(variable);
                    }

                    for (variable, direction) in activity_counters_to_update {
                        self.update_activity_counters(variable, direction);
                    }

                    (true, false)
                }
            };
            *made_change = Change::Meaningful;
            // Constraint gets removed outside this function
            Ok((remove_constraint, apply_variable_part))
        } else {
            Ok((false, true))
        }
    }

    /// Whether a constraint can be removed.
    ///
    /// We start with `<a, x> (== =r= <= >=) b`, `l <= <a, x>` or `<a, x> <= u` and `a (= < >) b`
    /// where `a` is one of `l` and `u`. This results in `4x2x3 = 24` combinations, all of which
    /// are cases below.
    ///
    /// Note that always `l < u`, because otherwise there there would be fixed variables (and they
    /// are presolved in a rule with higher priority). This is not directly relevant though, as we
    /// treat only one side at the time.
    ///
    /// # Arguments
    ///
    /// * `constraint`: Index of the constraint under consideration.
    /// * `value`: Value of the activity lower- or upperbound.
    /// * `direction`: Whether it's the activity lower- or upperbound.
    ///
    /// # Return value
    ///
    /// Whether a change needs to be made, and if so, which one.
    ///
    /// # Errors
    ///
    /// If the constraint can't be satisfied anymore.
    fn constraint_update(
        &self,
        constraint: usize,
        bound_value: &OF,
        direction: BoundDirection,
    ) -> Result<Option<ConstraintUpdate<OF>>, LinearProgramType<OF>> {
        use RangedConstraintRelation::{Equal, Range, Less, Greater};

        let rhs = self.updates.b(constraint);
        let constraint_type = self.updates.constraint_type(constraint);

        match (direction, constraint_type, rhs.cmp(bound_value)) {
            // `<a, x> (== =r= <=) b < l`
            (BoundDirection::Lower, Equal | Range(_) | Less, Ordering::Less) |
            // `<a, x> (== =r= >=) b > u`
            (BoundDirection::Upper, Equal | Greater, Ordering::Greater) => {
                Err(LinearProgramType::Infeasible)
            },

            // `l <= <a, x> (== <=) b == l`
            (BoundDirection::Lower, Equal | Less, Ordering::Equal) |
            // `u >= <a, x> (== >=) b = u`
            (BoundDirection::Upper, Equal | Greater, Ordering::Equal) => {
                Ok(Some(ConstraintUpdate::SetVariablesToBound))
            },

            // `<a, x> <= u (= <) b`
            (BoundDirection::Upper, Less, Ordering::Equal | Ordering::Greater) |
            // `<a, x> >= l (= >) b`
            (BoundDirection::Lower, Greater, Ordering::Less | Ordering::Equal) => {
                Ok(Some(ConstraintUpdate::Remove))
            },

            // `b - r <= <a, x> <= b = u`
            (BoundDirection::Upper, Range(range), Ordering::Equal) => {
                Ok(Some(ConstraintUpdate::Replace(InequalityRelation::Greater, -range)))
            },

            (BoundDirection::Lower, Range(range), Ordering::Greater) => {
                let lower_bound = rhs - range;
                Ok(match bound_value.cmp(&lower_bound) {
                    // `l < b - r <= <a, x> <= b`
                    Ordering::Less => None,
                    // `b - r <= l <= <a, x> <= b`
                    Ordering::Equal | Ordering::Greater => {
                        Some(ConstraintUpdate::Replace(InequalityRelation::Less, OF::zero()))
                    },
                })
            },
            (BoundDirection::Upper, Range(range), Ordering::Greater) => {
                let lower_bound = rhs - range;
                match bound_value.cmp(&lower_bound) {
                    Ordering::Less => Err(LinearProgramType::Infeasible),
                    // `u >= <a, x> <= b == u`
                    Ordering::Equal => Ok(Some(ConstraintUpdate::SetVariablesToBound)),
                    // `b - r <= <a, x> <= u < b`
                    Ordering::Greater => {
                        // Keep only lower bound
                        Ok(Some(ConstraintUpdate::Replace(InequalityRelation::Greater, -range)))
                    }
                }
            },

            // `<a, x> (== <= >=) b > l`
            (BoundDirection::Lower, Equal | Less | Greater, Ordering::Greater) |
            // `<a, x> (== =r= <= >=) b < u`
            (BoundDirection::Upper, Equal | Range(_) | Less | Greater, Ordering::Less) => {
                Ok(None)
            },

            // `l <= <a, x> =r= b == l`
            (BoundDirection::Lower, Range(_), Ordering::Equal) => {
                unreachable!("Would require the range to be zero, in which case it should be an \
                              equality bound instead.")
            },
        }
    }

    fn can_apply_variable_part(&self, constraint: usize, direction: BoundDirection) -> bool {
        matches!(
            (direction, self.updates.constraint_type(constraint)),
            (_, RangedConstraintRelation::Equal)
                | (BoundDirection::Lower, RangedConstraintRelation::Less)
                | (BoundDirection::Upper, RangedConstraintRelation::Greater),
        )
    }

    fn variable_part(
        &mut self,
        constraint: usize,
        right_hand_side: OF,
        activity_bound: &OF,
        activity_direction: BoundDirection,
        made_change: &mut Change,
    ) {
        let targets = self.counters.iter_active_row(constraint)
            .map(|(i, v)| (i, v.clone())) // TODO(ARCHITECTURE): Avoid this clone
            .collect::<Vec<_>>();
        for (variable, coefficient) in targets {
            let coefficient_sign = coefficient.signum();
            let new_variable_bound_direction = !activity_direction * coefficient_sign;

            let variable_bound_value = self.updates.variable_bound(variable, activity_direction * coefficient_sign)
                .unwrap(); // Was used to compute the activity bound
            let residual_activity_bound = activity_bound - &coefficient * variable_bound_value;

            let new_variable_bound_value = (&right_hand_side - residual_activity_bound) / coefficient;

            let change = self.updates.update_activity_variable_bound(
                variable,
                new_variable_bound_direction,
                new_variable_bound_value,
            );

            match change {
                BoundChange::None => {},
                BoundChange::NewBound => {
                    self.after_bound_change(variable, new_variable_bound_direction, None);
                    *made_change = Change::Meaningful;
                },
                BoundChange::BoundShift(change) => {
                    self.after_bound_change(variable, new_variable_bound_direction, Some(change));
                    if *made_change != Change::Meaningful {
                        *made_change = Change::NotMeaningful;
                    }
                },
            }
        }
    }

    /// Apply domain propagation through activation bounds for a constraint where all of the
    /// relevant variable bounds are known.
    ///
    /// # Arguments
    ///
    /// * `constraint`: Index of the constraint under consideration.
    /// * `variable`: Index of variable who's relevant variable bound is not known, with all other
    /// variable bounds in this constraint being known.
    /// * `coefficient`: Nonzero coefficient of the variable in the constraint.
    /// * `direction`: Whether the activity lower- or upperbound should be used.
    /// * `activity_bound`: Value of the activity lower- or upperbound.
    ///
    /// # Return value
    ///
    /// Tuple of:
    ///
    /// * Direction of the bound that was changed (lower or upper).
    /// * Change that was done (any at all, and if so, by how much or was there no bound before? See
    /// the `update_bound` method for more).
    fn try_tighten_variable_bound(
        &mut self,
        constraint: usize,
        variable: usize,
        coefficient: &OF,
        activity_direction: BoundDirection,
        activity_bound: &OF,
    ) -> (BoundDirection, BoundChange<OF>) {
        let bound_direction = activity_direction * coefficient.signum();
        let bound_value = self.updates.variable_bound(variable, bound_direction).unwrap();
        let bound = activity_bound - coefficient * bound_value;

        let variable_bound = (self.updates.b(constraint) - bound) / coefficient;
        let change = self.updates.update_activity_variable_bound(variable, !bound_direction, variable_bound);
        (!bound_direction, change)
    }

    /// Apply domain propagation through activation bounds for a constraint where all but one of the
    /// relevant variable bounds are known.
    ///
    /// # Arguments
    ///
    /// * `constraint`: Index of the constraint under consideration.
    /// * `direction`: Whether the activity lower- or upperbound should be used.
    fn create_variable_bound(
        &mut self,
        constraint: usize,
        activity_direction: BoundDirection,
    ) -> Change {
        // We should be exactly one variable bound away from being able to compute an activity
        // bound, otherwise, a different rule should be applied.
        debug_assert_eq!(
            match activity_direction {
                BoundDirection::Lower => self.counters.activity[constraint].0,
                BoundDirection::Upper => self.counters.activity[constraint].1,
            },
            1,
        );

        let right_hand_side = match self.can_variable_rule_be_applied(constraint, activity_direction) {
            None => return Change::None,
            Some(value) => value,
        };

        // Compute the activity bound and notice to which variable it can be applied (because it
        // is the only one that doesn't have the relevant bound yet).
        let total_activity = self.counters.iter_active_row(constraint)
            .filter_map(|(variable, coefficient)| {
                let bound_direction = activity_direction * coefficient.signum();
                self.updates.variable_bound(variable, bound_direction).map(|bound| coefficient * bound)
            })
            .sum::<OF>();
        let (target_column, target_coefficient) = self.counters.iter_active_row(constraint)
            .find(|&(variable, coefficient)| {
                let bound_direction = activity_direction * coefficient.signum();
                self.updates.variable_bound(variable, bound_direction).is_none()
            }).unwrap();

        // Compute the variable bound and apply the change.
        let value = (right_hand_side - total_activity) / target_coefficient;
        let bound_direction = !activity_direction * target_coefficient.signum();
        match self.updates.update_activity_variable_bound(target_column, bound_direction, value) {
            BoundChange::None => Change::None,
            BoundChange::NewBound => {
                self.after_bound_change(target_column, bound_direction, None);
                Change::Meaningful
            },
            BoundChange::BoundShift(shift) => {
                self.after_bound_change(target_column, bound_direction, Some(shift));
                Change::NotMeaningful
            },
        }
    }

    fn can_variable_rule_be_applied(&self, constraint: usize, activity_direction: BoundDirection) -> Option<OF> {
        use RangedConstraintRelation::{Equal, Range, Less, Greater};

        let right_hand_side = self.updates.b(constraint).clone();
        match self.updates.constraint_type(constraint) {
            Equal => Some(right_hand_side),
            Range(range) => match activity_direction {
                BoundDirection::Lower => Some(right_hand_side),
                BoundDirection::Upper => Some(right_hand_side - range),
            },
            Less => match activity_direction {
                BoundDirection::Lower => Some(right_hand_side),
                BoundDirection::Upper => None,
            },
            Greater => match activity_direction {
                BoundDirection::Lower => None,
                BoundDirection::Upper => Some(right_hand_side),
            },
        }
    }
}

enum ConstraintUpdate<F> {
    Remove,
    Replace(InequalityRelation, F),
    SetVariablesToBound,
}
