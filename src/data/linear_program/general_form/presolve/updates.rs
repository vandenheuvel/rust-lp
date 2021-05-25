//! # Updates
//!
//! Changes derived are stored in an `Updates` struct. Not all derived information will necessarily
//! be returned, as e.g. redundant variable bounds don't make the problem simpler.
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

use relp_num::{OrderedField, OrderedFieldRef};

use crate::data::linear_algebra::traits::SparseElement;
use crate::data::linear_program::elements::{BoundDirection, LinearProgramType, Objective, RangedConstraintRelation};
use crate::data::linear_program::general_form::{GeneralForm, RemovedVariable};
use crate::data::linear_program::general_form::presolve::counters::Counters;
use crate::data::linear_program::general_form::presolve::is_empty_constraint_feasible;

/// Collection of changes that can be applied to the `GeneralForm` without changing its solution.
///
/// This struct doesn't just collect changes, it also presents the latest "version" of the problem
/// to the presolving logic, and as such it is read from continuously.
///
/// Some value are replacements, others are changes, so be careful. The one attribute that isn't a
/// change that can be made, is the borrow of the original problem to which these change can be
/// applied. This borrow is there only to enable a few methods to be shifted down from the `Index`
/// struct to this struct, such that borrows on `self` within `Index` block less.
///
/// All indices in this struct are relative to the problem without anything removed (couldn't really
/// be different, but for clarity).
pub struct Updates<'a, F> {
    // Permanent changes
    /// New constraint values (replace the current value), by constraint index.
    pub b: HashMap<usize, F>,
    /// New constraint types.
    ///
    /// Note that these are always inequalities.
    pub constraints: HashMap<usize, RangedConstraintRelation<F>>,
    /// Change to the fixed cost.
    pub fixed_cost: F,
    /// New or tightened bounds for variables, replace old bounds.
    pub bounds: HashMap<(usize, BoundDirection), F>,
    /// Variables that were eliminated, and their solution.
    ///
    /// This solution might depend on the solution of other variables, see the documentation of the
    /// `RemovedVariable` struct. All variables for which a solution was found should be part of the
    /// `removed_variables` field.
    pub removed_variables: Vec<(usize, RemovedVariable<F>)>,
    /// Collecting the indices, relative to the active part of the problem, that should be removed.
    ///
    /// This collection is here only to avoid an extra scan over the `column_counters`.
    pub constraints_marked_removed: Vec<usize>,

    // Changes during presolving, might be discarded when computing final changes
    /// Activity variable bounds.
    ///
    /// Bounds derived using constraint activities should only be saved if they actually make the
    /// problem simpler easier to solve. For example by having been used to remove a constraint
    /// (then the item is moved to the `bounds` field during solving) or to avoid having a free
    /// variable (decided during conversion to the `Changes` struct).
    ///
    /// Bounds in this field should not appear in the other bounds field.
    ///
    /// TODO(ENHANCEMENT): When exactly should these be kept, and when not? In case of free
    ///  variables, definitely yes.
    pub activity_bounds: HashMap<(usize, BoundDirection), F>,

    // Reference for architectural convenience
    /// Original problem being solved. Included as an attribute to allow nicer code, see above.
    general_form: &'a GeneralForm<F>,
}

/// What kind of change was made to a bound.
///
/// Used within the methods to determine what needs to happen after bound was (maybe) updated.
#[derive(Eq, PartialEq)]
pub (super) enum BoundChange<F> {
    /// No change was made.
    ///
    /// There must have already been a bound that is at least as tight.
    None,
    /// The bound was added as a new bound, there was no bound before.
    NewBound,
    /// A bound already existed, but the new bound was strictly tighter, so it was updated.
    ///
    /// The field is the difference between the old bound and the new bound.
    BoundShift(F),
}

impl<'a, OF> Updates<'a, OF>
where
    OF: OrderedField + SparseElement<OF>,
    for<'r> &'r OF: OrderedFieldRef<OF>,
{
    /// Create a new instance.
    ///
    /// Requires iteration over all counters.
    ///
    /// TODO(OPTIMIZATION): Are these iterations efficiently done after compiler optimizations?
    ///
    /// # Arguments
    ///
    /// * `general_form`: Problem being presolved.
    /// * `counters`: Counters initialized on the problem in the other argument.
    ///
    /// # Return value
    ///
    /// A new instance.
    pub(super) fn new(
        general_form: &'a GeneralForm<OF>,
        counters: &Counters<OF>,
    ) -> Result<Self, LinearProgramType<OF>> {
        let mut fixed_cost = OF::zero();
        let removed_variables = counters.variable.iter().enumerate()
            .filter(|&(_, &count)| count == 0)
            .map(|(j, _)| {
                let variable = &general_form.variables[j];
                // There should be a feasible value, otherwise, reading of the problem should have
                // been ended due to a trivial infeasibility.
                debug_assert!(variable.has_feasible_value());

                let value = if variable.cost.is_zero() {
                    variable.get_feasible_value().unwrap()
                } else {
                    let value = optimize_independent_column(
                        general_form.objective,
                        &variable.cost,
                        (variable.lower_bound.as_ref(), variable.upper_bound.as_ref()),
                    )?;
                    fixed_cost += &variable.cost * &value;
                    value
                };

                Ok((j, RemovedVariable::Solved(value)))
            }).collect::<Result<_, _>>()?;

        let constraints_marked_removed = counters.constraint.iter().enumerate()
            .filter(|&(_, &count)| count == 0)
            .map(|(constraint, _)| {
                let right_hand_side = &general_form.b[constraint];
                let constraint_type = &general_form.constraint_types[constraint];
                if is_empty_constraint_feasible(right_hand_side, constraint_type) {
                    Ok(constraint)
                } else {
                    Err(LinearProgramType::Infeasible)
                }
            })
            .collect::<Result<_, _>>()?;

        Ok(Self {
            b: HashMap::new(),
            constraints: HashMap::new(),
            fixed_cost,
            bounds: HashMap::new(),
            // Variables that don't appear in any constraint are removed right away.
            removed_variables,
            constraints_marked_removed,

            activity_bounds: HashMap::new(),

            general_form,
        })
    }

    /// Get the latest constraint value.
    pub(super) fn b(&self, constraint: usize) -> &OF {
        // There is no reason to want to know the latest value for a constraint already removed.
        debug_assert!(!self.constraints_marked_removed.contains(&constraint));

        self.b.get(&constraint).unwrap_or(&self.general_form.b[constraint])
    }
    /// Change the constraint value by a certain amount.
    pub(super) fn change_b(&mut self, constraint: usize, change: OF) {
        // There is no reason to change the value for a constraint already removed.
        debug_assert!(!self.constraints_marked_removed.contains(&constraint));

        if let Some(existing_value) = self.b.get_mut(&constraint) {
            *existing_value += change;
        } else {
            self.b.insert(constraint, &self.general_form.b[constraint] + change);
        }
    }

    /// Get the latest constraint type on a constraint.
    ///
    /// When the constraint has been removed (that is, the index is in the
    /// `constraints_marked_removed` field), this function will still return some value (that is
    /// then meaningless).
    pub(super) fn constraint_type(&self, constraint: usize) -> &RangedConstraintRelation<OF> {
        // There is no reason to want to know the latest value for a constraint already removed.
        debug_assert!(!self.constraints_marked_removed.contains(&constraint));

        self.constraints.get(&constraint)
            .unwrap_or(&self.general_form.constraint_types[constraint])
    }

    /// Whether the latest version of a variable has exactly one feasible value.
    pub(super) fn is_variable_fixed(&self, variable: usize) -> Option<&OF> {
        // There is no reason to want to know the latest value for a variable already removed.
        debug_assert!(!self.removed_variables.iter().any(|&(j, _)| j == variable));

        let lower = self.variable_bound(variable, BoundDirection::Lower);
        let upper = self.variable_bound(variable, BoundDirection::Upper);
        match (lower, upper) {
            (Some(lower), Some(upper)) if lower == upper => Some(lower),
            _ => None,
        }
    }
    /// Whether a the latest version of variable has any feasible value.
    pub(super) fn variable_feasible_value(&self, variable: usize) -> Option<OF> {
        // There is no reason to want to know the latest value for a variable already removed.
        debug_assert!(!self.removed_variables.iter().any(|&(j, _)| j == variable));

        let lower = self.variable_bound(variable, BoundDirection::Lower);
        let upper = self.variable_bound(variable, BoundDirection::Upper);

        // It is often positive slacks that get removed here, so we pick the upper bound if
        // there is one, because it makes it more likely that a meaningful variable is zero,
        // which in turn might lead to sparser solutions.
        match (lower, upper) {
            (None, None) => Some(OF::zero()),
            (None, Some(bound)) => Some(bound.clone()),
            (Some(bound), None) => Some(bound.clone()),
            (Some(lower), Some(upper)) => {
                if lower <= upper {
                    Some(upper.clone())
                } else {
                    None
                }
            }
        }
    }

    /// Latest version of a variable bound.
    pub(super) fn variable_bound(&self, variable: usize, direction: BoundDirection) -> Option<&OF> {
        debug_assert!(self.removed_variables.iter().all(|&(j, _)| j != variable));
        debug_assert!({
            // Variable bounds in the `activity_bounds` should always be tighter than the bounds
            // that will be adopted for certain.
            let certain = self.bounds.get(&(variable, direction));
            let activity = self.activity_bounds.get(&(variable, direction));
            !matches!((certain, activity), (Some(_), Some(_)))
        });

        self.activity_bounds.get(&(variable, direction))
            .or_else(|| self.bounds.get(&(variable, direction)))
            .or_else(|| {
                let variable = &self.general_form.variables[variable];
                match direction {
                    BoundDirection::Lower => &variable.lower_bound,
                    BoundDirection::Upper => &variable.upper_bound,
                }.as_ref()
            })
    }

    /// Update a variable bound.
    ///
    /// This is the type of bound that should always be exported, because other changes might depend
    /// on it.
    ///
    /// # Arguments
    ///
    /// * `variable`: Index of the variable to be changed.
    /// * `direction`: Whether the lower or upper bound should be changed.
    /// * `new`: New value for the bound.
    ///
    /// # Return value
    ///
    /// Type of change made: none, a shift, or added a new bound.
    pub(super) fn update_bound(
        &mut self,
        variable: usize,
        direction: BoundDirection,
        new: OF,
    ) -> BoundChange<OF> {
        debug_assert!(self.removed_variables.iter().all(|&(j, _)| j != variable));

        // The bounds in `self.bounds` are disjoint from the bounds in `self.activity_bounds`.
        let compare_with = match self.bounds.get(&(variable, direction)) {
            None => {
                match self.activity_bounds.remove(&(variable, direction)) {
                    None => {
                        let original = &self.general_form.variables[variable];
                        let original = match direction {
                            BoundDirection::Lower => &original.lower_bound,
                            BoundDirection::Upper => &original.upper_bound,
                        };
                        match original {
                            None => {
                                self.bounds.insert((variable, direction), new);
                                return BoundChange::NewBound
                            },
                            Some(original) => original,
                        }
                    },
                    Some(existing_activity) => {
                        self.bounds.insert((variable, direction), existing_activity);
                        self.bounds.get(&(variable, direction)).unwrap()
                    }
                }
            },
            Some(existing) => existing,
        };
        // TODO(ARCHITECTURE): Avoid this clone
        let compare_with = compare_with.clone();

        // We add to the `self.bounds` struct in every case
        bound_compare_and_update(variable, direction, new, &compare_with, &mut self.bounds)
    }

    /// Update an activity variable bound.
    ///
    /// This is the type of bound that may not be exported. That happens when the derived bound has
    /// not been used to make the problem easier to solve by an algorithm.
    ///
    /// # Arguments
    ///
    /// * `variable`: Index of the variable to be changed.
    /// * `direction`: Whether the lower or upper bound should be changed.
    /// * `new`: New value for the bound.
    ///
    /// # Return value
    ///
    /// Type of change made: none, a shift, or added a new bound.
    pub(super) fn update_activity_variable_bound(
        &mut self,
        variable: usize,
        direction: BoundDirection,
        new: OF,
    ) -> BoundChange<OF> {
        debug_assert!(self.removed_variables.iter().all(|&(j, _)| j != variable));

        // TODO(ARCHITECTURE): Avoid these clones
        match self.activity_bounds.get(&(variable, direction)) {
            None => match self.bounds.get(&(variable, direction)) {
                None => {
                    let original = &self.general_form.variables[variable];
                    let bound = match direction {
                        BoundDirection::Lower => &original.lower_bound,
                        BoundDirection::Upper => &original.upper_bound,
                    };

                    match bound {
                        None => {
                            self.activity_bounds.insert((variable, direction), new);
                            BoundChange::NewBound
                        },
                        Some(original) => bound_compare_and_update(
                            variable, direction, new,
                            original,
                            &mut self.activity_bounds,
                        ),
                    }
                },
                Some(derived) => bound_compare_and_update(
                    variable, direction, new,
                    &derived.clone(),
                    // If the bound will already be adopted for sure, change it here
                    &mut self.bounds,
                ),
            },
            Some(activity_derived) => bound_compare_and_update(
                variable, direction, new,
                &activity_derived.clone(),
                &mut self.activity_bounds,
            ),
        }
    }

    /// Sets variables that can be optimized independently of all others to their optimal values.
    pub(super) fn optimize_column_independently(
        &mut self,
        variable: usize,
    ) -> Result<RemovedVariable<OF>, LinearProgramType<OF>> {
        debug_assert_ne!(self.general_form.variables[variable].cost, OF::zero());

        let new_value = self.get_optimized_value(variable)?;
        self.fixed_cost += &self.general_form.variables[variable].cost * &new_value;

        Ok(RemovedVariable::Solved(new_value))
    }

    /// Gets the value that maximizes a column's contribution in the objective direction.
    fn get_optimized_value(&self, variable: usize) -> Result<OF, LinearProgramType<OF>> {
        let objective = self.general_form.objective;
        let cost = &self.general_form.variables[variable].cost;
        let lower_bound = self.variable_bound(variable, BoundDirection::Lower);
        let upper_bound = self.variable_bound(variable, BoundDirection::Upper);

        optimize_independent_column(objective, cost, (lower_bound, upper_bound))
    }

    /// Get all changes that the presolving resulted in.
    ///
    /// This method only exists because this struct also holds a borrow to the original problem,
    /// the caller of this method.
    ///
    /// # Return value
    ///
    /// Almost all values of this struct; just not the borrow to the original problem.
    pub fn into_changes(mut self) -> Changes<OF> {
        // Remove changes from rows and columns that are deleted anyway
        for constraint in &self.constraints_marked_removed {
            self.b.remove(constraint);
            self.constraints.remove(constraint);
        }
        for variable in self.removed_variables.iter().map(|&(j, _)| j) {
            self.bounds.remove(&(variable, BoundDirection::Lower));
            self.bounds.remove(&(variable, BoundDirection::Upper));

            self.activity_bounds.remove(&(variable, BoundDirection::Lower));
            self.activity_bounds.remove(&(variable, BoundDirection::Upper));
        }

        // Keep derived variable bounds only when they help eliminate free variables
        let free_to_be_restricted = self.activity_bounds.iter()
            .map(|(&(variable, _), _)| variable)
            .filter(|&variable| {
                // Was initially free?
                self.general_form.variables[variable].is_free() &&
                    // No new bounds derived?
                    self.bounds.get(&(variable, BoundDirection::Lower)).is_none() &&
                    self.bounds.get(&(variable, BoundDirection::Upper)).is_none()
            })
            .collect::<HashSet<_>>();
        for ((variable, direction), bound_value) in self.activity_bounds {
            if free_to_be_restricted.contains(&variable) {
                self.bounds.insert((variable, direction), bound_value);
            }
        }

        // Remove changes that are not actually changes (might have been changed back and forth)
        let existing_b = &self.general_form.b;
        let b = self.b.into_iter()
            .filter(|&(i, ref new_b)| new_b != &existing_b[i])
            .collect();

        let existing_constraint_types = &self.general_form.constraint_types;
        let constraints = self.constraints.into_iter()
            .filter(|(constraint, constraint_type)| {
                &existing_constraint_types[*constraint] != constraint_type
            })
            .collect();

        self.removed_variables.sort_unstable_by_key(|&(j, _)| j);
        self.constraints_marked_removed.sort_unstable();

        Changes {
            b,
            constraints,
            fixed_cost: self.fixed_cost,
            bounds: self.bounds,
            removed_variables: self.removed_variables,
            constraints_marked_removed: self.constraints_marked_removed,
        }
    }

    pub fn nr_variables_remaining(&self) -> usize {
        self.general_form.nr_active_variables() - self.removed_variables.len()
    }

    pub fn nr_constraints_remaining(&self) -> usize {
        self.general_form.nr_active_constraints() - self.constraints_marked_removed.len()
    }
}


/// Fields in this struct are mostly the same as in the `Updates` struct.
#[derive(Eq, PartialEq, Debug)]
pub struct Changes<F> {
    /// New constraint values (replace the current value), by constraint index.
    pub b: HashMap<usize, F>,
    /// New constraint types.
    ///
    /// Not sorted.
    pub constraints: Vec<(usize, RangedConstraintRelation<F>)>,
    /// Change to the fixed cost.
    pub fixed_cost: F,
    /// New or tightened bounds for variables, replace old bounds.
    pub bounds: HashMap<(usize, BoundDirection), F>,
    /// Variables that were eliminated, and their solution.
    ///
    /// This solution might depend on the solution of other variables, see the documentation of the
    /// `RemovedVariable` struct. All variables for which a solution was found should be part of the
    /// `removed_variables` field.
    pub removed_variables: Vec<(usize, RemovedVariable<F>)>,
    /// Collecting the indices, relative to the active part of the problem, that should be removed.
    ///
    /// This collection is here only to avoid an extra scan over the `column_counters`.
    pub constraints_marked_removed: Vec<usize>,
}

/// Helper function to update bounds.
///
/// Compares with an existing bound if there is one, and updates if the new bound is tighter.
///
/// # Arguments
///
/// * `variable`: Index of the variable who's bound might be changed.
/// * `direction`: Whether an upper or lower bound might be changed.
/// * `new`: New bound value
/// * `existing`: Existing value to compare with
/// * `bounds`: Collection to add the new bound to in case it is tighter.
fn bound_compare_and_update<OF>(
    variable: usize,
    direction: BoundDirection,
    new: OF,
    existing: &OF,
    bounds: &mut HashMap<(usize, BoundDirection), OF>,
) -> BoundChange<OF>
where
    OF: OrderedField,
    for<'r> &'r OF: OrderedFieldRef<OF>,
{
    if match direction {
        BoundDirection::Lower => &new > existing,
        BoundDirection::Upper => &new < existing,
    } {
        let difference = &new - existing;
        bounds.insert((variable, direction), new);

        BoundChange::BoundShift(difference)
    } else {
        BoundChange::None
    }
}

/// Determine the optimal, or just any, feasible value for a variable.
///
/// # Arguments
///
/// * `objective`: Optimization direction.
/// * `cost`: Coefficient of the variable in the cost function.
/// * `bounds`: The lower and / or upper bound, if there are any.
///
/// # Returns
///
/// The optimal value if the variable appears in the objective function. Otherwise, any feasible
/// value.
///
/// # Errors
///
/// The problem might be unbounded.
fn optimize_independent_column<OF>(
    objective: Objective,
    cost: &OF,
    (lower_bound, upper_bound): (Option<&OF>, Option<&OF>),
) -> Result<OF, LinearProgramType<OF>>
where
    OF: OrderedField,
    for<'r> &'r OF: OrderedFieldRef<OF>,
{
    let cost_sign = cost.cmp(&OF::zero());

    match (objective, cost_sign) {
        (Objective::Minimize, Ordering::Greater) | (Objective::Maximize, Ordering::Less) => {
            lower_bound.ok_or(LinearProgramType::Unbounded)
        },
        (Objective::Maximize, Ordering::Greater) | (Objective::Minimize, Ordering::Less) => {
            upper_bound.ok_or(LinearProgramType::Unbounded)
        },
        (_, Ordering::Equal) => panic!("Should not be called if there is no cost"),
    }.map(OF::clone)
}
