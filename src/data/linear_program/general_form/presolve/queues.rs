//! # Queues
//!
//! Rules are not directly applied, but instead queued.
use std::collections::{HashSet, VecDeque};

use relp_num::{OrderedField, OrderedFieldRef};

use crate::data::linear_program::elements::BoundDirection;
use crate::data::linear_program::general_form::GeneralForm;
use crate::data::linear_program::general_form::presolve::counters::Counters;

/// Which rules still need to be applied to which constraint or variable.
///
/// Note that a value could be removed from a queue before the rule is actually applied. A variable
/// might be a slack, for example, but is then removed because a bound has been adjusted.
pub(super) struct Queues {
    /// Variables that are fixed and need substitution.
    ///
    /// A fixed variable has the upper bound equal to the lower bound.
    ///
    /// * Start of process: Column is determined fixed.
    /// * Insertion: Column becomes fixed after a bound was adjusted, and it is equal to the other
    /// bound.
    /// * State: Column is fixed. When a bound gets adjusted again, the variable is infeasible, and
    /// the solving procedure should stop, making this queue irrelevant.
    /// * Removal: Fixed value for the column is about to be substituted in the problem, effecting
    /// the right-hand side.
    /// * Causes: All rows that have a small number of elements left in them to be added to the
    /// `empty_row` or `bound` queue, perhaps the `activity` queue if sufficient variables that are
    /// left in the column have the relevant bound.
    ///
    /// This is a stack because their order of removal does not matter because they will only be in
    /// this queue only once.
    pub(crate) substitution: Vec<usize>,

    /// Constraints to check to see whether they are a bound (without a slack).
    ///
    /// * Start of process: Row count is 1.
    /// * Insertion: Row count becomes 1.
    /// * State: Row count is at most one, and row could have been removed because it became empty
    /// before it is reached in this queue.
    /// * Removal: Fixed value for the column is about to be substituted in the problem, effecting
    /// the right-hand side.
    /// * Causes: Row will be marked removed and the implied bound will be added to the updates.
    ///
    /// This is a stack because their order of removal does not matter because they will only be in
    /// this queue only once. Check whether the constraint is still active at removal.
    pub(crate) bound: Vec<usize>,

    /// Variables to check for slack variables.
    ///
    /// * Start of process: Does not appear in objective and column count is 1.
    /// * Insertion: Does not appear in objective and column count becomes 1
    /// * State: Does not appear in objective. Column count is at most 1, column could have been
    /// removed because it became empty due to the constraint being recognized as a bound (that rule
    /// has higher priority).
    /// * Removal: Just before rule is attempted to be applied.
    /// * Causes: Nothing, column gets removed, or column and row get removed.
    ///
    /// Columns can be inserted into this queue at most two times (they can only be re-added to this
    /// queue when the constraint type changes from equality to inequality, which can happen only
    /// once). Because of that, a stack is used (check whether the variable is still active at
    /// removal).
    pub(crate) slack: Vec<usize>,

    /// Constraints to check for activity bound tightening.
    ///
    /// * Start of process: Number of variable bounds missing to compute the activation bounds is at
    /// most 1.
    /// * Insertion: Number of variable bounds missing to compute the activation bounds is at
    /// reaches 1 or 0. Note that this rule might add this constraint to the queue again.
    /// * State: Number of variable bounds missing to compute the activation bounds is at most 1.
    /// The constraint might have already been removed, so checking that it's still active is
    /// necessary.
    /// * Removal: Just before rule is attempted to be applied (note that it might be added to the
    /// queue again by the same rule).
    /// * Causes: Nothing, a constraint to be removed, or a column to be removed.
    ///
    /// The relevant activity counters are `.0` for `Lower`, `.1` for `Upper`.
    /// TODO(ENHANCEMENT): Consider tracking which indices are still missing to avoid a scan when
    ///  count is 1.
    pub(crate) activity: ActivityQueue,
}

/// A queue with unique values that has quick insertion and removal of an arbitrary element.
///
/// Sets like stdlib's default `BTreeSet` do in-order traversal. Instead, we want to consider these
/// constraints that haven't been considered the longest, such that hopefully, as much as possible
/// has changed in the available information.
///
/// Insertion and removal is performant. Duplicate elements will not be added to the queue again.
pub(super) struct ActivityQueue {
    vec_deque: VecDeque<(usize, BoundDirection)>,
    set: HashSet<(usize, BoundDirection)>,
}
impl ActivityQueue {
    fn new(vec_deque: VecDeque<(usize, BoundDirection)>) -> Self {
        let set = vec_deque.iter().copied().collect();

        Self {
            vec_deque,
            set,
        }
    }
    pub(crate) fn insert(&mut self, constraint: usize, direction: BoundDirection) {
        if self.set.insert((constraint, direction)) {
            self.vec_deque.push_back((constraint, direction));
        }
    }
    pub(crate) fn pop(&mut self) -> Option<(usize, BoundDirection)> {
        let chosen_value = self.vec_deque.pop_front();
        if let Some(value) = chosen_value {
            self.set.remove(&value);
        }
        chosen_value
    }
    fn is_empty(&self) -> bool {
        self.vec_deque.is_empty()
    }
}

impl Queues {
    /// Create a new instance.
    ///
    /// Requires iteration over all counters. This method also demonstrates, which constraints or
    /// variables get added to which queue under which condition.
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
    pub fn new<OF: OrderedField>(general_form: &GeneralForm<OF>, counters: &Counters<OF>) -> Self
    where
        for<'r> &'r OF: OrderedFieldRef<OF>,
    {
        Self {
            bound: counters.constraint.iter().enumerate()
                .filter(|&(_, &count)| count == 1)
                .map(|(i, _)| i).collect(),
            activity: ActivityQueue::new(counters.activity.iter().enumerate()
                .flat_map(|(i, &(lower_count, upper_count))| {
                    let mut sides = Vec::with_capacity(2);
                    if counters.constraint[i] > 1 {
                        if lower_count <= 1 {
                            sides.push((i, BoundDirection::Lower));
                        }
                        if upper_count <= 1 {
                            sides.push((i, BoundDirection::Upper));
                        }
                    }
                    sides
                }).collect()),

            slack: counters.variable.iter().enumerate()
                .filter(|&(_, &count)| count == 1)
                .filter(|&(j, _)| general_form.variables[j].cost.is_zero())
                .map(|(j, _)| j).collect(),
            substitution: general_form.variables.iter().enumerate()
                // Columns with no values get substituted right away, as they removal doesn't lead
                // to any further actions. The only values that should be put in the queue are those
                // that lead to other actions, that is, have a nonzero counter.
                .filter(|&(j, _)| counters.variable[j] > 0)
                .filter_map(|(j, variable)| variable.is_fixed().map(|_| j))
                .collect(),
        }
    }

    /// Whether all queues are empty.
    ///
    /// This indicates whether the repeated application of reduction rules can be stopped.
    pub(crate) fn are_empty(&self) -> bool {
        // Note the reverse order w.r.t. the order in which these queues are tested in the main loop
        self.activity.is_empty()
            && self.slack.is_empty()
            && self.bound.is_empty()
            && self.substitution.is_empty()
    }
}
