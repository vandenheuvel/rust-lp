//! # Linear programs in "general form"
//!
//! Data structure for manipulation of linear programs.
use std::cmp::Ordering;
use std::collections::HashMap;
use std::iter::Sum;
use std::mem;
use std::ops::{Add, AddAssign, Mul, Neg, Sub};

use daggy::{Dag, WouldCycle};
use daggy::petgraph::data::Element;
use itertools::repeat_n;
use num::Zero;

use crate::algorithm::two_phase::matrix_provider::matrix_data::MatrixData;
use crate::algorithm::utilities::remove_indices;
use crate::data::linear_algebra::matrix::{ColumnMajor, Order};
use crate::data::linear_algebra::matrix::Sparse;
use crate::data::linear_algebra::SparseTupleVec;
use crate::data::linear_algebra::traits::{SparseComparator, SparseElement as LinearAlgebraElement};
use crate::data::linear_algebra::vector::{DenseVector, SparseVector, Vector};
use crate::data::linear_program::elements::{BoundDirection, LinearProgramType, Objective, RangedConstraintRelation, RangedConstraintRelationKind, VariableType};
use crate::data::linear_program::general_form::presolve::{Change, Index as PresolveIndex};
use crate::data::linear_program::general_form::presolve::scale::Scalable;
use crate::data::linear_program::general_form::presolve::updates::Changes;
use crate::data::linear_program::general_form::RemovedVariable::{FunctionOfOthers, Solved};
use crate::data::linear_program::solution::Solution;
use crate::data::number_types::traits::{OrderedField, OrderedFieldRef};

mod presolve;

/// A linear program in general form.
///
/// This structure is used as a first storage independent representation format for different
/// parse results to be transformed to.
///
/// Can be checked for consistency by the `is_consistent` method in this module. That method can be
/// viewed as documentation for the requirements on the variables in this data structure.
#[derive(Debug, Eq, PartialEq)]
pub struct GeneralForm<F> {
    /// Which direction does the objective function go?
    objective: Objective,

    /// Constant in the cost function.
    fixed_cost: F,

    // Constraint related
    /// All constraint coefficients.
    ///
    /// Has size `constraint_types.len()` in the row direction, size `variables.len()` in the column
    /// direction.
    constraints: Sparse<F, F, ColumnMajor>,
    /// The equation type of all rows, ordered by index.
    ///
    /// These are read "from constraint to constraint value", meaning:
    /// * When a constraint is `ConstraintType::Equal`, the equation is `<a, x> == b`
    /// * When a constraint is `ConstraintType::Range(r)`, the equation is `b - r <= <a, x> <= b`
    /// * When a constraint is `ConstraintType::Less`, the equation is `<a, x> <= b`
    /// * When a constraint is `ConstraintType::Greater`, the equation is `<a, x> >= b`
    constraint_types: Vec<RangedConstraintRelation<F>>,
    /// All right-hands sides of equations.
    b: DenseVector<F>,

    // Variable related
    /// Information about all *active* variables, that is, variables that are not yet presolved.
    variables: Vec<Variable<F>>,
    /// For all variables, presolved or not, a placeholder with a potential solution or method to
    /// derive the solution once the other active variables are solved.
    ///
    /// It can be that a single original variable is represented in two active variables because it
    /// is free. See the `OriginalVariable::ActiveFree` variant for more info.
    original_variables: Vec<(String, OriginalVariable<F>)>,
    /// Mapping indices of unsolved variables to their index in the original problem.
    ///
    /// Generally, only a single active variable points to an original variable. In the special case
    /// of a free original variable, two active variables might point to the original variable.
    from_active_to_original: Vec<usize>,
}

/// Whether a variable from the original problem has been eliminated.
#[derive(Debug, Eq, PartialEq)]
enum OriginalVariable<F> {
    /// The variable is still active.
    ///
    /// Either presolve has not yet been attempted, or was unsuccessful in eliminating the variable
    /// from the problem.
    ///
    /// Contains the index of the variable in the active part of the problem. This index is used to
    /// index into `GeneralForm.variables`.
    Active(usize),

    /// The variable free, active and split up in two bounded variables.
    ///
    /// If a variable is free (not bounded below, not bounded above), it is unclear how to directly
    /// treat it with the simplex method. So we split the variable up: `x = x+ - x-` with both
    /// `x+ >= 0` and `x- >= 0`. The first field of this variant refers to the index of `x+`, the
    /// latter to the index of `x-`.
    ActiveFree(usize, usize),

    /// Variable was removed.
    ///
    /// Probably by a presolve operation. An explicit value might not be known.
    Removed(RemovedVariable<F>),
}

/// A variable from the original problem that was removed.
#[derive(Debug, Eq, PartialEq)]
pub enum RemovedVariable<F> {
    /// Variable was determined to an explicit value.
    Solved(F),
    /// Variable was determined as a function of other variables.
    ///
    /// Affine function of the form `constant - <coefficients, x>` where some of the `x` might be
    /// unknown (or at least not explicitly known) at this point in the solution process.
    // TODO(PERFORMANCE): Storing the coefficient of this slack instead of normalizing the
    //  constraint before saving it here.
    FunctionOfOthers {
        /// Value of the constraint the slack was removed from.
        ///
        /// This value was normalized by the coefficient the slack variable had in the objective
        /// function.
        ///
        /// The slack will be chosen such that the equation, as it looked like when this variable
        /// was removed, "activates" to this value (up to normalization).
        constant: F,
        /// Coefficients of the constraint the slack was removed from.
        ///
        /// These coefficients were normalized by the coefficient of the slack in the constraint.
        coefficients: SparseTupleVec<F>,
    },
}

/// Check whether the dimensions of the `GeneralForm` are consistent.
///
/// This method might be expensive, use it in debugging only. It can be viewed as a piece of
/// documentation on the requirements of a `GeneralForm` struct.
fn is_consistent<F: LinearAlgebraElement<F> + SparseComparator>(general_form: &GeneralForm<F>) -> bool {
    // Reference values
    let nr_active_constraints = general_form.nr_active_constraints();
    let nr_active_variables = general_form.nr_active_variables();
    let nr_free_variables = general_form.original_variables.iter().filter(|(_, variable)| {
        matches!(variable, OriginalVariable::ActiveFree(_, _))
    }).count();

    let b = general_form.b.len() == nr_active_constraints;
    let constraints = general_form.constraint_types.len() == nr_active_constraints;
    let rows = general_form.constraints.nr_rows() == nr_active_constraints;

    let variables = general_form.variables.len() == nr_active_variables;
    let columns = general_form.constraints.nr_columns() == nr_active_variables;
    let original_variables = {
        let nr_original_variables = general_form.original_variables.len();
        let size = nr_original_variables >= nr_active_variables - nr_free_variables;
        let kept_increasing = general_form.original_variables.iter()
            .filter_map(|(_, variable)| match variable {
                OriginalVariable::Active(index) => Some(*index),
                OriginalVariable::ActiveFree(positive, _) => Some(*positive),
                _ => None,
            }).is_sorted();
        let no_cycles = {
            let nodes = repeat_n(Element::Node { weight: (), }, nr_original_variables);
            let edges = general_form.original_variables
                .iter().enumerate()
                .filter_map(|(target, (_, variable))| match variable {
                    OriginalVariable::Removed(FunctionOfOthers{ constant: _, coefficients: els, }) => Some((els, target)),
                    _ => None,
                })
                .flat_map(|(els, target)| {
                    els.iter().map(move |&(j, _)| Element::<(), ()>::Edge { source: j, target, weight: (), })
                });

            !matches!(Dag::<(), (), usize>::from_elements(nodes.chain(edges)), Err(WouldCycle(_)))
        };

        [size, kept_increasing, no_cycles].iter().all(|v| *v)
    };
    let from_active_to_original = {
        let size = general_form.from_active_to_original.len() == nr_active_variables;
        let sorted_non_free = general_form.from_active_to_original[..(nr_active_variables - nr_free_variables)].is_sorted();
        let sorted_free = general_form.from_active_to_original[(nr_active_variables - nr_free_variables)..].is_sorted();
        let max = if nr_active_variables > 0 {
            general_form.from_active_to_original[nr_active_variables - 1] < general_form.nr_original_variables()
        } else { true };

        [
            size,
            sorted_non_free,
            sorted_free,
            max,
        ].iter().all(|v| *v)
    };

    [
        b,
        constraints,
        rows,
        variables,
        columns,
        original_variables,
        from_active_to_original,
    ].iter().all(|v| *v)
}


impl<F: LinearAlgebraElement<F> + SparseComparator> GeneralForm<F> {
    /// Create a new linear program in general form.
    ///
    /// Simple constructor except for two indices that get created.
    pub fn new(
        objective: Objective,
        constraints: Sparse<F, F, ColumnMajor>,
        constraint_types: Vec<RangedConstraintRelation<F>>,
        b: DenseVector<F>,
        variables: Vec<Variable<F>>,
        variable_names: Vec<String>,
        fixed_cost: F,
    ) -> Self {
        let nr_active_variables = variables.len();

        let general_form = Self {
            objective,
            constraints,
            constraint_types,
            b,
            fixed_cost,
            variables,
            original_variables: variable_names.into_iter().enumerate()
                .map(|(j, name)| (name, OriginalVariable::Active(j))).collect(),
            from_active_to_original: (0..nr_active_variables).collect(),
        };

        debug_assert!(is_consistent(&general_form));

        general_form
    }
}

impl<OF: 'static> GeneralForm<OF>
where
    OF: OrderedField + LinearAlgebraElement<OF>,
    for<'r> &'r OF: OrderedFieldRef<OF>,
{
    /// Modify this linear problem such that it is representable by a `MatrixData` structure.
    ///
    /// The problem gets transformed into standard form, which also includes a presolve operation.
    /// Note that this call might be expensive.
    ///
    /// TODO(ENHANCEMENT): Make sure that presolving can be skipped.
    ///
    /// See also the documentation of the `GeneralForm::standardize` method.
    ///
    /// # Return value
    ///
    /// A `Result` containing either the `MatrixData` form of the presolved and standardized
    /// problem.
    ///
    /// # Errors
    ///
    /// In case the linear program gets solved during this presolve operation, a solution.
    pub fn derive_matrix_data(&mut self) -> Result<MatrixData<OF>, LinearProgramType<OF>> {
        self.standardize()?;

        let (
            nr_equality,
            nr_range,
            nr_upper,
            nr_lower,
        ) = self.reorder_constraints_by_type();

        let ranges = self.constraint_types[nr_equality..(nr_equality + nr_range)].iter()
            .map(|constraint_type| match constraint_type {
                RangedConstraintRelation::Range(r) => r,
                _ => unreachable!("Constraint types were just sorted, and there should only be ranges here."),
            }).collect();

        Ok(MatrixData::new(
            &self.constraints,
            &self.b,
            ranges,
            nr_equality,
            nr_range,
            nr_upper,
            nr_lower,
            &self.variables,
        ))
    }

    /// Convert this `GeneralForm` problem to a form closer to the standard form representation.
    ///
    /// This involves:
    ///
    /// * Determining which rows and columns can be removed or differently represented to reduce the
    /// problem size and increase the reading speed from a `MatrixData` structure.
    /// * Determining which variables are implicitly fixed, such that their only feasible value can
    /// be substituted into the problem and the column eliminated.
    /// * Modifying variables such that they are either free or bounded below by zero (with possibly
    /// an upper bound).
    /// * Multiplying some rows such that the constraint value is non-negative.
    ///
    /// To do the above, a column major representation of the constraint data is built. This
    /// requires copying all constraint data once.
    ///
    /// TODO(ENHANCEMENT): Make sure that presolving can be skipped.
    ///
    /// # Errors
    ///
    /// In case the linear program gets solved during this presolve operation, a solution.
    pub fn standardize(&mut self) -> Result<(), LinearProgramType<OF>> {
        self.presolve()?;
        self.transform_variables();
        self.make_b_non_negative();
        self.make_minimization_problem();
        self.scale();

        Ok(())
    }

    /// Recursively analyse constraints and variable bounds and eliminating or tightning these.
    ///
    /// In order to make the linear program easier to solve, a set of rules is applied. These rules
    /// are cheaper than the full simplex algorithm and are aimed at making the program easier to
    /// solve by other algorithms.
    ///
    /// A set of queues containing constraint and variable indices are maintained. Each presolve
    /// step attempts to apply a presolve rule to either a constraint or variable, as indicated by
    /// these indices. After a rules is applied and a change occurred, relevant constraint or
    /// variable indices might be added to queues, because a rule might be applicable.
    ///
    /// TODO(ENHANCEMENT): Normalization for numerical stability of floating point types.
    ///
    /// # Return value
    ///
    /// If the linear program gets solved during this presolve operation, a `Result::Err` return
    /// value containing the solution.
    pub(crate) fn presolve(&mut self) -> Result<(), LinearProgramType<OF>> {
        let Changes {
            b,
            constraints,
            fixed_cost,
            bounds,
            removed_variables,
            constraints_marked_removed,
        } = self.compute_presolve_changes()?;

        let variable_indices_only = removed_variables.iter().map(|&(j, _)| j).collect();
        self.update_values_that_remain(
            b,
            constraints,
            fixed_cost,
            bounds,
            removed_variables,
        );
        self.remove_rows_and_columns(constraints_marked_removed, variable_indices_only);
        debug_assert!(is_consistent(&self));

        self.compute_solution_where_possible();
        self.get_solution()
            .map_or(Ok(()), |solution| Err(LinearProgramType::FiniteOptimum(solution)))
    }

    /// Run the presolve operation using a `PresolveIndex` and return all the proposed reductions.
    ///
    /// # Return value
    ///
    /// Updates for almost all fields, see method body.
    ///
    /// # Errors
    ///
    /// In case the linear program gets solved during this presolve operation, a solution.
    fn compute_presolve_changes(
        &self
    ) -> Result<Changes<OF>, LinearProgramType<OF>> {
        let mut index = PresolveIndex::new(&self)?;

        let mut iterations_without_meaningful_change = 0;
        while !index.are_queues_empty() && iterations_without_meaningful_change < {
            // It can happen that a linear program is such that it can be recursively "tightened" by
            // presolving, arbitrarily precise, with the process never terminating. As this happens,
            // the numbers in the problem will have increasingly large representations.
            //
            // To make sure that doesn't happen, this heuristic rule decides when to stop.
            // TODO(ENHANCEMENT): When should we stop presolving?
            index.updates.nr_variables_remaining() + index.updates.nr_constraints_remaining()
        } {
            match index.presolve_step()? {
                Change::Meaningful => iterations_without_meaningful_change = 0,
                Change::None => {},
                Change::NotMeaningful => iterations_without_meaningful_change += 1,
            }
        }

        Ok(index.updates.into_changes())
    }

    /// Modify the values that remain after the presolve operation.
    ///
    /// All indices are relative to the data as it looks before the presolving started. That is,
    /// this method should be called before any values in this problem are deleted.
    ///
    /// # Arguments
    ///
    /// * `b`: Map from constraint indices to their new values.
    /// * `constraints`: Map from constraint indices to the new constraint type.
    /// * `fixed_cost`: Change to the fixed cost.
    /// * `bounds`: New or tightened bounds, by value (as opposed to change).
    /// * `solved_variables`: Values to remove from the problem.
    fn update_values_that_remain(
        &mut self,
        b: HashMap<usize, OF>,
        constraints: Vec<(usize, RangedConstraintRelation<OF>)>,
        fixed_cost: OF,
        bounds: HashMap<(usize, BoundDirection), OF>,
        solved_variables: Vec<(usize, RemovedVariable<OF>)>,
    ) {
        for (i, change) in b {
            self.b[i] = change;
        }

        for (i, constraint) in constraints {
            self.constraint_types[i] = constraint;
        }

        self.fixed_cost += fixed_cost;

        for (j, variable) in solved_variables {
            self.original_variables[j].1 = OriginalVariable::Removed(variable);
        }

        for ((j, direction), value) in bounds {
            let variable = &mut self.variables[j];
            *match direction {
                BoundDirection::Lower => &mut variable.lower_bound,
                BoundDirection::Upper => &mut variable.upper_bound,
            } = Some(value);
        }
    }

    /// Remove a set of rows and columns from the constraint data.
    ///
    /// Constraints might have been determined redundant, or perhaps they represented a variable
    /// bound. Those can also be represented in the `self.variables` property of `GeneralForm`. This
    /// method is used to clean up those rows and indices that are left after the removal of those
    /// constraints and / or variables.
    ///
    /// Note that this method is somewhat expensive because it removes columns from a row-major
    /// ordered matrix.
    fn remove_rows_and_columns(&mut self, mut constraints: Vec<usize>, mut variables: Vec<usize>) {
        constraints.sort_unstable();
        variables.sort_unstable();

        self.constraints.remove_columns(&variables);
        remove_indices(&mut self.variables, &variables);
        if !variables.is_empty() {
            // Update the `from_active_to_original` map.
            let mut skipped = 1;
            let new_length = self.from_active_to_original.len() - variables.len();
            for j in variables[0]..new_length {
                while skipped < variables.len() && variables[skipped] == j + skipped {
                    skipped += 1;
                }
                self.from_active_to_original[j] = self.from_active_to_original[j + skipped];
            }
            self.from_active_to_original.drain(new_length..);

            // Update the reverse map
            for (new_index, &variable) in self.from_active_to_original.iter().enumerate() {
                match &mut self.original_variables[variable].1 {
                    OriginalVariable::Active(index) => *index = new_index,
                    OriginalVariable::ActiveFree(_, _) => unimplemented!(),
                    _ => panic!("Should still be in the problem."),
                }
            }
        }

        self.constraints.remove_rows_although_this_matrix_is_column_major(&constraints);
        remove_indices(&mut self.constraint_types, &constraints);
        self.b.remove_indices(&constraints);

        debug_assert!(is_consistent(&self));
    }

    /// Split free variables and shift all variables, such that the lower bound is zero.
    ///
    /// This allows the removal of those lower bounds afterwards; this lower bound is the only lower
    /// bound for problems in standard form. When working with a simplex tableau, this form allows
    /// us to eliminate all rows which describe a lower bound.
    ///
    /// If later on, such as during branch and bound, an extra lower bound needs to be inserted,
    /// this information can be stored regardless in a separate data structure.
    fn transform_variables(&mut self) {
        // After this call, all active variables have at least one bound.
        self.split_free_variables();
        debug_assert!(is_consistent(&self));

        for j in 0..self.variables.len() {
            let variable = &mut self.variables[j];

            // Flip such that there is at least a lower bound.
            if let (None, Some(upper)) = (&variable.lower_bound, &variable.upper_bound) {
                variable.flipped = !variable.flipped;
                variable.shift *= -OF::one();
                variable.cost *= -OF::one();

                variable.lower_bound = Some(-upper);
                variable.upper_bound = None;

                for (_, coefficient) in &mut self.constraints.data[j] {
                    *coefficient *= -OF::one();
                }
            };

            // Shift such that any lower bound is zero.
            if let Some(ref mut lower) = variable.lower_bound {
                variable.shift -= &*lower;
                if let Some(ref mut upper) = variable.upper_bound {
                    *upper -= &*lower;
                }

                self.fixed_cost += &*lower * &variable.cost;
                // Update the bounds
                for (i, coefficient) in &self.constraints.data[j] {
                    self.b[*i] -= coefficient * &*lower;
                }

                *lower = OF::zero();
            }
        }

        debug_assert!(is_consistent(&self));
    }

    /// Split all free variables.
    ///
    /// This is done to be able to call the simplex algorithm, which can't easily handle free
    /// variables directly.
    ///
    /// See also the documentation on the `OriginalVariable::ActiveFree` variant.
    fn split_free_variables(&mut self) {
        // Select all free variables and "negative-duplicate" the column.
        let (indices_duplicated, columns) = (0..self.nr_variables())
            .filter(|&j| {
                let variable = &self.variables[j];
                variable.lower_bound.is_none() && variable.upper_bound.is_none()
            })
            .map(|j| (j, self.constraints.iter_column(j).map(|(i, v)| (*i, -v)).collect()))
            .unzip::<_, _, Vec<_>, _>();

        // Concatenate the existing constraints with the newly derived columns.
        let as_matrix = ColumnMajor::new(columns, self.nr_constraints(), indices_duplicated.len());
        self.constraints = mem::replace(
                &mut self.constraints,
                ColumnMajor::new(Vec::with_capacity(0), 0, 0),
            )
            .concatenate_horizontally(as_matrix);

        // Create the new active variables and make sure the old variable is non negative.
        for index in indices_duplicated {
            let original_index = self.from_active_to_original[index];
            self.original_variables[original_index].1 = OriginalVariable::ActiveFree(index, self.from_active_to_original.len());
            self.from_active_to_original.push(original_index);
            self.variables.push(Variable {
                variable_type: self.variables[index].variable_type,
                cost: -&self.variables[index].cost,
                lower_bound: Some(OF::zero()),
                upper_bound: None,
                shift: OF::zero(),
                flipped: false,
            });
            self.variables[index].lower_bound = Some(OF::zero());
        }
    }

    /// Multiply the constraints by a constant such that the constraint value is >= 0.
    ///
    /// This is a step towards representing a `GeneralForm` problem in standard form.
    fn make_b_non_negative(&mut self) {
        let rows_to_negate = self.b.iter_values().enumerate()
            .filter(|&(_, v)| v < &OF::zero())
            .map(|(i, _)| i)
            .collect();

        self.constraints.change_row_signs(&rows_to_negate);
        for row in rows_to_negate {
            match self.constraint_types[row] {
                RangedConstraintRelation::Less => {
                    self.constraint_types[row] = RangedConstraintRelation::Greater;
                    self.b[row] *= -OF::one();
                },
                RangedConstraintRelation::Equal => {
                    self.b[row] *= -OF::one();
                },
                RangedConstraintRelation::Greater => {
                    self.constraint_types[row] = RangedConstraintRelation::Less;
                    self.b[row] *= -OF::one();
                },
                RangedConstraintRelation::Range(ref range) => {
                    let minus_old_lower_bound = range - &self.b[row];
                    self.b[row] = minus_old_lower_bound;
                },
            };
        }

        debug_assert!(is_consistent(&self));
    }

    /// Make this a minimization problem by multiplying the cost function by `-1`.
    fn make_minimization_problem(&mut self) {
        if self.objective == Objective::Maximize {
            self.objective = Objective::Minimize;

            for variable in &mut self.variables {
                variable.cost *= -OF::one();
            }
        }
    }

    /// Split the constraints out per type.
    ///
    /// The constraints in a `GeneralForm` linear program are mixed; the of the constraint is saved
    /// in `self.constraint_types`. A standard form linear program has a separate data structure
    /// for each constraint type. This to facilitate the easy creation of a `MatrixData` data
    /// struct, which "simulates" the presence of slack variables based on those different
    /// constraint types.
    ///
    /// The constraints are sorted as:
    ///
    /// 1. Equality,
    /// 2. Range,
    /// 3. Upper bounded,
    /// 4. Lower bounded.
    ///
    /// # Return value
    ///
    /// The number of equality, range, upper and lower bounds.
    fn reorder_constraints_by_type(&mut self) -> (usize, usize, usize, usize) {
        let (mut e_counter, mut r_counter, mut l_counter, mut g_counter) = (0, 0, 0, 0);
        let map = self.constraint_types.iter().map(|constraint_type| {
            match constraint_type {
                RangedConstraintRelation::Equal => { e_counter += 1; e_counter - 1 },
                RangedConstraintRelation::Range(_) => { r_counter += 1; r_counter - 1 }
                RangedConstraintRelation::Less => { l_counter += 1; l_counter - 1 },
                RangedConstraintRelation::Greater => { g_counter += 1; g_counter - 1 },
            }
        }).collect::<Vec<usize>>();

        let old_constraint_types = self.constraint_types.iter()
            .map(From::from)
            .collect::<Vec<RangedConstraintRelationKind>>();
        self.constraint_types.sort_by(|a, b| {
            // Ordering is: ==  <  =r= < <=  <  >=
            // Sorting needs to be stable for the range constraints.
            match a {
                RangedConstraintRelation::Equal => match b {
                    RangedConstraintRelation::Equal => Ordering::Equal,
                    _ => Ordering::Less,
                },
                RangedConstraintRelation::Range(_) => match b {
                    RangedConstraintRelation::Range(_) => Ordering::Equal,
                    RangedConstraintRelation::Equal => Ordering::Greater,
                    _ => Ordering::Less,
                },
                RangedConstraintRelation::Less => match b {
                    RangedConstraintRelation::Less => Ordering::Equal,
                    RangedConstraintRelation::Greater => Ordering::Less,
                    _ => Ordering::Greater,
                },
                RangedConstraintRelation::Greater => match b {
                    RangedConstraintRelation::Greater => Ordering::Equal,
                    _ => Ordering::Greater,
                },
            }
        });

        let get_destination = |source| {
            match old_constraint_types[source] {
                RangedConstraintRelationKind::Equal => map[source],
                RangedConstraintRelationKind::Range => e_counter + map[source],
                RangedConstraintRelationKind::Less => e_counter + r_counter + map[source],
                RangedConstraintRelationKind::Greater => e_counter + r_counter + l_counter + map[source],
            }
        };

        let mut new_b = vec![None; self.b.len()];
        for i in 0..self.b.len() {
            // TODO(ARCHITECTURE): Avoid this mem::replace
            new_b[get_destination(i)] = Some(mem::replace(&mut self.b[i], OF::zero()));
        }
        self.b = DenseVector::new(new_b.into_iter().collect::<Option<Vec<_>>>().unwrap(), self.b.len());

        for column in &mut self.constraints.data {
            for (i, _) in column.iter_mut() {
                *i = get_destination(*i);
            }
            column.sort_unstable_by_key(|&(i, _)| i);
        }

        debug_assert!(is_consistent(&self));

        (e_counter, r_counter, l_counter, g_counter)
    }

    /// Get the known solution value for a variable, if there is one.
    ///
    /// # Arguments
    ///
    /// * `variable`: Index of the variable to get the solution for with respect to the *original,
    /// non-presolved* problem.
    ///
    /// # Return value
    ///
    /// `None` if the variable is still in the problem, `Some` if not. In the latter case, if the
    /// solution is a specific value known at this point in time, it contains another `Some` with
    /// value and `None` otherwise.
    ///
    /// ## Note
    ///
    /// This is only an actual solution for the variable if the problem turns out to be feasible.
    fn is_variable_presolved(&self, variable: usize) -> bool {
        debug_assert!(variable < self.variables.len());

        match self.original_variables[variable].1 {
            OriginalVariable::Active(_) => false,
            OriginalVariable::ActiveFree(_, _) => false,
            OriginalVariable::Removed(RemovedVariable::Solved(_)) => true,
            OriginalVariable::Removed(RemovedVariable::FunctionOfOthers { .. }) => true,
        }
    }

    /// Compute explicit solution from slack variables where possible.
    fn compute_solution_where_possible(&mut self) {
        // To avoid recomputations, we store all computed intermediate values in this collection.
        let mut new_solutions = vec![None; self.original_variables.len()];
        let mut changed = Vec::new();
        for (j, (_, variable)) in self.original_variables.iter().enumerate() {
            if let OriginalVariable::Removed(FunctionOfOthers { .. }) = variable {
                if self.compute_solution_value(j, &mut new_solutions).is_some() {
                    changed.push(j);
                }
            }
        }

        // Save all values
        for (variable, solution) in new_solutions.into_iter().enumerate() {
            if let Some(value) = solution {
                self.original_variables[variable].1 = OriginalVariable::Removed(Solved(value));
            }
        }

        // TODO: Include updating of fixed cost in case a presolved variable as a function of others
        //  had a nonzero cost.
    }

    /// Compute the solution value for a single variable.
    ///
    /// This is a helper method of the `compute_solution_where_possible` function.
    ///
    /// # Arguments
    ///
    /// * `variable`: Index of the variable for which we try to determine an explicit solution.
    /// * `new_solutions`: Collection containing solutions previously computed, that the computed
    /// solution value will also be written in to (if it is determined).
    ///
    /// # Return value
    ///
    /// The solution value, if it could be determined.
    fn compute_solution_value<'a>(
        &'a self,
        variable: usize,
        new_solutions: &'a mut Vec<Option<OF>>,
    ) -> Option<&'a OF> {
        match &self.original_variables[variable].1 {
            OriginalVariable::Active(_) => None,
            OriginalVariable::ActiveFree(_, _) => None,
            OriginalVariable::Removed(Solved(value)) => Some(value),
            OriginalVariable::Removed(FunctionOfOthers { constant, coefficients }) => {
                if new_solutions[variable].is_none() {
                    new_solutions[variable] = coefficients.iter().map(|(j, coefficient)| {
                        self.compute_solution_value(*j, new_solutions).map(|v| coefficient * v)
                    })
                        .sum::<Option<OF>>()
                        .map(|inner_product| constant - inner_product);
                }

                new_solutions[variable].as_ref()
            }
        }
    }

    /// If this problem is fully solved (probably by presolving), get the solution.
    ///
    /// All original variables need to have an explicit solution. Hint: Try calling
    /// `compute_solution_where_possible`.
    ///
    /// # Return value
    ///
    /// If one of the variables is still undetermined, `None`.
    pub fn get_solution(&self) -> Option<Solution<OF>> {
        let maybe_variable_values = self.original_variables.iter().map(|(name, variable)| {
            if let OriginalVariable::Removed(Solved(value)) = variable {
                Some((name.clone(), value.clone()))
            } else {
                None
            }
        }).collect::<Option<Vec<_>>>();
        maybe_variable_values.map(|variable_values| {
            Solution::new(self.fixed_cost.clone(), variable_values)
        })
    }

    /// Transform a solution back.
    ///
    /// When converting the problem to standard form, and during presolving, transformations to
    /// variables are applied that need to be undone once a solution was computed.
    ///
    /// # Arguments
    ///
    /// * `reduced_solution`: Solution values for all variables that are still marked as `Active` in
    /// the "original variables" of this problem.
    pub fn reshift_solution<G>(&self, reduced_solution: &mut SparseVector<G, G>)
    where
        G: LinearAlgebraElement<G> + SparseComparator + PartialEq<OF> + AddAssign<OF> + From<OF>,
        for<'r> &'r OF: Neg<Output=OF>,
        for<'r> &'r G: Neg<Output=G>,
    {
        debug_assert_eq!(reduced_solution.len(), self.variables.len());

        for (j, variable) in self.variables.iter().enumerate() {
            reduced_solution.shift_value::<OF>(j, -&variable.shift);
            if variable.flipped {
                if reduced_solution.get(j).is_some() {
                    let flipped = -reduced_solution.get(j).unwrap();
                    reduced_solution.set(j, flipped);
                }
            }
        }
    }

    /// Extend a reduced solution with already known values.
    ///
    /// Presolving might have determined variables explicitly or as a function of others. This
    /// method augments the solution of the open part of the problem to a complete solution.
    ///
    /// # Arguments
    ///
    /// * `reduced_solution`: A value for each of the variables in this general form that are still
    /// marked as unsolved (in the `original_variables` field).
    ///
    /// # Returns
    ///
    /// A complete solution.
    pub fn compute_full_solution_with_reduced_solution<G>(
        self,
        mut reduced_solution: SparseVector<G, G>,
    ) -> Solution<G>
    where
        // TODO: Find a suitable trait alias to avoid this many trait bounds.
        G: Sum + Neg<Output=G> + AddAssign<OF> + PartialEq<OF> + LinearAlgebraElement<G> + From<OF> + Eq + Ord + Zero,
        for<'r> G: Add<&'r OF, Output=G>,
        for<'r> &'r G: Neg<Output=G> + Mul<&'r OF, Output=G> + Add<&'r OF, Output=G> + Sub<&'r G, Output=G>,
    {
        debug_assert_eq!(reduced_solution.len(), self.variables.len());

        let cost = reduced_solution.iter_values()
            .map(|(j, v)| v * &self.variables[*j].cost)
            .sum::<G>() + &self.fixed_cost;
        self.reshift_solution(&mut reduced_solution);

        // To avoid recomputations, we store all computed intermediate values in this collection.
        let mut new_solutions = vec![None; self.original_variables.len()];
        for j in 0..self.original_variables.len() {
            self.compute_solution_value_with_bfs::<G>(j, &mut new_solutions, &reduced_solution);
        }
        debug_assert!(new_solutions.iter().all(Option::is_some));

        Solution::new(
            cost,
            self.original_variables.into_iter().zip(new_solutions.into_iter())
                .map(|((name, _), value)| (name, value.unwrap()))
                .collect(),
        )
    }

    /// Compute the solution value for a single variable.
    ///
    /// This is a helper method of the `compute_full_solution_with_reduced_solution` function.
    ///
    /// # Arguments
    ///
    /// * `variable`: Index of the variable for which we try to determine an explicit solution.
    /// Relative to the original problem with variables solved and unsolved.
    /// * `new_solutions`: Collection containing solutions previously computed, that the computed
    /// solution value will also be written in to (if it is determined).
    /// * `reduced_solution`: Solution to the reduced problem (for all variables still active).
    ///
    /// # Return value
    ///
    /// The solution value for that variable. Can always be determined if the general form is
    /// consistent (that is, there are no cycles in the original variables, see the
    /// `is_consistent` function).
    fn compute_solution_value_with_bfs<'a, G>(
        &self,
        variable: usize,
        new_solutions: &'a mut Vec<Option<G>>,
        reduced_solution: &SparseVector<G, G>,
    ) -> &'a G
    where
        G: LinearAlgebraElement<G> + Zero + From<OF> + Sum + Neg<Output=G>,
        for<'r> G: Add<&'r OF, Output=G>,
        for<'r> &'r G: Neg<Output=G> + Mul<&'r OF, Output=G> + Add<&'r OF, Output=G> + Sub<&'r G, Output=G>,
    {
        debug_assert!(variable < new_solutions.len());
        debug_assert_eq!(new_solutions.len(), self.original_variables.len());
        debug_assert_eq!(reduced_solution.len(), self.variables.len());

        if new_solutions[variable].is_some() {
            return new_solutions[variable].as_ref().unwrap();
        }

        let new_solution = match &self.original_variables[variable].1 {
            OriginalVariable::Active(j) => match reduced_solution.get(*j) {
                Some(value) => value.clone(),
                None => G::zero(),
            },
            OriginalVariable::ActiveFree(j_positive, j_negative) => {
                match (reduced_solution.get(*j_positive), reduced_solution.get(*j_negative)) {
                    (None, None) => G::zero(),
                    (Some(positive), None) => positive.clone(),
                    (None, Some(negative)) => -negative,
                    (Some(positive), Some(negative)) => positive - negative,
                }
            },
            OriginalVariable::Removed(Solved(ref v)) => v.clone().into(),
            OriginalVariable::Removed(FunctionOfOthers { constant, coefficients }) => {
                -coefficients.iter()
                    .map(|(j, coefficient)| {
                        self.compute_solution_value_with_bfs::<G>(*j, new_solutions, reduced_solution) * coefficient
                    })
                    .sum::<G>() + constant
            }
        };

        new_solutions[variable] = Some(new_solution);
        new_solutions[variable].as_ref().unwrap()
    }
}

impl<OF> GeneralForm<OF>
where
    OF: LinearAlgebraElement<OF> + SparseComparator,
{
    /// Number of constraints that have not been eliminated after a presolving operation.
    /// 
    /// During a presolving operation, the number of variables that is both active and not yet 
    /// marked for elimination in the `PresolveIndex` can only be derived using that data structure.
    fn nr_active_constraints(&self) -> usize {
        self.constraints.nr_rows()
    }

    /// The number of constraints in this linear program.
    ///
    /// # Return value
    ///
    /// The number of constraints, which excludes any variable bounds.
    pub fn nr_constraints(&self) -> usize {
        self.nr_active_constraints()
    }
    
    /// Number of variables that have not been eliminated after a presolving operation.
    /// 
    /// During a presolving operation, the number of variables that is both active and not yet 
    /// marked for elimination in the `PresolveIndex` can only be derived using that data structure.
    fn nr_active_variables(&self) -> usize {
        self.constraints.nr_columns()
    }
    
    /// Number of variables at the time this object was first created.
    fn nr_original_variables(&self) -> usize {
        self.original_variables.len()
    }

    /// The number of variables in this linear program.
    ///
    /// # Return value
    ///
    /// The number of columns / variables, which includes the slack columns / variables.
    pub fn nr_variables(&self) -> usize {
        self.nr_active_variables()
    }
}
/// A variable as part of a linear problem without restrictions (as opposed to for a `MatrixData` variable).
///
/// A variable is named, of continuous or integer type and may be shifted and flipped w.r.t. how it
/// was originally present in the problem.
///
/// The upper bound is relative to the shift; that is, the lower bound is `lower_bound - shift`, the
/// upper bound is `upper_bound - shift`. For example, the range stays the same, regardless of the
/// shift.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Variable<F> {
    /// Whether the variable is integer or not.
    pub variable_type: VariableType,
    /// Coefficient in the objective function.
    pub cost: F,
    /// Describing the accepted values for this variable
    ///
    /// Lower bound should be set to 0 when a variable is nonnegative. If it is `None`, the variable
    /// is considered to be in (-oo, upper_bound).
    pub lower_bound: Option<F>,
    /// Describing the accepted values for this variable
    ///
    /// If it is `None`, the variable is considered to be in (lower_bound, oo).
    pub upper_bound: Option<F>,
    /// How much this variable was shifted to have a zero lower bound.
    ///
    /// To find the "true" solution value, one needs to subtract this shift from the solution value
    /// produced by an optimization routine using the lower bound of 0.
    pub shift: F,
    /// Whether this variable was originally negative.
    ///
    /// To find the "true" solution value, one needs to multiply the solutionvalue found by -1, and
    /// then shift the value by the `shifted_by` field value.
    pub flipped: bool,
}

impl<OF> Variable<OF>
where
    OF: OrderedField,
    for<'r> &'r OF: OrderedFieldRef<OF>,
{
    /// Whether the variable allows only a single value.
    ///
    /// # Return value
    ///
    /// `Some` with the value if so, `None` otherwise.
    fn is_fixed(&self) -> Option<&OF> {
        match (&self.lower_bound, &self.upper_bound) {
            (Some(lower), Some(upper)) if lower == upper => Some(lower),
            _ => None,
        }
    }

    /// Whether a variable is unconstrained (has no bounds).
    fn is_free(&self) -> bool {
        self.lower_bound.is_none() && self.upper_bound.is_none()
    }

    /// Whether the variable admits a feasible value (and the upper bound is not below the lower
    /// bound).
    fn has_feasible_value(&self) -> bool {
        match (&self.lower_bound, &self.upper_bound) {
            (Some(lower), Some(upper)) => lower <= upper,
            _ => true,
        }
    }

    /// Get any feasible value, if there is one.
    fn get_feasible_value(&self) -> Option<OF> {
        if self.has_feasible_value() {
            self.upper_bound.as_ref().or_else(|| self.lower_bound.as_ref()).cloned().or_else(|| Some(OF::zero()))
        } else { None }
    }

    /// Change the lower bound if the given bound is higher.
    ///
    /// # Arguments
    ///
    /// * `new`: New value to compare the existing bound against (if there is any)
    ///
    /// # Return value
    ///
    /// `None` if the bound was not updated, `Some` if it was. If there was no bound before, a
    /// `None`, otherwise, a `Some` with the difference between the old and new bound. This
    /// difference is always strictly positive, as the lower bound can only be increased.
    fn update_lower_bound(&mut self, new: &OF) -> Option<Option<OF>> {
        Self::update_bound(&mut self.lower_bound, new, |new, existing| new > existing)
    }

    /// Change the upper bound if the given bound is lower.
    ///
    /// # Arguments
    ///
    /// * `new`: New value to compare the existing bound against (if there is any)
    ///
    /// # Return value
    ///
    /// `None` if the bound was not updated, `Some` if it was. If there was no bound before, a
    /// `None`, otherwise, a `Some` with the difference between the old and new bound. This
    /// difference is always strictly negative, as the upper bound can only be decreased.
    fn update_upper_bound(&mut self, new: &OF) -> Option<Option<OF>> {
        Self::update_bound(&mut self.upper_bound, new, |new, existing| new < existing)
    }

    /// Update either the lower or upper bound, if it makes the bound tighter.
    ///
    /// This is a helper method.
    ///
    /// # Arguments
    ///
    /// * `existing_bound`: Bound that could be updated.
    /// * `new`: New bound value.
    /// * `is_better`: A predicate indicating whether the new bound is better than the old bound, or
    /// not.
    ///
    /// # Return value
    ///
    /// `None` if the bound was not updated, `Some` if it was. If there was no bound before, a
    /// `None`, otherwise, a `Some` with the difference between the old and new bound.
    fn update_bound<P: Fn(&OF, &OF) -> bool>(
        existing_bound: &mut Option<OF>,
        new: &OF,
        is_better: P,
    ) -> Option<Option<OF>> {
        match existing_bound {
            Some(existing) => {
                if is_better(new, existing) {
                    let difference = new - &*existing;
                    *existing = new.clone();
                    Some(Some(difference))
                } else {
                    None
                }
            }
            None => {
                *existing_bound = Some(new.clone());
                Some(None)
            }
        }
    }
}

#[cfg(test)]
mod test {
    use num::FromPrimitive;

    use crate::data::linear_algebra::matrix::{ColumnMajor, Order};
    use crate::data::linear_algebra::vector::DenseVector;
    use crate::data::linear_algebra::vector::test::TestVector;
    use crate::data::linear_program::elements::{Objective, RangedConstraintRelation, VariableType};
    use crate::data::linear_program::general_form::{GeneralForm, OriginalVariable, Variable};
    use crate::data::number_types::rational::Rational32;
    use crate::R32;

    /// Shifting a variable.
    #[test]
    fn shift_variables() {
        type T = Rational32;

        let bound_value = 2.5_f64;

        let data = vec![
            vec![1, 0],
            vec![2, 1],
        ];
        let constraints = ColumnMajor::from_test_data::<T, T, _>(&data, 2);
        let b = DenseVector::from_test_data(vec![
            2,
            8,
        ]);
        let constraint_types = vec![
            RangedConstraintRelation::Greater,
            RangedConstraintRelation::Less,
        ];
        let variables = vec![Variable {
            variable_type: VariableType::Continuous,
            cost: R32!(1),
            lower_bound: None,
            upper_bound: None,
            shift: R32!(0),
            flipped: false
        }, Variable {
            variable_type: VariableType::Continuous,
            cost: R32!(3),
            lower_bound: Some(R32!(bound_value)),
            upper_bound: None,
            shift: R32!(0),
            flipped: false
        }];
        let variable_names = vec!["XONE".to_string(), "XTWO".to_string()];
        let mut general_form = GeneralForm::new(
            Objective::Minimize,
            constraints,
            constraint_types,
            b,
            variables,
            variable_names,
            R32!(1),
        );
        general_form.transform_variables();

        let expected = GeneralForm {
            objective: Objective::Minimize,
            fixed_cost: R32!(1) + R32!(3) * R32!(bound_value),
            constraints: ColumnMajor::from_test_data(&[
                vec![1, 0, -1],
                vec![2, 1, -2],
            ], 3),
            constraint_types: vec![
                RangedConstraintRelation::Greater,
                RangedConstraintRelation::Less,
            ],
            b: DenseVector::from_test_data(vec![
                2f64 - bound_value * 0f64,
                8f64 - bound_value * 1f64,
            ]),
            variables: vec![
                Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(1),
                    lower_bound: Some(R32!(0)),
                    upper_bound: None,
                    shift: R32!(0),
                    flipped: false,
                },
                Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(3),
                    lower_bound: Some(R32!(0)),
                    upper_bound: None,
                    shift: -R32!(bound_value),
                    flipped: false,
                },
                Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(-1),
                    lower_bound: Some(R32!(0)),
                    upper_bound: None,
                    shift: R32!(0),
                    flipped: false,
                },
            ],
            original_variables: vec![
                ("XONE".to_string(), OriginalVariable::ActiveFree(0, 2)),
                ("XTWO".to_string(), OriginalVariable::Active(1)),
            ],
            from_active_to_original: vec![0, 1, 0],
        };

        assert_eq!(general_form, expected);
    }

    #[test]
    fn make_b_non_negative() {
        type T = Rational32;

        let rows = ColumnMajor::from_test_data::<T, T, _>(&[vec![2]], 1);
        let b = DenseVector::from_test_data(vec![-1]);
        let variables = vec![
            Variable {
                variable_type: VariableType::Continuous,
                cost: R32!(1),
                lower_bound: None,
                upper_bound: None,
                shift: R32!(0),
                flipped: false
            },
        ];
        let variable_names = vec!["XONE".to_string()];
        let constraints = vec![RangedConstraintRelation::Equal];
        let mut result = GeneralForm::new(
            Objective::Minimize,
            rows,
            constraints,
            b,
            variables,
            variable_names,
            R32!(0),
        );
        result.make_b_non_negative();

        let data = ColumnMajor::from_test_data(&[vec![-2]], 1);
        let b = DenseVector::from_test_data(vec![1]);
        let constraints = vec![RangedConstraintRelation::Equal];
        let variables = vec![Variable {
            variable_type: VariableType::Continuous,
            cost: R32!(1),
            lower_bound: None,
            upper_bound: None,
            shift: R32!(0),
            flipped: false
        }];
        let variable_names = vec!["XONE".to_string()];
        let expected = GeneralForm::new(
            Objective::Minimize,
            data,
            constraints,
            b,
            variables,
            variable_names,
            R32!(0),
        );

        assert_eq!(result, expected);
    }
}
