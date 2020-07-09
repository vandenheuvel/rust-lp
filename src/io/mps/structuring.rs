//! # Organizing data read
//!
//! Organizing data read in the `parsing` module, and checking the linear program for consistency.
//! Contains also the definition of the struct representing linear programs in MPS format.
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::convert::TryFrom;
use std::convert::TryInto;
use std::mem;

use crate::data::linear_algebra::matrix::{ColumnMajor, Order as MatrixOrder, Sparse};
use crate::data::linear_algebra::SparseTuple;
use crate::data::linear_algebra::vector::{Dense as DenseVector, Vector};
use crate::data::linear_program::elements::{ConstraintType, VariableType};
use crate::data::linear_program::elements::Objective;
use crate::data::linear_program::general_form::GeneralForm;
use crate::data::linear_program::general_form::Variable as ShiftedVariable;
use crate::data::number_types::traits::{Field, OrderedField};
use crate::io::error::Inconsistency;
use crate::io::mps::{Bound, Range};
use crate::io::mps::BoundType;
use crate::io::mps::Constraint;
use crate::io::mps::parsing::{UnstructuredBound, UnstructuredRange};
use crate::io::mps::parsing::UnstructuredColumn;
use crate::io::mps::parsing::UnstructuredMPS;
use crate::io::mps::parsing::UnstructuredRhs;
use crate::io::mps::parsing::UnstructuredRow;
use crate::io::mps::Rhs;
use crate::io::mps::Variable;

impl<'a, 'b, F> TryFrom<UnstructuredMPS<'a, F>> for MPS<F> {
    type Error = Inconsistency;

    /// Try to convert into an `MPS` instance.
    ///
    /// This method organizes and structures the unstructured information contained in the
    /// `UnstructuredMPS` instance.
    ///
    /// The instance owns all the data it contains, implying that data is copied here from the
    /// original program string to owned `String`s.
    ///
    /// # Arguments
    ///
    /// * `unstructured_mps`: A possibly inconsistent MPS instance.
    ///
    /// # Return value
    ///
    /// A "structured" MPS.
    ///
    /// # Errors
    ///
    /// An inconsistency error if the program is inconsistent. An example of an inconsistency is a
    /// bound for a column which is unknown.
    fn try_from(unstructured_mps: UnstructuredMPS<'a, F>) -> Result<Self, Self::Error> {
        // Indices are used, because we can not trust the order of the input data.

        let (row_names, row_index) = build_row_index(&unstructured_mps.rows)?;
        let rows = order_rows(unstructured_mps.rows, &row_index);
        let (cost_values, columns, column_names) = build_columns(
            unstructured_mps.columns,
            &unstructured_mps.cost_row_name,
            &row_index,
        )?;
        let rhss = build_rhss(unstructured_mps.rhss, &row_index)?;
        let ranges = build_ranges(unstructured_mps.ranges, &row_index)?;
        let column_index = build_column_index(&column_names);
        let bounds = build_bounds(unstructured_mps.bounds, &column_index)?;

        // Owning all data
        let name = unstructured_mps.name.to_string();
        let cost_row_name = unstructured_mps.cost_row_name.to_string();

        Ok(MPS::new(
            name,
            cost_row_name,
            cost_values,
            row_names,
            rows,
            column_names,
            columns,
            rhss,
            ranges,
            bounds,
        ))
    }
}

/// Extract all row names, and assign to each row a fixed row index.
///
/// This index will be used throughout building the `MPS`, and ordered as in the original problem.
/// It is needed because constraint values for columns might appear in an unsorted order, so then a
/// lookup is practical.
///
/// # Arguments
///
/// * `unstructured_rows`: Collection of unstructured rows, cost row should not be in here.
///
/// # Return value
///
/// Names of all rows (this excludes the cost row) and a map from the name to the index at which
/// they are stored. This arbitrary order is the order in which rows will be stored in the final
/// MPS data structure.
///
/// # Errors
///
/// If there are rows with duplicate names, the MPS is considered inconsistent.
fn build_row_index<'a>(
    unstructured_rows: &Vec<UnstructuredRow<'a>>
) -> Result<(Vec<String>, HashMap<&'a str, usize>), Inconsistency> {
    // Check whether all defined rows have a unique name
    let mut row_index = HashMap::with_capacity(unstructured_rows.len());
    let unique = unstructured_rows.iter().enumerate()
        // Create a map from the row names to their index
        .all(|(i, row)| row_index.insert(row.name, i).is_none());
    if !unique {
        return Err(Inconsistency::new("Not all row names are unique"));
    }

    // Copy the names from the original file
    let row_names = unstructured_rows.iter()
        .map(|&UnstructuredRow { name, .. }| name.to_string())
        .collect();

    Ok((row_names, row_index))
}

/// Order the rows according the provided index.
///
/// The name of the row is match against the index.
///
/// # Arguments
///
/// * `unstructured_rows`: Collection of unstructured rows sorted according to the index (excludes
/// cost row).
/// * `index`: Corresponds to the ordering of `unstructured_rows`. Only used for debug.
///
/// # Return value
///
/// `Vec` of `Constraint`s ordered like the input.
#[cfg_attr(not(debug_assertions), allow(dead_code))] // The index is not used outside the debug assert
fn order_rows<'a>(
    unstructured_rows: Vec<UnstructuredRow<'a>>,
    index: &HashMap<&'a str, usize>,
) -> Vec<Constraint> {
    // Assume that the unstructured rows are sorted correctly: rows are at their correct index
    debug_assert!(unstructured_rows.iter().enumerate().all(|(i, ur)| {
        index.get(ur.name).map_or(false, |&ii| ii == i)
    }));

    // Assign the index based on order, although it should be identical to the order given by
    // index.
    unstructured_rows.iter().enumerate()
        .map(|(i, row)| {
            Constraint {
                name_index: i,
                constraint_type: row.constraint_type,
            }
        })
        .collect()
}

/// Collect the column values.
///
/// # Arguments
///
/// * `unstructured_columns`: Collection of `UnstructuredColumn`s.
/// * `cost_row_name`: Name of the cost row, used to identify which values belong to the objective
/// function.
/// * `row_index`: Index providing an ordering of the row names (excluding cost row name).
///
/// # Return value
///
/// If successful, the column ind
///
/// # Note
///
/// Assumes that the unstructured columns are sorted together. In the MPS file, all data for a
/// column should be provided together.
///
/// TODO: Generalize the method to relax the above assumption.
fn build_columns<'a, 'b, F>(
    mut unstructured_columns: Vec<UnstructuredColumn<F>>,
    cost_row_name: &'a str,
    row_index: &HashMap<&'a str, usize>,
) -> Result<(Vec<SparseTuple<F>>, Vec<Variable<F>>, Vec<String>), Inconsistency> {
    // Sort all values by column such that a change of column name indicates that the final value of
    // a column has been processed.
    //
    // This sorting is a bit expensive, but it simplifies the code that comes after.
    //
    // TODO: The ordering for equal column names doesn't need to stay the same. Could unstable
    //  sorting be beneficial? These keys need to be sorted again after, but then again in an
    //  unstable manner (rows keys should be unique).
    //  Some of this also applies to the methods below.
    unstructured_columns.sort_by_key(|&UnstructuredColumn { name, .. }| name);

    // Collector for cost values, (column_index, coefficient) tuples, is returned from function.
    // Others are assumed zero. The column index corresponds to the location in the `columns` Vec.
    let mut cost_values = Vec::new();
    // Collector for the columns, returned from function.
    let mut columns = Vec::new();
    // Collector for cloned column names, returned from function.
    let mut column_names = Vec::new();

    // Implementation uses first element to initialize values that test for change, so need at least
    // one value to be able to index with 0 on the unstructured_column vec.
    if !unstructured_columns.is_empty() {
        let UnstructuredColumn { name: mut current_name, variable_type: mut current_type, .. } = unstructured_columns[0];
        // Values per column, again as sparse (row_index, coefficient) tuples.
        let mut values = Vec::new();
        for UnstructuredColumn { name, variable_type, row_name, value, } in unstructured_columns {
            // Because unstructured_columns is sorted, a name change implies all values of a column
            // can be processed after each other. A name change implies that the previous column is
            // done.
            if name != current_name {
                // Sorted by row id
                // TODO: Unstable sorting? See sorting of unstructured_columns.
                values.sort_by_key(|&(i, _)| i);
                // Check for multiple row entries (even if two entries have the same value, it is
                // still considered an inconsistency).
                let unduplicated_length = values.len();
                values.dedup_by_key(|&mut (i, _)| i);
                if values.len() != unduplicated_length {
                    return Err(Inconsistency::new(
                        format!("Duplicate row value for column \"{}\".", current_name),
                    ));
                }

                // Save the value at the latest index (effectively alphabetically).
                columns.push(Variable {
                    name_index: column_names.len(),
                    variable_type: current_type,
                    values,
                });
                column_names.push(current_name.to_string());

                // Update the values for the new column that will now be processed.
                values = Vec::new();
                current_name = name;
                current_type = variable_type;
            }

            // Save coefficients in the cost functions separately.
            if row_name == cost_row_name {
                cost_values.push((column_names.len(), value));
            } else {
                let index = match row_index.get(row_name) {
                    Some(&index) => index,
                    None => return Err(Inconsistency::new(format!("Row name \"{}\" not known.", row_name))),
                };
                values.push((index, value));
            }
        }
        // Save the values from the last column.
        values.sort_by_key(|&(i, _)| i);
        columns.push(Variable { name_index: column_names.len(), variable_type: current_type, values, });
        column_names.push(current_name.to_string());
    }

    // Test for multiple cost values
    let unduplicated_length = cost_values.len();
    cost_values.dedup_by_key(|&mut (i, _)| i);
    if cost_values.len() != unduplicated_length {
        return Err(Inconsistency::new("Duplicate cost row value for a column."));
    }

    Ok((cost_values, columns, column_names))
}

/// Structure the right-hand side data.
///
/// Multiple right hand sides might be specified.
///
/// # Arguments
///
/// * `unstructured_rhss`: Collection of unstructured right-hand side data.
/// * `row_index`: Assignment of rows (by name) to an index.
///
/// # Return value
///
/// Vector of right-hands side constraints .
///
/// # Errors
///
/// If a constraint name has not already been encountered in the ROWS section.
///
/// TODO: Generalize the method as to relax the above requirement.
fn build_rhss<'a, 'b, F>(
    mut unstructured_rhss: Vec<UnstructuredRhs<'a, F>>,
    row_index: &HashMap<&'a str, usize>,
) -> Result<Vec<Rhs<F>>, Inconsistency> {
    // See the build_columns method for a description of the below logic.
    unstructured_rhss.sort_by_key(|&UnstructuredRhs { name, .. }| name);

    let mut rhss = Vec::new();
    if !unstructured_rhss.is_empty() {
        let UnstructuredRhs { name: mut current_name, .. } = unstructured_rhss[0];
        let mut values = Vec::new();
        for UnstructuredRhs { name, row_name, value, } in unstructured_rhss {
            if name != current_name {
                rhss.push(Rhs { name: current_name.to_string(), values, });

                current_name = name;
                values = Vec::new()
            }
            values.push((match row_index.get(row_name) {
                Some(&index) => index,
                None => return Err(Inconsistency::new(format!("Row name \"{}\" not known.", row_name))),
            }, value));
        }
        values.sort_by_key(|&(i, _)| i);
        rhss.push(Rhs { name: current_name.to_string(), values, });
    }

    Ok(rhss)
}

/// Structure the ranges data.
///
/// # Arguments
///
/// * `unstructured_ranges`: Collection of unstructured ranges.
/// * `row_index`: Assignment of rows (by name) to an index.
///
/// # Return value
///
/// Vector of ranges.
///
/// # Errors
///
/// If a constraint name has not already been encountered in the ROWS section, a
/// `LinearProgramError`.
///
/// # Note
///
/// Requires that right-hand side data is sorted by the name of the right-hand side.
///
/// TODO: Generalize the method as to relax the above requirement.
fn build_ranges<'a, F>(
    mut unstructured_ranges: Vec<UnstructuredRange<'a, F>>,
    row_index: &HashMap<&'a str, usize>,
) -> Result<Vec<Range<F>>, Inconsistency> {
    // See the build_columns method for a description of the below logic.
    let names_ok = {
        let mut unique = HashSet::with_capacity(unstructured_ranges.len());
        unstructured_ranges.iter()
            .all(|range| {
                row_index.get(range.row_name).map_or(false, |index| unique.insert(index))
            })
    };
    if !names_ok {
        return Err(Inconsistency::new("A row name is unknown or a row has multiple ranges."));
    }

    unstructured_ranges.sort_by_key(|&UnstructuredRange { name, .. }| name);

    let mut ranges = Vec::new();
    if !unstructured_ranges.is_empty() {
        let UnstructuredRange { name: mut current_name, .. } = unstructured_ranges[0];
        let mut values = Vec::new();
        for UnstructuredRange { name, row_name, value, } in unstructured_ranges {
            if name != current_name {
                ranges.push(Range { name: current_name.to_string(), values, });

                current_name = name;
                values = Vec::new()
            }
            values.push((match row_index.get(row_name) {
                Some(&index) => index,
                None => return Err(Inconsistency::new(format!("Row name \"{}\" not known.", row_name))),
            }, value));
        }
        values.sort_by_key(|&(i, _)| i);
        ranges.push(Range { name: current_name.to_string(), values, });
    }

    Ok(ranges)
}

/// Assign to each column a fixed column index.
///
/// This index will be used to organize the bounds information.
///
/// # Arguments
///
/// * `column_names`: Collection of names of the columns
///
/// # Return value
///
/// A map assigning to each column name a value.
///
/// # Note
///
/// This assignment is not a specific order.
fn build_column_index(column_names: &Vec<String>) -> HashMap<&str, usize> {
    column_names.iter()
        .enumerate()
        .map(|(index, name)| (name.as_str(), index))
        .collect()
}

/// Structure the bound data
///
/// # Arguments
///
/// * `unstructured_bounds`: Collection of unstructured bound data.
/// * `column_index`: Assignment of columns/variables (by name) to an index.
///
/// # Return value
///
/// Vector of bound constraints if successful, `LinearProgramError` if not.
///
/// # Note
///
/// Requires that bound data is sorted.
///
/// TODO: Generalize the method as to relax the above requirement.
fn build_bounds<F>(
    mut unstructured_bounds: Vec<UnstructuredBound<F>>,
    column_index: &HashMap<&str, usize>,
) -> Result<Vec<Bound<F>>, Inconsistency> {
    // See the build_columns method for a description of the below logic.
    unstructured_bounds.sort_by_key(|&UnstructuredBound { name, .. }| name);

    let mut bounds = Vec::new();
    if !unstructured_bounds.is_empty() {
        let UnstructuredBound { name: mut bound_name, .. } = unstructured_bounds[0];
        let mut values = Vec::new();
        // TODO: Get the innver bound type value fro m &f to f
        for UnstructuredBound { name, bound_type, column_name, } in unstructured_bounds {
            if name != bound_name {
                bounds.push(Bound { name: bound_name.to_string(), values, });

                bound_name = name;
                values = Vec::new();
            }
            values.push((
                bound_type,
                match column_index.get(column_name) {
                    Some(&index) => index,
                    None => return Err(Inconsistency::new(format!("Variable \"{}\" not known", column_name))),
                },
            ));
        }
        bounds.push(Bound { name: bound_name.to_string(), values, });
    }

    Ok(bounds)
}

/// Represents the contents of a MPS file in a structured manner.
///
/// `usize` variables in contained structs refer to the index of the cost and row names.
#[derive(Debug, PartialEq)]
pub struct MPS<F> {
    /// Name of the linear program.
    name: String,
    /// Name of the cost row.
    cost_row_name: String,
    /// Variable index and value tuples, describing how the variables appear in the objective
    /// function.
    ///
    /// The column indices in these tuples correspond to the names in column_names.
    cost_values: Vec<SparseTuple<F>>,
    /// Name of every constraint row.
    row_names: Vec<String>,
    /// All named constraint types (see the ConstraintType enum).
    ///
    /// Ordering corresponds to the row_names field.
    rows: Vec<Constraint>,
    /// Name of every variable.
    column_names: Vec<String>,
    /// Constraint name and variable name combinations.
    ///
    /// Ordering in each variable corresponds to the row_names field.
    columns: Vec<Variable<F>>,
    /// Right-hand side constraint values.
    ///
    /// Ordering in each right hand side corresponds to the row_names field.
    rhss: Vec<Rhs<F>>,
    /// Limiting constraint activations two-sidedly.
    ranges: Vec<Range<F>>,
    /// Extra bounds on variables.
    bounds: Vec<Bound<F>>,
}

impl<F> MPS<F> {
    /// Collect structured information into a `MPS` instance.
    ///
    /// # Arguments
    ///
    /// * `name`: Name of the linear program.
    /// * `cost_row_name`: Name of the cost row / objective function.
    /// * `cost_values`: Column (by index) and coefficient combinations for the objective function.
    /// * `row_names`: Names of all rows. The ordering corresponds with the order of the data in
    /// the `rows` argument.
    /// * `rows`: Constraint types. Names of the constraints are in the `row_names` argument, with
    /// corresponding order.
    /// * `column_names`: Names of all columns / variables. The ordering corresponds with the order
    /// of the data in the `columns` argument.
    /// * `columns`: Constraint data by column. Names of the variables are in `column_names`, with
    /// corresponding order.
    /// * `rhss`: Constraint values.
    /// * `ranges`: Flexibility on constraint activations, see the struct documentation.
    /// * `bounds`: Separate set of constraints, applying to the variables.
    ///
    /// # Return value
    ///
    /// All argument data collected in the `MPS` type.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn new(
        name: String,
        cost_row_name: String,
        cost_values: Vec<SparseTuple<F>>,
        row_names: Vec<String>,
        rows: Vec<Constraint>,
        column_names: Vec<String>,
        columns: Vec<Variable<F>>,
        rhss: Vec<Rhs<F>>,
        ranges: Vec<Range<F>>,
        bounds: Vec<Bound<F>>,
    ) -> Self {
        Self {
            name,
            cost_row_name,
            cost_values,
            row_names,
            rows,
            column_names,
            columns,
            rhss,
            ranges,
            bounds,
        }
    }
}

impl<OF: OrderedField, NZ> TryInto<GeneralForm<OF, NZ>> for MPS<OF>
where
{
    type Error = Inconsistency;

    /// Convert an `MPS` into a `GeneralForm` linear program.
    ///
    /// # Arguments
    ///
    /// * `self`: `MPS` instance.
    ///
    /// # Return value
    ///
    /// A linear program in general form.
    /// 
    /// # Errors
    /// 
    /// TODO: When can errors occur?
    fn try_into(self) -> Result<GeneralForm<OF, NZ>, Self::Error> {
        let (variable_info, variable_names) = compute_variable_info(
            &self.columns,
            self.column_names,
            self.cost_values,
            self.bounds,
        )?;
        let (columns, constraint_types, b) = compute_constraint_info(
            self.rows,
            self.rhss,
            self.columns,
            self.ranges,
        )?;

        Ok(GeneralForm::new(
            Objective::Minimize,
            columns,
            constraint_types,
            b,
            variable_info,
            variable_names,
            OF::zero(),
        ))
    }
}

/// Convert the variable oriented information of a MPS into `GeneralForm` fields.
///
/// # Arguments
///
/// * `columns`: MPS columns used to get the name index and cost value from.
/// * `column_names`: Names that get moved into the general form.
/// * `cost_values`: Coefficients in cost function.
/// * `bounds`: Variable bounds.
///
/// # Return value
///
/// Collection of variable info's as used in the GeneralForm and a collection of the column names.
///
/// # Errors
///
/// If there is an inconsistency in bound information, such as a trivial infeasibility.
fn compute_variable_info<OF: OrderedField, NZ>(
    columns: &Vec<Variable<OF>>,
    mut column_names: Vec<String>,
    cost_values: Vec<SparseTuple<OF>>,
    bounds: Vec<Bound<OF>>,
) -> Result<(Vec<ShiftedVariable<OF, NZ>>, Vec<String>), Inconsistency> {
    // Reorder the column names.
    {
        // This set should be a permutation.
        let nr_columns = columns.len();
        let as_set = columns.iter().map(|variable| variable.name_index).collect::<HashSet<_>>();
        debug_assert!(as_set.len() == nr_columns);
        debug_assert!(columns.iter().all(|variable| variable.name_index < nr_columns));
    }
    let variable_names = columns.iter().map(|variable| {
        mem::take(&mut column_names[variable.name_index])
    }).collect::<Vec<_>>();
    {
        // Each value should have been taken out exactly once, and each location should have been
        // filled exactly once.
        let value_substituted = String::default();
        debug_assert!(column_names.iter().all(|name| name == &value_substituted));
        drop(column_names);
        debug_assert!(variable_names.iter().all(|name| name != &value_substituted));
    }

    // Read in the sparse cost values
    debug_assert!(cost_values.is_sorted_by_key(|&(j, _)| j));
    let mut cost_values = cost_values.into_iter().peekable();
    let mut variable_info = columns.iter().enumerate().map(|(j, variable)| {
        let cost = match cost_values.peek() {
            None => OF::zero(),
            Some(&(jj, _)) => if jj == j {
                cost_values.next().unwrap().1
            } else {
                OF::zero()
            },
        };
        ShiftedVariable {
            variable_type: variable.variable_type,
            cost,
            upper_bound: None,
            lower_bound: None,
            shift: OF::zero(),
            flipped: false
        }
    }).collect();
    process_bounds(&mut variable_info, bounds)?;

    Ok((variable_info, variable_names))
}

/// Modify the variable info bounds to contain the bound information.
///
/// # Arguments
///
/// * `variable_info`: Variables to add the bounds to.
/// * `bounds`: Bounds organized by name (this organization is not relevant and discarded).
///
/// # Errors
///
/// If there is a trivial infeasibility (a variable has no feasible values).
/// TODO: Consider changing this into an "this LP is infeasible" return type
fn process_bounds<OF: OrderedField, NZ>(
    variable_info: &mut Vec<ShiftedVariable<OF, NZ>>,
    bounds: Vec<Bound<OF>>,
) -> Result<(), Inconsistency> {
    // Variables should not have existing bounds in them, because a default bound will be substituted
    // in this function. Note that this debug statement doesn't entirely cover that, in theory, all
    // variables could be completely free (unlikely).
    debug_assert!(variable_info.iter().all(|variable| {
        variable.lower_bound.is_none() && variable.upper_bound.is_none()
    }));

    // Variables that "have been touched" will not get a default bound substituted.
    let mut needs_default_lower_bound = vec![true; variable_info.len()];
    let mut variable_is_free = vec![false; variable_info.len()];

    // Bound names are irrelevant, treat them all the same
    for bound in bounds {
        for (bound_type, variable_index) in bound.values {
            let variable = &mut variable_info[variable_index];
            let (needs_default_lower, is_free) = process_bound(bound_type, variable)?;
            variable_is_free[variable_index] |= is_free;
            needs_default_lower_bound[variable_index] &= needs_default_lower;
        }
    }

    let any_free_has_bound = variable_info.iter().enumerate().any(|(j, variable)| {
        variable_is_free[j] && {
            variable.lower_bound.is_some() || variable.upper_bound.is_some()
        }
    });
    if any_free_has_bound {
        return Err(Inconsistency::new("A variable is both free and bounded."));
    }

    fill_in_default_lower_bounds(variable_info, needs_default_lower_bound);

    Ok(())
}

/// Update the variable with this bound.
///
/// Variables start out completely unconstrained and repeated calling of this method constrains them
/// increasingly.
///
/// # Arguments
///
/// * `bound_type`: One of 8 or so bound types that the MPS describes.
/// * `variable`: Variable to which this bound will be applied.
///
/// # Return value
///
/// Whether the variable is free, and whether the bound still needs a default zero lower bound.
///
/// # Errors
///
/// Inconsistency error if this variable is no longer be feasible after adding the bound.
/// TODO: Consider changing this into an "this LP is infeasible" return type
fn process_bound<OF: OrderedField, NZ>(
    bound_type: BoundType<OF>,
    variable: &mut ShiftedVariable<OF, NZ>,
) -> Result<(bool, bool), Inconsistency> {
    match bound_type {
        BoundType::LowerContinuous(value) => {
            replace_existing_with(&mut variable.lower_bound, value, Ordering::Greater);
            Ok((false, false))
        },
        BoundType::UpperContinuous(value) => {
            // The implied zero lower bound gets filled in only if no other lower bound is present.
            // This behavior is copied from GLPK.
            replace_existing_with(&mut variable.upper_bound, value, Ordering::Less);
            Ok((true, false))
        },
        BoundType::Fixed(value) => {
            replace_existing_with(&mut variable.lower_bound, value.clone(), Ordering::Greater);
            replace_existing_with(&mut variable.upper_bound, value, Ordering::Less);
            Ok((false, false))
        }
        BoundType::Free => {
            // This check is not enough, because these bounds might be set later. The caller checks
            // as well after.
            if variable.lower_bound.is_some() || variable.upper_bound.is_some() {
                return Err(Inconsistency::new("Variable can't be bounded and free"))
            }
            Ok((false, true))
        },
        // Infinity bounds only say which sign a variable should have; the implied zero bound is
        // taken as an explicit one (it doesn't matter if another bound is present as well).
        BoundType::LowerMinusInfinity => {
            // Lower bound minus infinity is implied by variable.lower_bound.is_none()
            replace_existing_with(&mut variable.upper_bound, OF::zero(), Ordering::Less);
            Ok((false, false))
        },
        BoundType::UpperInfinity => {
            replace_existing_with(&mut variable.lower_bound, OF::zero(), Ordering::Greater);
            // Upper bound infinity is implied by variable.upper_bound.is_none()
            Ok((false, false))
        },
        BoundType::Binary => {
            replace_existing_with(&mut variable.lower_bound, OF::zero(), Ordering::Greater);
            replace_existing_with(&mut variable.upper_bound, OF::one(), Ordering::Less);
            variable.variable_type = VariableType::Integer;
            Ok((false, false))
        }
        BoundType::LowerInteger(value) => {
            replace_existing_with(&mut variable.lower_bound, value, Ordering::Greater);
            variable.variable_type = VariableType::Integer;
            Ok((false, false))
        },
        BoundType::UpperInteger(value) => {
            replace_existing_with(&mut variable.upper_bound, value, Ordering::Less);
            variable.variable_type = VariableType::Integer;
            Ok((true, false))
        },
        BoundType::SemiContinuous(_, _) => unimplemented!(),
    }
}

/// Tighten a bound.
///
/// # Arguments
///
/// * `option`: Bound value that will potentially be added or changed.
/// * `new_value`: Candidate value.
/// * `ordering`: Whether the new value should be smaller or larger than the existing bound value in
/// order for the bound to be changed (new w.r.t. old).
fn replace_existing_with<OF: OrderedField>(option: &mut Option<OF>, new_value: OF, ordering: Ordering) {
    // Nothing would change if they would need to be equal, so this doesn't make sense.
    debug_assert_ne!(ordering, Ordering::Equal);

    if let Some(ref mut existing_value) = option {
        if new_value.cmp(existing_value) == ordering {
            *existing_value = new_value;
        }
    } else {
        *option = Some(new_value);
    }
}

/// Fill in default lower bound for those variables that need it.
///
/// # Arguments
///
/// * `variables`: Variables with some bounds processed and default bounds not yet substituted.
/// * `needs_lower_bound`: Whether the variable at that index needs a lower bound.
fn fill_in_default_lower_bounds<OF: OrderedField, NZ>(
    variables: &mut Vec<ShiftedVariable<OF, NZ>>,
    needs_lower_bound: Vec<bool>,
) {
    debug_assert_eq!(variables.len(), needs_lower_bound.len());
    // Lower bounds should not have been touched yet for those variables that this method will modify.
    debug_assert!(variables.iter().zip(needs_lower_bound.iter()).all(|(variable, needs_lower)| {
        !needs_lower || variable.lower_bound.is_none()
    }));

    for (j, needs_lower_bound) in needs_lower_bound.into_iter().enumerate() {
        if needs_lower_bound {
            variables[j].lower_bound = Some(OF::zero());
        }
    }
}

/// Convert the constraint related information of a MPS into `GeneralForm` fields.
///
/// # Arguments
///
/// * `rhss`: Right-hand sides to be converted into a `b`.
/// * `columns`: Variables containing the constraint coefficients.
/// * `rows`: Direction of the constraint (name is not used).
/// * `range`: Flexibility for constraints.
fn compute_constraint_info<OF: OrderedField>(
    rows: Vec<Constraint>,
    rhss: Vec<Rhs<OF>>,
    columns: Vec<Variable<OF>>,
    ranges: Vec<Range<OF>>,
) -> Result<(Sparse<OF, OF, ColumnMajor>, Vec<ConstraintType>, DenseVector<OF>), Inconsistency> {
    let original_nr_rows = rows.len();

    // Flatten, we don't care about the different range names
    let mut range_rows = ranges.into_iter()
        .flat_map(|range| range.values.into_iter())
        .collect::<Vec<_>>();
    // We process them by row.
    // TODO: Order doesn't matter, use unstable sort?
    range_rows.sort_by_key(|&(i, _)| i);
    let unduplicated_length = range_rows.len();
    range_rows.dedup_by_key(|&mut (i, _)| i);
    if range_rows.len() < unduplicated_length {
        return Err(Inconsistency::new("Only one range per row can be specified."));
    }

    let columns = compute_columns(columns, original_nr_rows, &range_rows);
    let constraint_types = compute_constraint_types(&rows, &range_rows);
    let b = compute_b(rhss, &constraint_types, &rows, original_nr_rows, range_rows)?;

    Ok((columns, constraint_types, b))
}

/// Duplicate values within the columns when necessary according to the ranges.
///
/// # Arguments
///
/// * `columns`: MPS variables with sparse (row index, value) tuples.
/// * `original_nr_rows`: Number of rows when discarding the ranges.
/// * `ranges`: Tuples with (row index, r value) indicating where an extra range constraint should
/// be created.
///
/// # Return value
///
/// Column-major sparse matrix of constraint coefficients.
fn compute_columns<F: Field>(
    columns: Vec<Variable<F>>,
    original_nr_rows: usize,
    ranges: &Vec<SparseTuple<F>>,
) -> Sparse<F, F, ColumnMajor> {
    debug_assert!(ranges.is_sorted_by_key(|&(i, _)| i));
    debug_assert_eq!(ranges.iter().map(|&(i, _)| i).collect::<HashSet<_>>().len(), ranges.len());
    debug_assert!(columns.iter().all(|variable| {
        variable.values.is_sorted_by_key(|&(i, _)| i)
    }));
    debug_assert!(columns.iter().all(|variable| {
        variable.values.iter().all(|&(i, _)| i < original_nr_rows)
    }));
    debug_assert!(columns.iter().all(|variable| {
        variable.values.iter().map(|&(i, _)| i).collect::<HashSet<_>>().len() == variable.values.len()
    }));
    let nr_columns = columns.len();

    let mut new_columns = vec![Vec::new(); columns.len()];
    for (j, column) in columns.into_iter().enumerate() {
        let mut extra_done = 0;
        for (i, value) in column.values {
            while extra_done < ranges.len() && ranges[extra_done].0 < i {
                extra_done += 1;
            }
            new_columns[j].push((i + extra_done, value));

            if extra_done < ranges.len() && ranges[extra_done].0 == i {
                extra_done += 1;
                let value_copy = new_columns[j].last().unwrap().1.clone();
                new_columns[j].push((i + extra_done, value_copy));
            }
        }
    }

    ColumnMajor::new(
        new_columns,
        original_nr_rows + ranges.len(),
        nr_columns,
    )
}

/// Compute the constraint types by integrating bounds.
///
/// # Arguments
///
/// * `rows`: Contains the constraint types without ranges applied.
/// * `ranges`: Tuples with (row index, r value) indicating where an extra range constraint should
/// be created.
///
/// # Return value
///
/// Extended constraint types. See the documentation of the `UnstructuredRange` for more.
fn compute_constraint_types<F: Field>(
    rows: &Vec<Constraint>,
    ranges: &Vec<SparseTuple<F>>,
) -> Vec<ConstraintType> {
    debug_assert!(ranges.is_sorted_by_key(|&(i, _)| i));
    debug_assert!(ranges.iter().all(|&(i, _)| i < rows.len()));

    let mut constraint_types = Vec::with_capacity(rows.len() + ranges.len());
    let mut extra_done = 0;
    for (i, &constraint) in rows.iter().enumerate() {
        if extra_done < ranges.len() && ranges[extra_done].0 == i {
            while extra_done < ranges.len() && ranges[extra_done].0 == i {
                constraint_types.push(ConstraintType::Greater);
                constraint_types.push(ConstraintType::Less);
                extra_done += 1;
            }
        } else {
            constraint_types.push(constraint.constraint_type);
        }
    }

    constraint_types
}

/// Combine all constraint values.
///
/// # Arguments
///
/// * `rhss`: Right hand sides (often only one), b's values.
/// * `constraints`: Constraint directions (relevant for ranges).
/// * `rows`: Original constraint types.
/// * `original_nr_rows`: Number of constraints without ranges.
/// * `ranges`: Tuples with (row index, r value) indicating where an extra range constraint should
/// be created.
///
/// # Return value
///
/// A single right hand side.
///
/// # Errors
///
/// When there is a trivial infeasibility due to multiple equality bounds being specified with
/// different values.
#[allow(unreachable_patterns)]
fn compute_b<OF: OrderedField>(
    rhss: Vec<Rhs<OF>>,
    constraints: &Vec<ConstraintType>,
    rows: &Vec<Constraint>,
    original_nr_rows: usize,
    ranges: Vec<SparseTuple<OF>>,
) -> Result<DenseVector<OF>, Inconsistency> {
    let new_nr_rows = original_nr_rows + ranges.len();
    debug_assert!(rhss.iter().all(|rhs| rhs.values.is_sorted_by_key(|&(i, _)| i)));
    debug_assert!(rhss.iter().all(|rhs| rhs.values.iter().all(|&(i, _)| i < original_nr_rows)));

    // We fill be with options, and then explicitly substitute the default value later.
    let mut b = vec![None; new_nr_rows];
    for rhs in rhss {
        let mut extra_done = 0;
        for (i, value) in rhs.values {
            while extra_done < ranges.len() && ranges[extra_done].0 < i {
                extra_done += 1;
            }
            if extra_done < ranges.len() && ranges[extra_done].0 == i {
                let r = &ranges[extra_done].1;
                // See the documentation of `UnstructuredRhs` for the below logic.
                let r_abs = r.clone().abs();
                let (h, u) = match (rows[i].constraint_type, r.cmp(&OF::zero())) {
                    (ConstraintType::Greater, _) => (value.clone(), value + r_abs),
                    (ConstraintType::Less, _) => (value.clone() - r_abs, value),
                    (ConstraintType::Equal, Ordering::Greater | Ordering::Equal) => (value.clone(), value + r_abs),
                    (ConstraintType::Equal, Ordering::Less | Ordering::Equal) => (value.clone() - r_abs, value),
                };

                b[i + extra_done] = Some(h);
                extra_done += 1;
                b[i + extra_done] = Some(u);
            } else if let Some(current) = &mut b[i + extra_done] {
                match constraints[i + extra_done] {
                    ConstraintType::Equal => if &value != &*current {
                        return Err(Inconsistency::new(
                            format!("Trivial infeasibility: a constraint can't equal both {} and {}", current, value),
                        ))
                    },
                    ConstraintType::Greater => if &value > current {
                        *current = value;
                    },
                    ConstraintType::Less => if &value < current {
                        *current = value;
                    },
                }
            } else {
                b[i + extra_done] = Some(value);
            }
        }
    }

    // Substitute the default value.
    Ok(DenseVector::new(
        b.into_iter().map(|value| value.unwrap_or(OF::zero())).collect(),
        original_nr_rows + ranges.len(),
    ))
}

#[cfg(test)]
#[allow(clippy::shadow_unrelated)]
mod test {
    use num::rational::Ratio;

    use crate::data::linear_algebra::matrix::{ColumnMajor, Order};
    use crate::data::linear_algebra::vector::{Dense, Vector};
    use crate::data::linear_program::elements::{ConstraintType, VariableType};
    use crate::io::mps::{Constraint, Rhs, Variable};
    use crate::io::mps::structuring::{compute_b, compute_columns};
    use crate::num::FromPrimitive;
    use crate::R32;

    type T = Ratio<i32>;

    #[test]
    fn test_compute_columns() {
        // No ranges, no values
        let columns = vec![Variable {
            name_index: 0,
            variable_type: VariableType::Continuous,
            values: vec![],
        }];
        let original_nr_rows = 0;
        let ranges = vec![];
        let columns = compute_columns::<T, T>(columns, original_nr_rows, &ranges);
        assert_eq!(columns, ColumnMajor::new(vec![vec![]], 0, 1));

        // No ranges, some values
        let columns = vec![Variable {
            name_index: 0,
            variable_type: VariableType::Continuous,
            values: vec![(0, R32!(123))],
        }];
        let original_nr_rows = 2;
        let ranges = vec![];
        let columns = compute_columns::<T, T>(columns, original_nr_rows, &ranges);
        assert_eq!(columns, ColumnMajor::new(vec![vec![(0, R32!(123))]], 2, 1));
        let columns = vec![Variable {
            name_index: 0,
            variable_type: VariableType::Continuous,
            values: vec![(1, R32!(123))],
        }];
        let original_nr_rows = 2;
        let ranges = vec![];
        let columns = compute_columns::<T, T>(columns, original_nr_rows, &ranges);
        assert_eq!(columns, ColumnMajor::new(vec![vec![(1, R32!(123))]], 2, 1));

        // One range, no values
        let columns = vec![Variable {
            name_index: 0,
            variable_type: VariableType::Continuous,
            values: vec![],
        }];
        let original_nr_rows = 1;
        let ranges = vec![(0, R32!(1))];
        let columns = compute_columns::<T, T>(columns, original_nr_rows, &ranges);
        assert_eq!(columns, ColumnMajor::new(vec![vec![]], 2, 1));

        // One range, some values
        let columns = vec![Variable {
            name_index: 0,
            variable_type: VariableType::Continuous,
            values: vec![(0, R32!(1))],
        }];
        let original_nr_rows = 1;
        let ranges = vec![(0, R32!(1))];
        let columns = compute_columns::<T, T>(columns, original_nr_rows, &ranges);
        assert_eq!(columns, ColumnMajor::new(vec![vec![(0, R32!(1)), (1, R32!(1))]], 2, 1));

        // One range, value before range row
        let columns = vec![Variable {
            name_index: 0,
            variable_type: VariableType::Continuous,
            values: vec![(0, R32!(1))],
        }];
        let original_nr_rows = 2;
        let ranges = vec![(1, R32!(1))];
        let columns = compute_columns::<T, T>(columns, original_nr_rows, &ranges);
        assert_eq!(columns, ColumnMajor::new(vec![vec![(0, R32!(1))]], 3, 1));

        // One range, value after range row
        let columns = vec![Variable {
            name_index: 0,
            variable_type: VariableType::Continuous,
            values: vec![(1, R32!(1))],
        }];
        let original_nr_rows = 2;
        let ranges = vec![(0, R32!(1))];
        let columns = compute_columns::<T, T>(columns, original_nr_rows, &ranges);
        assert_eq!(columns, ColumnMajor::new(vec![vec![(2, R32!(1))]], 3, 1));
    }

    #[test]
    fn test_compute_b() {
        // No ranges, no data
        let rhss = vec![];
        let constraints = vec![];
        let rows = vec![];
        let original_nr_rows = 0;
        let ranges = vec![];
        let b = compute_b(rhss, &constraints, &rows, original_nr_rows, ranges);
        assert_eq!(b, Ok(Dense::<T>::new(vec![], 0)));

        // No ranges, one rhs
        let rhss = vec![Rhs { name: "R".to_string(), values: vec![(0, R32!(1))], }];
        let constraints = vec![ConstraintType::Equal];
        let rows = vec![Constraint { name_index: 0, constraint_type: ConstraintType::Equal}];
        let original_nr_rows = 1;
        let ranges = vec![];
        let b = compute_b(rhss, &constraints, &rows, original_nr_rows, ranges);
        assert_eq!(b, Ok(Dense::<T>::new(vec![R32!(1)], 1)));

        // No ranges, two rhses
        let rhss = vec![
            Rhs { name: "R1".to_string(), values: vec![(0, R32!(1))], },
            Rhs { name: "R2".to_string(), values: vec![(0, R32!(2))], },
        ];
        let constraints = vec![ConstraintType::Greater];
        let rows = vec![Constraint { name_index: 0, constraint_type: ConstraintType::Greater}];
        let original_nr_rows = 1;
        let ranges = vec![];
        let b = compute_b(rhss, &constraints, &rows, original_nr_rows, ranges);
        assert_eq!(b, Ok(Dense::<T>::new(vec![R32!(2)], 1)));

        // One range with data before
        let rhss = vec![
            Rhs { name: "R".to_string(), values: vec![(0, R32!(1)), (1, R32!(5))], },
        ];
        let constraints = vec![
            ConstraintType::Greater,
            ConstraintType::Equal,
        ];
        let rows = vec![
            Constraint { name_index: 0, constraint_type: ConstraintType::Greater,},
            Constraint { name_index: 1, constraint_type: ConstraintType::Equal,},
        ];
        let original_nr_rows = 2;
        let ranges = vec![(1, R32!(2))];
        let b = compute_b(rhss, &constraints, &rows, original_nr_rows, ranges);
        assert_eq!(b, Ok(Dense::<T>::new(vec![R32!(1), R32!(5), R32!(7)], 3)));

        // One range with data after
        let rhss = vec![
            Rhs { name: "R".to_string(), values: vec![(0, R32!(1)), (1, R32!(5))], },
        ];
        let constraints = vec![
            ConstraintType::Greater,
            ConstraintType::Equal,
        ];
        let rows = vec![
            Constraint { name_index: 0, constraint_type: ConstraintType::Greater,},
            Constraint { name_index: 1, constraint_type: ConstraintType::Equal,},
        ];
        let original_nr_rows = 2;
        let ranges = vec![(0, R32!(2))];
        let b = compute_b(rhss, &constraints, &rows, original_nr_rows, ranges);
        assert_eq!(b, Ok(Dense::<T>::new(vec![R32!(1), R32!(3), R32!(5)], 3)));
    }
}
