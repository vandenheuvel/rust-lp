//! # Organizing data read
//!
//! Organizing data read in the `parsing` module, and checking the linear program for consistency.
//! Contains also the definition of the struct representing linear programs in MPS format.
use std::cmp::Ordering;
use std::collections::HashSet;
use std::convert::TryInto;

use num_traits::{One, Zero};

use crate::data::linear_algebra::{SparseTuple, SparseTupleVec};
use crate::data::linear_algebra::matrix::{ColumnMajor, Order as MatrixOrder, Sparse};
use crate::data::linear_algebra::vector::{Dense as DenseVector, Vector};
use crate::data::linear_program::elements::{ConstraintType, VariableType};
use crate::data::linear_program::elements::Objective;
use crate::data::linear_program::general_form::GeneralForm;
use crate::data::linear_program::general_form::Variable as ShiftedVariable;
use crate::io::error::Inconsistency;
use crate::io::mps::{Bound, MPS, Range};
use crate::io::mps::BoundType;
use crate::io::mps::Column;
use crate::io::mps::number::traits::Field;
use crate::io::mps::number::traits::NonZero;
use crate::io::mps::Rhs;
use crate::io::mps::Row;

impl<F, NZ> TryInto<GeneralForm<F, NZ>> for MPS<F, NZ>
where
    F: Zero + One,
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
    fn try_into(self) -> Result<GeneralForm<F, NZ>, Self::Error> {
        let (variable_info, variable_values, variable_names) = compute_variable_info(
            self.columns,
            self.cost_values,
            self.bounds,
        )?;
        let (columns, constraint_types, b) = compute_constraint_info(
            self.rows,
            self.rhss,
            variable_values,
            self.ranges,
        )?;

        Ok(GeneralForm::new(
            Objective::Minimize,
            columns,
            constraint_types,
            b,
            variable_info,
            variable_names,
            F::zero(),
        ))
    }
}

/// Convert the variable oriented information of a MPS into `GeneralForm` fields.
///
/// # Arguments
///
/// * `columns`: MPS columns used to get the name index and cost value from.
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
fn compute_variable_info<F: Zero + One, NZ>(
    columns: Vec<Column<F>>,
    cost_values: Vec<SparseTuple<F>>,
    bounds: Vec<Bound<F>>,
) -> Result<(Vec<ShiftedVariable<F, NZ>>, Vec<SparseTupleVec<F>>, Vec<String>), Inconsistency> {
    // Reorder the column names.
    debug_assert_eq!(
        columns.iter().map(|variable| variable.name.as_str()).collect::<HashSet<_>>().len(),
        columns.len(),
    );

    // Read in the sparse cost values
    debug_assert!(cost_values.is_sorted_by_key(|&(j, _)| j));
    let mut cost_values = cost_values.into_iter().peekable();
    let (mut info, values_names): (Vec<_>, Vec<_>) = columns.into_iter()
        .enumerate()
        .map(|(j, variable)| {
            let cost = match cost_values.peek() {
                None => F::zero(),
                Some(&(jj, _)) => if jj == j {
                    cost_values.next().unwrap().1
                } else {
                    F::zero()
                },
            };

            (ShiftedVariable {
                variable_type: variable.variable_type,
                cost,
                upper_bound: None,
                lower_bound: None,
                shift: F::zero(),
                flipped: false,
            }, (variable.values, variable.name))
        }).unzip();
    process_bounds(&mut info, bounds)?;
    let (values, names): (Vec<_>, Vec<_>) = values_names.into_iter().unzip();

    Ok((info, values, names))
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
fn process_bounds<N: Zero + One, NZ>(
    variable_info: &mut Vec<ShiftedVariable<N, NZ>>,
    bounds: Vec<Bound<N>>,
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
        for (variable_index, bound_type) in bound.values {
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
fn process_bound<N: Zero + One, NZ>(
    bound_type: BoundType<N>,
    variable: &mut ShiftedVariable<N, NZ>,
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
            replace_existing_with(&mut variable.upper_bound, N::zero(), Ordering::Less);
            Ok((false, false))
        },
        BoundType::UpperInfinity => {
            replace_existing_with(&mut variable.lower_bound, N::zero(), Ordering::Greater);
            // Upper bound infinity is implied by variable.upper_bound.is_none()
            Ok((false, false))
        },
        BoundType::Binary => {
            replace_existing_with(&mut variable.lower_bound, N::zero(), Ordering::Greater);
            replace_existing_with(&mut variable.upper_bound, N::one(), Ordering::Less);
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
fn replace_existing_with<OF>(option: &mut Option<OF>, new_value: OF, ordering: Ordering) {
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
fn fill_in_default_lower_bounds<F: Zero, NZ>(
    variables: &mut Vec<ShiftedVariable<F, NZ>>,
    needs_lower_bound: Vec<bool>,
) {
    debug_assert_eq!(variables.len(), needs_lower_bound.len());
    // Lower bounds should not have been touched yet for those variables that this method will modify.
    debug_assert!(variables.iter().zip(needs_lower_bound.iter()).all(|(variable, needs_lower)| {
        !*needs_lower || variable.lower_bound.is_none()
    }));

    for (j, needs_lower_bound) in needs_lower_bound.into_iter().enumerate() {
        if needs_lower_bound {
            variables[j].lower_bound = Some(F::zero());
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
fn compute_constraint_info<N: Field, NZ: NonZero>(
    rows: Vec<Row>,
    rhss: Vec<Rhs<N>>,
    columns: Vec<SparseTupleVec<N>>,
    ranges: Vec<Range<N>>,
) -> Result<(Sparse<NZ, NZ, ColumnMajor>, Vec<ConstraintType>, DenseVector<N>), Inconsistency> {
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
fn compute_columns<N: Field, NZ: NonZero>(
    columns: Vec<SparseTupleVec<N>>,
    original_nr_rows: usize,
    ranges: &Vec<SparseTuple<N>>,
) -> Sparse<NZ, NZ, ColumnMajor> {
    debug_assert!(ranges.is_sorted_by_key(|&(i, _)| i));
    debug_assert_eq!(ranges.iter().map(|&(i, _)| i).collect::<HashSet<_>>().len(), ranges.len());
    debug_assert!(columns.iter().all(|variable| {
        variable.is_sorted_by_key(|&(i, _)| i)
    }));
    debug_assert!(columns.iter().all(|variable| {
        variable.iter().all(|&(i, _)| i < original_nr_rows)
    }));
    debug_assert!(columns.iter().all(|variable| {
        variable.iter().map(|&(i, _)| i).collect::<HashSet<_>>().len() == variable.len()
    }));
    let nr_columns = columns.len();

    let mut new_columns = vec![Vec::new(); columns.len()];
    for (j, column) in columns.into_iter().enumerate() {
        let mut extra_done = 0;
        for (i, value) in column {
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
    rows: &Vec<Row>,
    ranges: &Vec<SparseTuple<F>>,
) -> Vec<ConstraintType> {
    debug_assert!(ranges.is_sorted_by_key(|&(i, _)| i));
    // TODO: How about uniqueness?
    debug_assert!(ranges.iter().all(|&(i, _)| i < rows.len()));

    let mut constraint_types = Vec::with_capacity(rows.len() + ranges.len());
    let mut extra_done = 0;
    for (i, constraint) in rows.iter().enumerate() {
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
fn compute_b<OF: Field>(
    rhss: Vec<Rhs<OF>>,
    constraints: &Vec<ConstraintType>,
    rows: &Vec<Row>,
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

// #[cfg(test)]
// #[allow(clippy::shadow_unrelated)]
// mod test {
//     use num_traits::FromPrimitive;
//
//     use crate::data::linear_algebra::matrix::{ColumnMajor, Order};
//     use crate::data::linear_algebra::vector::{Dense, Vector};
//     use crate::data::linear_program::elements::ConstraintType;
//     use crate::io::mps::{Rhs, Row};
//     use crate::io::mps::convert::{compute_b, compute_columns};
//     use crate::R64;
//
//     type T = Ratio<i32>;
//
//     #[test]
//     fn test_compute_columns() {
//         // No ranges, no values
//         let columns = vec![vec![]];
//         let original_nr_rows = 0;
//         let ranges = vec![];
//         let columns = compute_columns::<T, T>(columns, original_nr_rows, &ranges);
//         assert_eq!(columns, ColumnMajor::new(vec![vec![]], 0, 1));
//
//         // No ranges, some values
//         let columns = vec![vec![(0, R64!(123))]];
//         let original_nr_rows = 2;
//         let ranges = vec![];
//         let columns = compute_columns::<T, T>(columns, original_nr_rows, &ranges);
//         assert_eq!(columns, ColumnMajor::new(vec![vec![(0, R64!(123))]], 2, 1));
//         let columns = vec![vec![(1, R64!(123))]];
//         let original_nr_rows = 2;
//         let ranges = vec![];
//         let columns = compute_columns::<T, T>(columns, original_nr_rows, &ranges);
//         assert_eq!(columns, ColumnMajor::new(vec![vec![(1, R64!(123))]], 2, 1));
//
//         // One range, no values
//         let columns = vec![vec![]];
//         let original_nr_rows = 1;
//         let ranges = vec![(0, R64!(1))];
//         let columns = compute_columns::<T, T>(columns, original_nr_rows, &ranges);
//         assert_eq!(columns, ColumnMajor::new(vec![vec![]], 2, 1));
//
//         // One range, some values
//         let columns = vec![vec![(0, R64!(1))]];
//         let original_nr_rows = 1;
//         let ranges = vec![(0, R64!(1))];
//         let columns = compute_columns::<T, T>(columns, original_nr_rows, &ranges);
//         assert_eq!(columns, ColumnMajor::new(vec![vec![(0, R64!(1)), (1, R64!(1))]], 2, 1));
//
//         // One range, value before range row
//         let columns = vec![vec![(0, R64!(1))]];
//         let original_nr_rows = 2;
//         let ranges = vec![(1, R64!(1))];
//         let columns = compute_columns::<T, T>(columns, original_nr_rows, &ranges);
//         assert_eq!(columns, ColumnMajor::new(vec![vec![(0, R64!(1))]], 3, 1));
//
//         // One range, value after range row
//         let columns = vec![vec![(1, R64!(1))]];
//         let original_nr_rows = 2;
//         let ranges = vec![(0, R64!(1))];
//         let columns = compute_columns::<T, T>(columns, original_nr_rows, &ranges);
//         assert_eq!(columns, ColumnMajor::new(vec![vec![(2, R64!(1))]], 3, 1));
//     }
//
//     #[test]
//     fn test_compute_b() {
//         // No ranges, no data
//         let rhss = vec![];
//         let constraints = vec![];
//         let rows = vec![];
//         let original_nr_rows = 0;
//         let ranges = vec![];
//         let b = compute_b(rhss, &constraints, &rows, original_nr_rows, ranges);
//         assert_eq!(b, Ok(Dense::<T>::new(vec![], 0)));
//
//         // No ranges, one rhs
//         let rhss = vec![Rhs { name: "R".to_string(), values: vec![(0, R64!(1))], }];
//         let constraints = vec![ConstraintType::Equal];
//         let rows = vec![Row { name: "".to_string(), constraint_type: ConstraintType::Equal}];
//         let original_nr_rows = 1;
//         let ranges = vec![];
//         let b = compute_b(rhss, &constraints, &rows, original_nr_rows, ranges);
//         assert_eq!(b, Ok(Dense::<T>::new(vec![R64!(1)], 1)));
//
//         // No ranges, two rhses
//         let rhss = vec![
//             Rhs { name: "R1".to_string(), values: vec![(0, R64!(1))], },
//             Rhs { name: "R2".to_string(), values: vec![(0, R64!(2))], },
//         ];
//         let constraints = vec![ConstraintType::Greater];
//         let rows = vec![Row { name: "".to_string(), constraint_type: ConstraintType::Greater}];
//         let original_nr_rows = 1;
//         let ranges = vec![];
//         let b = compute_b(rhss, &constraints, &rows, original_nr_rows, ranges);
//         assert_eq!(b, Ok(Dense::<T>::new(vec![R64!(2)], 1)));
//
//         // One range with data before
//         let rhss = vec![
//             Rhs { name: "R".to_string(), values: vec![(0, R64!(1)), (1, R64!(5))], },
//         ];
//         let constraints = vec![
//             ConstraintType::Greater,
//             ConstraintType::Equal,
//         ];
//         let rows = vec![
//             Row { name: "".to_string(), constraint_type: ConstraintType::Greater,},
//             Row { name: "".to_string(), constraint_type: ConstraintType::Equal,},
//         ];
//         let original_nr_rows = 2;
//         let ranges = vec![(1, R64!(2))];
//         let b = compute_b(rhss, &constraints, &rows, original_nr_rows, ranges);
//         assert_eq!(b, Ok(Dense::<T>::new(vec![R64!(1), R64!(5), R64!(7)], 3)));
//
//         // One range with data after
//         let rhss = vec![
//             Rhs { name: "R".to_string(), values: vec![(0, R64!(1)), (1, R64!(5))], },
//         ];
//         let constraints = vec![
//             ConstraintType::Greater,
//             ConstraintType::Equal,
//         ];
//         let rows = vec![
//             Row { name: "".to_string(), constraint_type: ConstraintType::Greater,},
//             Row { name: "".to_string(), constraint_type: ConstraintType::Equal,},
//         ];
//         let original_nr_rows = 2;
//         let ranges = vec![(0, R64!(2))];
//         let b = compute_b(rhss, &constraints, &rows, original_nr_rows, ranges);
//         assert_eq!(b, Ok(Dense::<T>::new(vec![R64!(1), R64!(3), R64!(5)], 3)));
//     }
// }
