//! # Organizing data read
//!
//! Organizing data read in the `parsing` module, and checking the linear program for consistency.
//! Contains also the definition of the struct representing linear programs in MPS format.
use std::cmp::Ordering;
use std::collections::HashSet;
use std::convert::TryInto;
use std::fmt::{Debug, Display};
use std::ops::{Add, Neg, Sub};

use num_traits::{One, Zero};
use relp_num::Abs;
use relp_num::NonZero;

use crate::data::linear_algebra::matrix::{ColumnMajor, Order as MatrixOrder, Sparse};
use crate::data::linear_algebra::SparseTuple;
use crate::data::linear_algebra::traits::Element;
use crate::data::linear_algebra::vector::{DenseVector, Vector};
use crate::data::linear_program::elements::{ConstraintRelation, RangedConstraintRelation, VariableType};
use crate::data::linear_program::general_form::GeneralForm;
use crate::data::linear_program::general_form::Variable as ShiftedVariable;
use crate::io::error::Inconsistency;
use crate::io::mps::{Bound, MPS, Range};
use crate::io::mps::BoundType;
use crate::io::mps::Column;
use crate::io::mps::Rhs;
use crate::io::mps::Row;

impl<FI, FO: From<FI>> TryInto<GeneralForm<FO>> for MPS<FI>
where
    FI: Sub<Output=FI> + Abs + Ord + Zero + Display + Clone,
    FO: NonZero + Zero + One + Neg<Output=FO> + Ord + Element,
    for<'r> FO: Add<&'r FO, Output=FO>,
    for<'r> &'r FO: Neg<Output=FO>,
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
    /// TODO(DOCUMENTATION): When can errors occur?
    fn try_into(self) -> Result<GeneralForm<FO>, Self::Error> {
        let (variable_info, columns, variable_names) = compute_variable_info(
            self.columns,
            self.cost_values,
            self.bounds,
            self.rows.len(),
        )?;
        let (constraint_types, b) = compute_constraint_info(
            self.rows,
            self.rhss,
            self.ranges,
        )?;

        Ok(GeneralForm::new(
            self.objective,
            columns,
            constraint_types,
            b,
            variable_info,
            variable_names,
            FO::zero(),
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
/// * `nr_rows`: Amount of rows, excluding the cost row, read.
///
/// # Return value
///
/// Collection of variable info's as used in the `GeneralForm` and a collection of the column names.
///
/// # Errors
///
/// If there is an inconsistency in bound information, such as a trivial infeasibility.
fn compute_variable_info<FI, FO: From<FI> + NonZero + Zero + One + Ord + Display + Debug + Clone>(
    columns: Vec<Column<FI>>,
    cost_values: Vec<SparseTuple<FI>>,
    bounds: Vec<Bound<FI>>,
    nr_rows: usize,
) -> Result<(Vec<ShiftedVariable<FO>>, Sparse<FO, FO, ColumnMajor>, Vec<String>), Inconsistency> {
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
            let cost = cost_values.next_if(|&(jj, _)| jj == j)
                .map_or(FO::zero(), |(_, v)| v.into());

            (ShiftedVariable {
                variable_type: variable.variable_type,
                cost,
                upper_bound: None,
                lower_bound: None,
                shift: FO::zero(),
                flipped: false,
            }, (variable.values, variable.name))
        }).unzip();
    process_bounds(&mut info, bounds)?;
    let (values, names): (Vec<_>, Vec<_>) = values_names.into_iter()
        .map(|(collection, names)| {
            (collection.into_iter().map(|(i, v)| (i, v.into())).collect(), names)
        })
        .unzip();

    let columns = ColumnMajor::new(values, nr_rows, names.len());

    Ok((info, columns, names))
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
/// TODO(CORRECTNESS): Consider changing this into an "this LP is infeasible" return type
fn process_bounds<FI, FO: From<FI> + Zero + One + Ord + Clone>(
    variable_info: &mut Vec<ShiftedVariable<FO>>,
    bounds: Vec<Bound<FI>>,
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
/// TODO(CORRECTNESS): Consider changing this into an "this LP is infeasible" return type
fn process_bound<FI, FO: From<FI> + Ord + Zero + One + Clone>(
    bound_type: BoundType<FI>,
    variable: &mut ShiftedVariable<FO>,
) -> Result<(bool, bool), Inconsistency> {
    match bound_type {
        BoundType::LowerContinuous(value) => {
            replace_existing_with(&mut variable.lower_bound, value.into(), Ordering::Greater);
            Ok((false, false))
        },
        BoundType::UpperContinuous(value) => {
            // The implied zero lower bound gets filled in only if no other lower bound is present.
            // This behavior is copied from GLPK.
            replace_existing_with(&mut variable.upper_bound, value.into(), Ordering::Less);
            Ok((true, false))
        },
        BoundType::Fixed(value) => {
            let converted: FO = value.into();
            replace_existing_with(&mut variable.lower_bound, converted.clone(), Ordering::Greater);
            replace_existing_with(&mut variable.upper_bound, converted, Ordering::Less);
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
            replace_existing_with(&mut variable.upper_bound, FO::zero(), Ordering::Less);
            Ok((false, false))
        },
        BoundType::UpperInfinity => {
            replace_existing_with(&mut variable.lower_bound, FO::zero(), Ordering::Greater);
            // Upper bound infinity is implied by variable.upper_bound.is_none()
            Ok((false, false))
        },
        BoundType::Binary => {
            replace_existing_with(&mut variable.lower_bound, FO::zero(), Ordering::Greater);
            replace_existing_with(&mut variable.upper_bound, FO::one(), Ordering::Less);
            variable.variable_type = VariableType::Integer;
            Ok((false, false))
        }
        BoundType::LowerInteger(value) => {
            replace_existing_with(&mut variable.lower_bound, value.into(), Ordering::Greater);
            variable.variable_type = VariableType::Integer;
            Ok((false, false))
        },
        BoundType::UpperInteger(value) => {
            replace_existing_with(&mut variable.upper_bound, value.into(), Ordering::Less);
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
fn replace_existing_with<OF: Ord>(option: &mut Option<OF>, new_value: OF, ordering: Ordering) {
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
fn fill_in_default_lower_bounds<F: Zero>(
    variables: &mut Vec<ShiftedVariable<F>>,
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
/// * `rows`: Direction of the constraint (name is not used).
/// * `rhss`: Right-hand sides to be converted into a `b`.
/// * `ranges`: Ranges read (will be flattened, grouping does not matter).
// TODO(ARCHITECTURE): Simplify these trait bounds
fn compute_constraint_info<
    FI: Sub<Output=FI> + Abs + Ord + Zero + Display + Clone,
    FO: From<FI> + Zero + Ord + Abs + Element + Clone,
>(
    rows: Vec<Row>,
    rhss: Vec<Rhs<FI>>,
    ranges: Vec<Range<FI>>,
) -> Result<(Vec<RangedConstraintRelation<FO>>, DenseVector<FO>), Inconsistency>
where
    for<'r> FO: Add<&'r FO, Output=FO>,
    for<'r> &'r FO: Neg<Output=FO>,
{
    let ranges = compute_ranges(&rhss, ranges, rows.len())?;
    let mut constraint_types = compute_constraint_types(&rows, ranges);
    let b = compute_b(rhss, &mut constraint_types, &rows, rows.len())?;

    Ok((constraint_types, b))
}

/// Flatten the ranges and ensure consistency related to them.
///
/// Consistency: e.g. there should be at most one range per row, and if there is one, there should
/// be only one right hand side (or they should be equal).
fn compute_ranges<F: PartialEq>(
    rhss: &[Rhs<F>],
    ranges: Vec<Range<F>>,
    nr_rows: usize,
) -> Result<Vec<(usize, F)>, Inconsistency> {
    // If there are no ranges, this is not an issue anyway
    if ranges.is_empty() {
        return Ok(Vec::with_capacity(0))
    }

    // Flatten, we don't care about the different range names
    let mut range_rows = ranges.into_iter()
        .flat_map(|range| range.values.into_iter())
        .collect::<Vec<_>>();
    // We process them by row.
    range_rows.sort_unstable_by_key(|&(i, _)| i);
    let unduplicated_length = range_rows.len();
    range_rows.dedup_by_key(|&mut (i, _)| i);
    if range_rows.len() < unduplicated_length {
        return Err(Inconsistency::new("Only one range per row can be specified."));
    }

    let mut already_seen = vec![false; nr_rows];
    // It's very unlikely that this code path is necessary
    let mut duplicates = Vec::with_capacity(0);
    for rhs in rhss {
        for &(i, _) in &rhs.values {
            if already_seen[i] {
                duplicates.push(i);
            } else {
                already_seen[i] = true;
            }
        }
    }

    // See if any of the duplicate value have a range
    for duplicate in duplicates {
        if range_rows.iter().any(|&(i, _)| i == duplicate) {
            let values = rhss.iter().flat_map(|rhs| rhs.values.iter())
                .filter_map(|(i, v)| if *i == duplicate { Some(v) } else { None })
                .collect::<Vec<_>>();
            if let Some(first) = values.first() {
                if values.iter().any(|v| v != first) {
                    return Err(Inconsistency::new("Multiple rhs values for a constraint with a range"));
                }
            }
        }
    }

    Ok(range_rows)
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
fn compute_constraint_types<FI: Zero + PartialEq, FO: From<FI>>(
    rows: &[Row],
    ranges: Vec<SparseTuple<FI>>,
) -> Vec<RangedConstraintRelation<FO>> {
    // Sorted and unique
    debug_assert!(ranges.windows(2).all(|w| w[0].0 < w[1].0));
    debug_assert!(ranges.iter().all(|&(i, _)| i < rows.len()));

    let mut ranges = ranges.into_iter().peekable();
    rows.iter().enumerate()
        .map(|(i, row)| {
            if let Some((_, r)) = ranges.next_if(|&(ii, _)| ii == i) {
                if r == FI::zero() {
                    RangedConstraintRelation::Equal
                } else {
                    RangedConstraintRelation::Range(r.into())
                }
            } else {
                row.constraint_type.into()
            }
        })
        .collect()
}

/// Combine all constraint values.
///
/// We mostly just take the tightest bound and process the ranges.
///
/// # Arguments
///
/// * `rhss`: Right hand sides (often only one), b's values.
/// * `constraints`: Constraint directions.
/// * `rows`: Original constraint types (used for processing ranges).
/// * `nr_rows`: Number of constraints without ranges.
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
fn compute_b<
    OFI: Sub<Output=OFI> + Abs + PartialOrd + Zero + Clone,
    OFO: From<OFI> + Zero + Abs + PartialOrd + Element + Display,
>(
    rhss: Vec<Rhs<OFI>>,
    constraints: &mut [RangedConstraintRelation<OFO>],
    rows: &[Row],
    nr_rows: usize,
) -> Result<DenseVector<OFO>, Inconsistency>
where
    for<'r> OFO: Add<&'r OFO, Output=OFO>,
    for<'r> &'r OFO: Neg<Output=OFO>,
{
    debug_assert!(rhss.iter().all(|rhs| rhs.values.is_sorted_by_key(|&(i, _)| i)));
    debug_assert!(rhss.iter().all(|rhs| rhs.values.iter().all(|&(i, _)| i < nr_rows)));
    // Assumption: if there are multiple rhs values for a row with a range, all these values are
    // equal

    // We fill b with options, and then explicitly substitute the default value later.
    let mut b: Vec<Option<OFO>> = vec![None; rows.len()];
    for (i, value) in rhss.into_iter().flat_map(|rhs| rhs.values.into_iter()) {
        match &mut b[i] {
            None => match &mut constraints[i] {
                RangedConstraintRelation::Range(range) => {
                    let range_sign = range.cmp(&&mut OFO::zero());
                    if range < &mut OFO::zero() {
                        *range = -&*range;
                    }
                    let value = value.into();
                    let bound  = match (rows[i].constraint_type, range_sign) {
                        (ConstraintRelation::Greater, _) => value + &*range,
                        (ConstraintRelation::Less, _) => value,
                        (ConstraintRelation::Equal, Ordering::Greater | Ordering::Equal) => value + &*range,
                        (ConstraintRelation::Equal, Ordering::Less | Ordering::Equal) => value,
                    };
                    b[i] = Some(bound);
                },
                _ => b[i] = Some(value.into()),
            }
            Some(current) => {
                // There can be at most one rhs value for each row that has a range.
                debug_assert!(!matches!(constraints[i], RangedConstraintRelation::Range(_)));

                let converted: OFO = value.into();
                match rows[i].constraint_type {
                    ConstraintRelation::Equal => if &converted != current {
                        return Err(Inconsistency::new(
                            format!("Trivial infeasibility: a constraint can't equal both {} and {}", current, converted),
                        ))
                    },
                    ConstraintRelation::Greater => if &converted > current {
                        *current = converted;
                    },
                    ConstraintRelation::Less => if &converted < current {
                        *current = converted;
                    },
                }
            }
        }
    }

    // Substitute the default value.
    let b = b.into_iter().map(|value| value.unwrap_or_else(OFO::zero)).collect();
    Ok(DenseVector::new(b, nr_rows))
}

#[cfg(test)]
#[allow(clippy::shadow_unrelated)]
mod test {
    use relp_num::R32;
    use relp_num::Rational32;

    use crate::data::linear_algebra::vector::{DenseVector, Vector};
    use crate::data::linear_program::elements::{ConstraintRelation, RangedConstraintRelation};
    use crate::io::mps::{Rhs, Row};
    use crate::io::mps::convert::compute_b;

    type T = Rational32;

    #[test]
    fn test_compute_b() {
        // No ranges, no data
        let rhss: Vec<Rhs<Rational32>> = vec![];
        let mut constraints = vec![];
        let rows = vec![];
        let original_nr_rows = 0;
        let b = compute_b(rhss, &mut constraints, &rows, original_nr_rows);
        assert_eq!(b, Ok(DenseVector::<T>::new(vec![], 0)));

        // No ranges, one rhs
        let rhss = vec![Rhs { name: "R".to_string(), values: vec![(0, R32!(1))], }];
        let mut constraints = vec![RangedConstraintRelation::Equal];
        let rows = vec![Row { name: "".to_string(), constraint_type: ConstraintRelation::Equal}];
        let original_nr_rows = 1;
        let b = compute_b(rhss, &mut constraints, &rows, original_nr_rows);
        assert_eq!(b, Ok(DenseVector::<T>::new(vec![R32!(1)], 1)));

        // No ranges, two rhses
        let rhss = vec![
            Rhs { name: "R1".to_string(), values: vec![(0, R32!(1))], },
            Rhs { name: "R2".to_string(), values: vec![(0, R32!(2))], },
        ];
        let mut constraints = vec![RangedConstraintRelation::Greater];
        let rows = vec![Row { name: "".to_string(), constraint_type: ConstraintRelation::Greater}];
        let original_nr_rows = 1;
        let b = compute_b(rhss, &mut constraints, &rows, original_nr_rows);
        assert_eq!(b, Ok(DenseVector::<T>::new(vec![R32!(2)], 1)));

        // One range with data before
        let rhss = vec![
            Rhs { name: "R".to_string(), values: vec![(0, R32!(1)), (1, R32!(5))], },
        ];
        let mut constraints = vec![
            RangedConstraintRelation::Greater,
            RangedConstraintRelation::Equal,
        ];
        let rows = vec![
            Row { name: "".to_string(), constraint_type: ConstraintRelation::Greater,},
            Row { name: "".to_string(), constraint_type: ConstraintRelation::Equal,},
        ];
        let original_nr_rows = 2;
        let b = compute_b(rhss, &mut constraints, &rows, original_nr_rows);
        assert_eq!(b, Ok(DenseVector::<T>::new(vec![R32!(1), R32!(5)], 2)));

        // One range with data after
        let rhss = vec![
            Rhs { name: "R".to_string(), values: vec![(0, R32!(1)), (1, R32!(5))], },
        ];
        let mut constraints = vec![
            RangedConstraintRelation::Greater,
            RangedConstraintRelation::Equal,
        ];
        let rows = vec![
            Row { name: "".to_string(), constraint_type: ConstraintRelation::Greater,},
            Row { name: "".to_string(), constraint_type: ConstraintRelation::Equal,},
        ];
        let original_nr_rows = 2;
        let b = compute_b(rhss, &mut constraints, &rows, original_nr_rows);
        assert_eq!(b, Ok(DenseVector::<T>::new(vec![R32!(1), R32!(5)], 2)));
    }
}
