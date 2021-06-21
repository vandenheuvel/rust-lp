use std::cmp::max;
use std::collections::HashSet;
use std::hash::Hash;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub};

use num_traits::{One, Zero};
use num_traits::Pow;
use relp_num::{Abs, NonZeroFactorizable, NonZeroFactorization, Sign, Signed};

use crate::data::linear_algebra::SparseTuple;
use crate::data::linear_algebra::traits::{SparseComparator, SparseElement};
use crate::data::linear_program::general_form::GeneralForm;
use crate::data::linear_program::general_form::presolve::scale::{Scalable, Scaling};
use std::iter;
impl<R> Scalable<R> for GeneralForm<R>
where
    for<'r> R:
        NonZeroFactorizable +
        MulAssign<R::Factor> +
        DivAssign<R::Factor> +
        MulAssign<&'r R::Factor> +
        DivAssign<&'r R::Factor> +
        MulAssign<R> +
        DivAssign<R> +
        MulAssign<&'r R> +
        DivAssign<&'r R> +
        Mul<&'r R, Output=R> +
        One +
        SparseElement<R> +
        SparseComparator +
    ,
    for<'r> &'r R: Mul<&'r R, Output=R>,
    R::Factor: Hash + Pow<R::Power, Output=R::Factor>,
    R::Power: Neg<Output=R::Power> + Sub<Output=R::Power> + Add<Output=R::Power> + Zero + One + AddAssign + Ord,
    for<'r> &'r R: Div<&'r R, Output=R>,
{
    fn scale(&mut self) -> Scaling<R> {
        // Factorize all numbers at once
        let factorization = self.factorize();
        // Optimize per factor
        let scale_per_factor = factorization.solve();
        // Aggregate scaling per factor
        let scaling: Scaling<R> = total_scaling(scale_per_factor);

        // Apply the scaling
        let Scaling {
            c_factor,
            constraint_row_factors,
            b_factor,
            constraint_column_factors,
        } = &scaling;
        for variable in &mut self.variables {
            variable.cost *= c_factor;
        }
        for rhs in &mut self.b.data {
            *rhs *= b_factor;
        }
        for (j, column) in self.constraints.data.iter_mut().enumerate() {
            for (i, value) in column {
                // The row_scales value is large, the column_scales value is small, so combine them
                // first, some factors might cancel
                *value *= &constraint_row_factors[*i] * &constraint_column_factors[j];
            }
        }

        scaling
    }

    fn scale_back(&mut self, scaling: Scaling<R>) {
        for variable in &mut self.variables {
            variable.cost /= &scaling.c_factor;
        }
        for rhs in &mut self.b.data {
            *rhs /= &scaling.b_factor;
        }
        for (j, column) in self.constraints.data.iter_mut().enumerate() {
            for (i, value) in column {
                // The row_scales value is large, the column_scales value is small
                *value /= &scaling.constraint_column_factors[j] * &scaling.constraint_row_factors[*i];
            }
        }
    }
}

#[derive(Eq, PartialEq, Clone, Debug)]
struct ColumnInfo<P> {
    /// Signed, needs at most 32 bits for arbitrary precision numbers
    lowest_value: P,
    /// Relates to the number of rows
    lowest_count: u32,
    /// Signed, needs at most 32 bits for arbitrary precision numbers
    highest_value: P,
}

#[derive(Eq, PartialEq, Debug)]
struct Factorization<R: NonZeroFactorizable> {
    /// No duplicate elements
    all_factors: Vec<R::Factor>,
    /// Zero values are None, others have a factorization that might be empty
    b: Vec<Option<Vec<(R::Factor, R::Power)>>>,
    c: Vec<Option<Vec<(R::Factor, R::Power)>>>,
    /// Column major
    A: Vec<Vec<SparseTuple<Vec<(R::Factor, R::Power)>>>>,
}

impl<R: NonZeroFactorizable<Power: Zero + Ord>> Factorization<R> {
    fn solve(mut self) -> Vec<(R::Factor, ((R::Power, Vec<R::Power>), (R::Power, Vec<R::Power>)))> {
        // Row major column indices
        let index = {
            let mut index = vec![Vec::new(); self.A.len()];
            for (j, column) in self.A.iter().enumerate() {
                for (data_index, &(i, _)) in column.iter().enumerate() {
                    index[i].push((j, data_index));
                }
            }
            index
        };

        let mut results = Vec::with_capacity(self.all_factors.len());
        debug_assert!(self.all_factors.is_sorted());
        while !self.all_factors.is_empty() {
            let (factor, (row_scaling, column_scaling)) = self.solve_single(&index);
            results.push((factor, (row_scaling, column_scaling)));
        }

        results
    }
    fn solve_single(&mut self, index: &Vec<Vec<(usize, usize)>>) -> (R::Factor, ((R::Power, Vec<R::Power>), (R::Power, Vec<R::Power>))) {
        let (column_info, gains) = self.solve_single_setup(index);
        let row_increments = self.solve_single_rows(gains, column_info, index);
        let column_changes = self.solve_single_columns(&row_increments);
        let factor = self.remove_factor_info();

        (factor, (row_increments, column_changes))
    }

    fn solve_single_setup(
        &self, index: &Vec<Vec<(usize, usize)>>,
    ) -> ((Option<ColumnInfo<R::Power>>, Vec<ColumnInfo<R::Power>>), ((usize, usize), Vec<(usize, usize)>)) {
        let b_info = {
            let mut extreme_values = None;
            let mut lowest_count = 0;

            for value in &self.b {
                if let Some(factorization) = value {
                    let exponent = self.initial_exponent(factorization);
                    match &mut extreme_values {
                        None => {
                            extreme_values = Some((exponent, exponent));
                            lowest_count = 1;
                        },
                        Some((lowest, highest)) => {
                            if exponent < *lowest {
                                *lowest = exponent;
                                lowest_count = 0;
                            }
                            if exponent == *lowest {
                                lowest_count += 1;
                            }
                            *highest = max(exponent, *highest);
                        }
                    }
                }
            }

            extreme_values.map(|(lowest_value, highest_value)| {
                ColumnInfo { lowest_value, lowest_count, highest_value }
            })
        };

        let infos = self.A.iter().enumerate().map(|(j, column)| {
                let first_power = self.initial_exponent(&column[0].1);
                let mut lowest_value = first_power;
                let mut lowest_count = 1;
                let mut highest_value = first_power;

                let mut update = |power| {
                    if power < lowest_value {
                        lowest_value = power;
                        lowest_count = 0;
                    }
                    if power == lowest_value {
                        lowest_count += 1;
                    }
                    highest_value = max(power, highest_value);
                };

                for (_row_index, factorization) in &column[1..] {
                    let power = self.initial_exponent(factorization);
                    update(power);
                }

                if let Some(factorization) = &self.c[j] {
                    let power = self.initial_exponent(factorization);
                    update(power);
                }

                ColumnInfo { lowest_value, lowest_count, highest_value }
            })
            .collect::<Vec<_>>();

        let (lower, upper): (Vec<_>, Vec<_>) = self.c.iter().enumerate()
            .filter_map(|(j, v)| v.as_ref().map(|x| (j, x)))
            .map(|(j, factorization)| {
                let value = self.initial_exponent(factorization);

                let info = &infos[j];

                let lower_relevance = value == info.lowest_value && info.lowest_count == 1;
                let upper_relevance = value == info.highest_value;

                (lower_relevance, upper_relevance)
            }).unzip();
        let lower_score = lower.into_iter().filter(|&v| v).count();
        let upper_score = upper.into_iter().filter(|&v| v).count();
        let cost_gains = (lower_score, upper_score);

        let (b_weight, normal_column_weight) = self.relative_column_weight();
        let row_gains = index.iter().enumerate().map(|(row_index, cols)| {
            let (lower_b_score, upper_b_score) = match &self.b[row_index] {
                Some(factorization) => {
                    let info = b_info.as_ref()
                        .expect("There is a nonzero value in b, so a maximum and minimum count should be present");

                    let power = self.initial_exponent(factorization);
                    let lower = if power == info.lowest_value && info.lowest_count == 1 { 1 } else { 0 };
                    let upper = if power == info.highest_value { 1 } else { 0 };

                    (lower, upper)
                }
                None => (0, 0),
            };

            let (lower, upper): (Vec<_>, Vec<_>) = cols.iter().map(|&(col_index, data_index)| {
                let exponent = self.initial_exponent(&self.A[col_index][data_index].1);
                let info = &infos[col_index];

                let lower_relevance = exponent == info.lowest_value && info.lowest_count == 1;
                let upper_relevance = exponent == info.highest_value;

                (lower_relevance, upper_relevance)
            }).unzip();
            let lower_column_score = lower.into_iter().filter(|&v| v).count();
            let upper_column_score = upper.into_iter().filter(|&v| v).count();

            let lower_score = lower_b_score * b_weight + lower_column_score * normal_column_weight;
            let upper_score = upper_b_score * b_weight + upper_column_score * normal_column_weight;

            (lower_score, upper_score)
        }).collect();

        ((b_info, infos), (cost_gains, row_gains))
    }

    /// Compute the relative importance of the rhs column w.r.t. the constraint columns.
    ///
    /// Call the number of constraints `m`, the number of variables `n`. A typical problem has
    /// (much) more columns than constraints, so `m << n`. Columns will on average be in the basis
    /// in `m` out of `n` times, so we weigh the cost of having prime factors in those columns
    /// accordingly: weight `n` for the rhs column and smaller weight `m` for the constraint
    /// columns.
    ///
    /// # Return value
    ///
    /// A tuple weights, first element for the rhs vector, second element for the constraint
    /// vectors.
    fn relative_column_weight(&self) -> (usize, usize) {
        let (nr_rows, nr_columns) = self.dimensions();

        let normal_column_weight = nr_rows;
        let b_weight = nr_columns;

        (b_weight, normal_column_weight)
    }

    fn solve_single_rows(
        &self,
        gains: ((usize, usize), Vec<(usize, usize)>),
        column_info: (Option<ColumnInfo<R::Power>>, Vec<ColumnInfo<R::Power>>),
        index: &Vec<Vec<(usize, usize)>>,
    ) -> (R::Power, Vec<R::Power>) {
        let (nr_rows, _) = self.dimensions();

        let ((mut cost_gain, mut cost_penalty), mut constraint_row_gains) = gains;
        let (mut rhs_info, mut constraint_column_info) = column_info;

        let mut cost_row_increment = R::Power::zero();
        let mut constraint_row_increments = vec![R::Power::zero(); nr_rows];
        let mut not_incremented = 1 + nr_rows;
        while let Some(to_increment) = Self::determine_next_row(
            not_incremented,
            cost_row_increment,
            &constraint_row_increments,
            &(cost_gain, cost_penalty),
            &constraint_row_gains,
        ) {
            match to_increment {
                RowToIncrement::CostRow => {
                    self.increment_cost_row(
                        &mut cost_gain, &mut cost_penalty,
                        &mut constraint_column_info,
                        cost_row_increment,
                    );

                    cost_row_increment += R::Power::one();
                }
                RowToIncrement::ConstraintRow(row_index) => {
                    self.increment_constraint_row(
                        row_index,
                        &mut rhs_info, &mut constraint_column_info,
                        &mut constraint_row_gains,
                        &constraint_row_increments,
                        index,
                    );

                    let increments = &mut constraint_row_increments[row_index];
                    if increments.is_zero() {
                        not_incremented -= 1;
                    }
                    *increments += R::Power::one();
                }
            }
        }

        (cost_row_increment, constraint_row_increments)
    }

    fn determine_next_row(
        not_incremented: usize,
        cost_row_increment: R::Power,
        constraint_row_increments: &Vec<R::Power>,
        cost_row_gain: &(usize, usize),
        constraint_row_gains: &Vec<(usize, usize)>,
    ) -> Option<RowToIncrement> {
        constraint_row_gains.iter().enumerate()
            .map(|(row_index, gain)| (RowToIncrement::ConstraintRow(row_index), gain))
            .chain(iter::once((RowToIncrement::CostRow, cost_row_gain)))
            .filter(|(_, (gain, penalty))| gain >= penalty)
            .max_by_key(|(row, (gain, penalty))| {
                (gain - penalty, match *row {
                    RowToIncrement::CostRow => cost_row_increment.is_positive(),
                    RowToIncrement::ConstraintRow(index) => constraint_row_increments[index].is_positive(),
                })
            })
            .filter(|(row, _)| {
                not_incremented > 1 || match *row {
                    RowToIncrement::CostRow => cost_row_increment.is_positive(),
                    RowToIncrement::ConstraintRow(index) => constraint_row_increments[index].is_positive(),
                }
            })
            .map(|(row, _)| row)
    }

    fn increment_cost_row(
        &self,
        gain: &mut usize, penalty: &mut usize,
        constraint_column_info: &mut Vec<ColumnInfo<R::Power>>,
        cost_row_increment: R::Power,
    ) {
        for (j, cost_value) in self.c.iter().enumerate() {
            if let Some(factorization) = cost_value {
                let initial_exponent = self.initial_exponent(factorization);
                let exponent = initial_exponent + cost_row_increment;
                self.update_info_with_exponent(
                    exponent,
                    self.A[j].iter().map(|(_, factorization)| factorization),
                    &mut constraint_column_info[j],
                    gain, penalty,
                );
            }
        }
    }

    fn increment_constraint_row(
        &self,
        row_index: usize,
        rhs_info: &mut Option<ColumnInfo<R::Power>>, constraint_column_info: &mut Vec<ColumnInfo<R::Power>>,
        constraint_row_gains: &mut Vec<(usize, usize)>,
        constraint_row_increments: &Vec<R::Power>,
        index: &Vec<Vec<(usize, usize)>>,
    ) {
        let (gain, penalty) = &mut constraint_row_gains[row_index];
        if let Some(factorization) = &self.b[row_index] {
            let initial_exponent = self.initial_exponent(factorization);
            let exponent = initial_exponent + constraint_row_increments[row_index];
            self.update_info_with_exponent(
                exponent,
                self.b.iter().filter_map(Option::as_ref),
                rhs_info.as_mut().unwrap(),
                gain, penalty,
            );
        }

        for &(column, row_data_index) in &index[row_index] {
            let (_, factorization) = &self.A[column][row_data_index];
            let initial_exponent = self.initial_exponent(factorization);
            let exponent = initial_exponent + constraint_row_increments[row_index];
            self.update_info_with_exponent(
                exponent,
                self.A[column].iter().map(|(_, factorization)| factorization),
                &mut constraint_column_info[column],
                gain, penalty,
            );
        }
    }

    fn initial_exponent(
        &self,
        factorization: &Vec<(R::Factor, R::Power)>,
    ) -> R::Power {
        let factor = self.all_factors.last().unwrap();

        factorization.last()
            .filter(|(f, _)| f == factor)
            .map(|&(_, power)| power)
            .unwrap_or(R::Power::zero())
    }

    fn update_info_with_exponent<'a>(
        &'a self,
        exponent: R::Power,
        column: impl Iterator<Item = &'a Vec<(R::Factor, R::Power)>>,
        info: &mut ColumnInfo<R::Power>,
        gain: &mut usize, penalty: &mut usize,
    ) {
        if exponent == info.lowest_value {
            info.lowest_count -= 1;
            if info.lowest_count == 0 {
                *gain -= 1;

                // A new lost value, scan for the count
                info.lowest_value += R::Power::one();
                for factorization in column {
                    let power = self.initial_exponent(factorization);

                    if power == info.lowest_value {
                        info.lowest_count += 1;
                    }
                }
            }
        }
        if exponent == info.highest_value - R::Power::one() {
            *penalty += 1;
        }
        if exponent == info.highest_value {
            info.highest_value += R::Power::one();
            // No need to adapt gains
        }
    }

    fn solve_single_columns(&self, row_increments: &(R::Power, Vec<R::Power>)) -> (R::Power, Vec<R::Power>) {
        let factor = self.all_factors.last().unwrap();
        let (cost_row_increment, constraint_row_increments) = row_increments;

        let mut rhs_powers = self.b.iter().enumerate()
            .filter(|(_, factorization)| factorization.is_some())
            .map(|(row_index, factorization)| (row_index, factorization.as_ref().unwrap()))
            .map(|(row_index, factorization)| {
                constraint_row_increments[row_index] + match factorization.last() {
                    Some((f, power)) if f == factor => *power,
                    _ => R::Power::zero(),
                }
            }).collect::<Vec<_>>();
        rhs_powers.sort_unstable();
        let middle = rhs_powers.len() / 2;
        debug_assert!(!rhs_powers.is_empty(), "At least one element per column (even if zero).");
        let rhs_middle_power = rhs_powers[middle];

        let constraint_columns = self.A.iter().enumerate().map(|(j, column)| {
            let mut powers = column.iter()
                .map(|(row_index, factorization)| (constraint_row_increments[*row_index], factorization))
                .chain(self.c[j].iter().map(|factorization| (*cost_row_increment, factorization)))
                .map(|(increment, factorization)| {
                increment + match factorization.last() {
                    Some((f, power)) if f == factor => *power,
                    _ => R::Power::zero(),
                }
            }).collect::<Vec<_>>();
            powers.sort_unstable();
            let middle = powers.len() / 2;
            debug_assert!(!powers.is_empty(), "At least one element per column (even if zero).");
            powers[middle]
        }).collect();

        (rhs_middle_power, constraint_columns)
    }

    fn remove_factor_info(&mut self) -> R::Factor {
        let factor = self.all_factors.pop().unwrap();

        for column in &mut self.A {
            for (_, factorization) in column {
                match factorization.last() {
                    Some((f, _)) if f == &factor => {factorization.pop();},
                    _ => {},
                }
            }
        }

        let remove_from = |data: &mut [Option<Vec<(R::Factor, R::Power)>>]| {
            for factorization in data {
                if let Some(factorization) = factorization {
                    match factorization.last() {
                        Some((f, _)) if f == &factor => {factorization.pop();},
                        _ => {},
                    }
                }
            }
        };

        remove_from(&mut self.b);
        remove_from(&mut self.c);

        factor
    }

    fn dimensions(&self) -> (usize, usize) {
        let nr_rows = self.b.len();
        let nr_columns = self.c.len();

        (nr_rows, nr_columns)
    }
}

enum RowToIncrement {
    CostRow,
    ConstraintRow(usize),
}

impl<R> GeneralForm<R>
where
    R: NonZeroFactorizable<Factor: Hash> + SparseElement<R> + SparseComparator,
{
    fn factorize(&self) -> Factorization<R> {
        let mut all_factors = HashSet::new();

        let b = self.b.data.iter()
            .map(|v| {
                if v.is_not_zero() {
                    let NonZeroFactorization { factors, .. } = v.factorize();
                    all_factors.extend(factors.iter().map(|(f, _)| f).cloned());
                    Some(factors)
                } else {
                    None
                }
            })
            .collect();

        let c = self.variables.iter()
            .map(|variable| {
                let NonZeroFactorization { factors, .. } = variable.cost.factorize();
                all_factors.extend(factors.iter().map(|(f, _)| f).cloned());
                Some(factors)
            })
            .collect();

        let A = self.constraints.iter_columns().map(|column| {
            column.iter().map(|(i, v)| {
                let NonZeroFactorization { factors, .. } = v.factorize();
                all_factors.extend(factors.iter().map(|(f, _)| f).cloned());
                (*i, factors)
            }).collect()
        }).collect();

        let mut all_factors = all_factors.into_iter().collect::<Vec<_>>();
        all_factors.sort_unstable();
        Factorization { all_factors, b, c, A }
    }
}

fn total_scaling<R>(
    scale_per_factor: Vec<(R::Factor, ((R::Power, Vec<R::Power>), (R::Power, Vec<R::Power>)))>,
) -> Scaling<R>
where
    R: NonZeroFactorizable<Power: Abs>,
    for<'r> R: MulAssign<&'r R::Factor> + DivAssign<&'r R::Factor>,
{
    debug_assert!(!scale_per_factor.is_empty());

    let nr_rows = scale_per_factor[0].1.0.1.len();
    let nr_columns = scale_per_factor[0].1.1.1.len();

    let mut c_factor = R::one();
    let mut row_scales = vec![R::one(); nr_rows];
    let mut b_factor = R::one();
    let mut column_scales = vec![R::one(); nr_columns];

    for (factor, ((c_change, row_changes), (b_change, column_changes))) in scale_per_factor {
        match c_change.signum() {
            Sign::Positive => {
                let mut iter = R::Power::zero();
                while iter < c_change.abs() {
                    c_factor *= &factor;
                    iter += R::Power::one();
                }
            }
            Sign::Zero => {}
            Sign::Negative => {
                let mut iter = R::Power::zero();
                while iter < c_change.abs() {
                    c_factor /= &factor;
                    iter += R::Power::one();
                }
            }
        }
        for (i, row_change) in row_changes.into_iter().enumerate() {
            match row_change.signum() {
                Sign::Positive => {
                    let mut iter = R::Power::zero();
                    while iter < row_change.abs() {
                        row_scales[i] *= &factor;
                        iter += R::Power::one();
                    }
                },
                Sign::Zero => {},
                Sign::Negative => {
                    let mut iter = R::Power::zero();
                    while iter < row_change.abs() {
                        row_scales[i] /= &factor;
                        iter += R::Power::one();
                    }
                },
            }
        }
        match b_change.signum() {
            Sign::Positive => {
                let mut iter = R::Power::zero();
                while iter < b_change.abs() {
                    b_factor *= &factor;
                    iter += R::Power::one();
                }
            }
            Sign::Zero => {}
            Sign::Negative => {
                let mut iter = R::Power::zero();
                while iter < b_change.abs() {
                    b_factor /= &factor;
                    iter += R::Power::one();
                }
            }
        }
        for (j, column_change) in column_changes.into_iter().enumerate() {
            match column_change.signum() {
                Sign::Positive => {
                    let mut iter = R::Power::zero();
                    while iter < column_change.abs() {
                        column_scales[j] *= &factor;
                        iter += R::Power::one();
                    }
                },
                Sign::Zero => (),
                Sign::Negative => {
                    let mut iter = R::Power::zero();
                    while iter < column_change.abs() {
                        column_scales[j] /= &factor;
                        iter += R::Power::one();
                    }
                },
            }
        }
    }

    Scaling {
        c_factor,
        constraint_row_factors: row_scales,
        b_factor,
        constraint_column_factors: column_scales,
    }
}

#[cfg(test)]
mod test {
    use relp_num::{R8, Rational16, Rational8};

    use crate::data::linear_algebra::matrix::{ColumnMajor, Order};
    use crate::data::linear_algebra::vector::{DenseVector, Vector};
    use crate::data::linear_algebra::vector::test::TestVector;
    use crate::data::linear_program::elements::{Objective, RangedConstraintRelation, VariableType};
    use crate::data::linear_program::general_form::{GeneralForm, Variable};
    use crate::data::linear_program::general_form::presolve::scale::{Scalable, Scaling};
    use crate::data::linear_program::general_form::presolve::scale::rational::{ColumnInfo, Factorization};

    #[test]
    fn test_factorize() {
        let mut general_form: GeneralForm<Rational8> = GeneralForm::new(
            Objective::Minimize,
            ColumnMajor::from_test_data::<_, _, u8>(&[
                vec![1, 2],
                vec![4, 6],
                vec![7, 14],
                vec![0, 11],
            ], 2),
            vec![
                RangedConstraintRelation::Equal,
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Greater,
                RangedConstraintRelation::Equal,
            ],
            DenseVector::from_test_data::<u8>(vec![3, 0, 21, 11]),
            vec![
                Variable {
                    variable_type: VariableType::Continuous,
                    cost: R8!(4),
                    lower_bound: Some(R8!(0)),
                    upper_bound: Some(R8!(6)),
                    shift: R8!(0),
                    flipped: false
                },
                Variable {
                    variable_type: VariableType::Continuous,
                    cost: R8!(11),
                    lower_bound: Some(R8!(1)),
                    upper_bound: Some(R8!(2)),
                    shift: R8!(0),
                    flipped: false
                },
            ],
            vec!["x".to_string(), "y".to_string()],
            R8!(16),
        );

        let factorization = general_form.factorize();

        let expected = Factorization {
            all_factors: vec![2, 3, 7, 11].into_iter().collect(),
            b: vec![Some(vec![(3, 1)]), None, Some(vec![(3, 1), (7, 1)]), Some(vec![(11, 1)])],
            c: vec![Some(vec![(2, 2)]), Some(vec![(11, 1)])],
            A: vec![
                vec![(0, vec![]), (1, vec![(2, 2)]), (2, vec![(7, 1)])],
                vec![(0, vec![(2, 1)]), (1, vec![(2, 1), (3, 1)]), (2, vec![(2, 1), (7, 1)]), (3, vec![(11, 1)])],
            ],
        };

        assert_eq!(factorization, expected);
    }

    #[test]
    fn test_solve_single_setup_without_b() {
        let factorization = Factorization::<Rational16> {
            all_factors: vec![2, 3, 7, 11].into_iter().collect(),
            b: vec![None, None, None, None],
            c: vec![Some(vec![(11, 1)]), Some(vec![(2, 2)])],
            A: vec![
                vec![(0, vec![]), (1, vec![(2, 2)]), (2, vec![(7, 1)])],
                vec![(0, vec![(2, 1)]), (1, vec![(2, 1), (3, 1)]), (2, vec![(2, 1), (7, 1)]), (3, vec![(11, 1)])],
            ],
        };

        let index = vec![
            vec![(0, 0), (1, 0)],
            vec![(0, 1), (1, 1)],
            vec![(0, 2), (1, 2)],
            vec![(1, 3)],
        ];
        let result = factorization.solve_single_setup(&index);
        // For factor `11`
        let expected = (
            (
                None,
                vec![
                    ColumnInfo {
                        lowest_value: 0,
                        lowest_count: 3,
                        highest_value: 1,
                    },
                    ColumnInfo {
                        lowest_value: 0,
                        lowest_count: 4,
                        highest_value: 1,
                    },
                ],
            ),
            (
                (0, 1),
                vec![
                    (0, 0),
                    (0, 0),
                    (0, 0),
                    (0, 4),
                ],
            )
        );
        assert_eq!(result, expected);
    }

    #[test]
    fn test_solve_single_setup_with_b() {
        let factorization = Factorization::<Rational16> {
            all_factors: vec![2, 3, 7, 11].into_iter().collect(),
            b: vec![Some(vec![(3, 1)]), None, Some(vec![(3, 1), (7, 1)]), Some(vec![(11, 1)])],
            c: vec![Some(vec![(11, 1)]), Some(vec![(2, 2)])],
            A: vec![
                vec![(0, vec![]), (1, vec![(2, 2)]), (2, vec![(7, 1)])],
                vec![(0, vec![(2, 1)]), (1, vec![(2, 1), (3, 1)]), (2, vec![(2, 1), (7, 1)]), (3, vec![(11, 1)])],
            ],
        };

        let index = vec![
            vec![(0, 0), (1, 0)],
            vec![(0, 1), (1, 1)],
            vec![(0, 2), (1, 2)],
            vec![(1, 3)],
        ];
        let result = factorization.solve_single_setup(&index);
        // For the highest factor, so factor `11`
        let expected = (
            (
                Some(ColumnInfo {
                    lowest_value: 0,
                    lowest_count: 2,
                    highest_value: 1,
                }),
                vec![
                    ColumnInfo {
                        lowest_value: 0,
                        lowest_count: 3,
                        highest_value: 1,
                    },
                    ColumnInfo {
                        lowest_value: 0,
                        lowest_count: 4,
                        highest_value: 1,
                    },
                ],
            ),
            (
                (0, 1),
                vec![
                    (0, 0),
                    (0, 0),
                    (0, 0),
                    (0, 1 * 2 + 1 * 4),
                ],
            )
        );
        assert_eq!(result, expected);
    }

    #[test]
    fn test_solve_single_rows() {

    }

    #[test]
    fn test_scale() {
        let mut general_form: GeneralForm<Rational8> = GeneralForm::new(
            Objective::Minimize,
            ColumnMajor::from_test_data::<_, _, u8>(&[
                vec![1, 2],
                vec![4, 6],
                vec![7, 14],
            ], 2),
            vec![
                RangedConstraintRelation::Equal,
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Greater,
            ],
            DenseVector::from_test_data::<u8>(vec![3, 12, 21]),
            vec![
                Variable {
                    variable_type: VariableType::Continuous,
                    cost: R8!(11),
                    lower_bound: Some(R8!(0)),
                    upper_bound: Some(R8!(6)),
                    shift: R8!(0),
                    flipped: false
                },
                Variable {
                    variable_type: VariableType::Continuous,
                    cost: R8!(4),
                    lower_bound: Some(R8!(1)),
                    upper_bound: Some(R8!(2)),
                    shift: R8!(0),
                    flipped: false
                },
            ],
            vec!["x".to_string(), "y".to_string()],
            R8!(16),
        );

        let scaling = general_form.scale();

        let expected = GeneralForm::new(
            Objective::Minimize,
            ColumnMajor::from_test_data::<_, _, u8>(&[
                vec![1, 1],
                vec![4, 3],
                vec![1, 1],
            ], 2),
            vec![
                RangedConstraintRelation::Equal,
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Greater,
            ],
            DenseVector::from_test_data::<u8>(vec![3, 6, 3]),
            vec![
                Variable {
                    variable_type: VariableType::Continuous,
                    cost: R8!(11),
                    lower_bound: Some(R8!(0)),
                    upper_bound: Some(R8!(6)),
                    shift: R8!(0),
                    flipped: false,
                },
                Variable {
                    variable_type: VariableType::Continuous,
                    cost: R8!(2),
                    lower_bound: Some(R8!(1, 2)),
                    upper_bound: Some(R8!(1)),
                    shift: R8!(0),
                    flipped: false,
                },
            ],
            vec!["x".to_string(), "y".to_string()],
            R8!(16),
        );

        assert_eq!(general_form, expected);
        assert_eq!(scaling, Scaling {
            c_factor: R8!(0),
            constraint_row_factors: vec![R8!(1), R8!(1, 2)],
            b_factor: R8!(0),
            constraint_column_factors: vec![R8!(1), R8!(1, 2), R8!(1, 7)]
        });
    }
}
