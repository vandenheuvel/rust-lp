use std::cmp::max;
use std::collections::HashSet;
use std::hash::Hash;
use std::ops::{Add, AddAssign, Div, DivAssign, MulAssign, Sub, Neg};

use num_traits::{One, Zero};
use num_traits::Pow;
use relp_num::{NonZeroFactorizable, NonZeroFactorization, Sign, Signed, Abs};

use crate::data::linear_algebra::SparseTuple;
use crate::data::linear_algebra::traits::{SparseComparator, SparseElement};
use crate::data::linear_program::general_form::GeneralForm;
use crate::data::linear_program::general_form::presolve::scale::{RowColumnScaling, Scalable};

impl<R> Scalable<R> for GeneralForm<R>
where
    R: NonZeroFactorizable + DivAssign<R::Factor> + MulAssign<R::Factor> + MulAssign<R> + One + SparseElement<R> + SparseComparator,
    R::Factor: Hash + Pow<R::Power, Output=R::Factor>,
    R::Power: Neg<Output=R::Power> + Sub<Output=R::Power> + Add<Output=R::Power> + Zero + One + AddAssign + Ord,
    for<'r> &'r R: Div<&'r R, Output=R>,
    for<'r> R: DivAssign<&'r R::Factor> + MulAssign<&'r R::Factor>,
{
    fn scale(&mut self) -> RowColumnScaling<R> {
        // Factorize all numbers at once
        let factorization = self.factorize();
        // Optimize per factor
        let scale_per_factor = factorization.solve();
        // Aggregate scaling per factor
        let (row_scales, column_scales) = total_scaling(scale_per_factor, self.nr_constraints(), self.nr_active_variables());

        // Apply the scaling
        for (j, column) in self.constraints.data.iter_mut().enumerate() {
            for (i, value) in column {
                // The row_scales value is large, the column_scales value is small, so combine them
                // first, some factors might cancel
                let total = &row_scales[*i] / &column_scales[j];
                *value *= total;
            }
        }

        RowColumnScaling { row_scales, column_scales }
    }

    fn scale_back(&mut self, scale_info: RowColumnScaling<R>) {
        for (j, column) in self.constraints.data.iter_mut().enumerate() {
            for (i, value) in column {
                // The row_scales value is large, the column_scales value is small
                let total = &scale_info.column_scales[j] / &scale_info.row_scales[*i];
                *value *= total;
            }
        }
    }
}

#[derive(Eq, PartialEq, Clone, Debug)]
struct PerFactor<P> {
    b: Vec<Option<P>>,
    c: Vec<Option<P>>,
    // Column major
    A: Vec<Vec<SparseTuple<P>>>,
    column_info: Vec<ColumnInfo<P>>,
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
    fn solve(mut self) -> Vec<(R::Factor, (Vec<R::Power>, Vec<R::Power>))> {
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
    fn solve_single(&mut self, index: &Vec<Vec<(usize, usize)>>) -> (R::Factor, (Vec<R::Power>, Vec<R::Power>)) {
        let (column_info, row_gains) = self.solve_single_setup(index);
        let row_increments = self.solve_single_rows(row_gains, column_info, index);
        let column_changes = self.solve_single_columns(&row_increments);
        let factor = self.remove_factor_info();

        (factor, (row_increments, column_changes))
    }
    fn solve_single_setup(
        &self, index: &Vec<Vec<(usize, usize)>>,
    ) -> (Vec<ColumnInfo<R::Power>>, Vec<(usize, usize)>) {
        let factor = self.all_factors.last().unwrap();
        let nr_columns = self.A.len();
        let mut infos = Vec::with_capacity(nr_columns);
        for column in &self.A {
            let get_power = |factorization: &[(_, R::Power)]| match factorization.last() {
                Some((f, power)) if f == factor => *power,
                _ => R::Power::zero(),
            };

            let first_power = get_power(&column[0].1);
            let mut lowest_value = first_power;
            let mut lowest_count = 1;
            let mut highest_value = first_power;

            for (_row_index, factorization) in &column[1..] {
                let power = get_power(factorization);

                if power < lowest_value {
                    lowest_value = power;
                    lowest_count = 0;
                }
                if power == lowest_value {
                    lowest_count += 1;
                }
                highest_value = max(power, highest_value);
            }

            infos.push(ColumnInfo { lowest_value, lowest_count, highest_value });
        }

        let row_gains = index.iter().enumerate().map(|(row_index, cols)| {
            let (lower, upper): (Vec<_>, Vec<_>) = cols.iter().map(|&(col_index, data_index)| {
                let value = match self.A[row_index][data_index].1.last() {
                    Some((f, power)) if f == factor => *power,
                    _ => R::Power::zero(),
                };
                let info = &infos[col_index];

                let lower_relevance = value == info.lowest_value && info.lowest_count == 1;
                let upper_relevance = value == info.highest_value;

                (lower_relevance, upper_relevance)
            }).unzip();
            let lower_score = lower.into_iter().filter(|&v| v).count();
            let upper_score = upper.into_iter().filter(|&v| v).count();

            // TODO: Adapt with b and c

            (lower_score, upper_score)
        }).collect::<Vec<_>>();

        (infos, row_gains)
    }

    fn solve_single_rows(
        &self,
        mut row_gains: Vec<(usize, usize)>,
        mut column_info: Vec<ColumnInfo<R::Power>>,
        index: &Vec<Vec<(usize, usize)>>,
    ) -> Vec<R::Power> {
        let factor = self.all_factors.last().unwrap();
        let nr_rows = index.len();
        let mut row_increments = vec![R::Power::zero(); nr_rows];
        let mut not_incremented = (0..nr_rows).collect::<HashSet<_>>();
        while let Some(row_index) = {
            let gain = |i: usize| row_gains[i].0 - row_gains[i].1;
            let candidate = (0..nr_rows).filter(|&index| gain(index) >= 0)
                .max_by_key(|index| {
                    (gain(*index), if not_incremented.contains(index) { 0 } else { 1 })
                });
            candidate.map(|index| {
                if not_incremented.len() == 1 && not_incremented.contains(&index) {
                    None
                } else {
                    Some(index)
                }
            }).flatten()
        } {
            not_incremented.remove(&row_index);

            let gains = &mut row_gains[row_index];

            for &(column, data_index) in &index[row_index] {
                let mut info = &mut column_info[column];
                let power = match self.A[column][data_index].1.last() {
                    Some((f, power)) if f == factor => *power,
                    _ => R::Power::zero(),
                };
                let value = power + row_increments[row_index];

                if value == info.lowest_value {
                    info.lowest_count -= 1;
                    if info.lowest_count == 0 {
                        gains.0 -= 1;

                        // A new lost value, scan for the count
                        info.lowest_value += R::Power::one();
                        for (_, value) in &self.A[column] {
                            let power = match value.last() {
                                Some((f, power)) if f == factor => *power,
                                _ => R::Power::zero(),
                            };

                            if power == info.lowest_value {
                                info.lowest_count += 1;
                            }
                        }
                    }
                }
                if value == info.highest_value - R::Power::one() {
                    gains.1 += 1;
                }
                if value == info.highest_value {
                    info.highest_value += R::Power::one();
                    // No need to adapt gains
                }
            }

            debug_assert_eq!(gains.0, 0);

            row_increments[row_index] += R::Power::one();
        }

        row_increments
    }

    fn solve_single_columns(&self, row_increments: &Vec<R::Power>) -> Vec<R::Power> {
        let factor = self.all_factors.last().unwrap();
        self.A.iter().map(|column| {
            let mut powers = column.iter().map(|(row_index, factorization)| {
                row_increments[*row_index] + match factorization.last() {
                    Some((f, power)) if f == factor => *power,
                    _ => R::Power::zero(),
                }
            }).collect::<Vec<_>>();
            powers.sort_unstable();
            let middle = powers.len() / 2;
            debug_assert!(!powers.is_empty(), "At least one element per column (even if zero).");
            powers[middle]
        }).collect()
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
    scale_per_factor: Vec<(R::Factor, (Vec<R::Power>, Vec<R::Power>))>,
    nr_rows: usize, nr_columns: usize,
) -> (Vec<R>, Vec<R>)
where
    R: NonZeroFactorizable<Power: Abs>,
    for<'r> R: MulAssign<&'r R::Factor> + DivAssign<&'r R::Factor>,
{
    debug_assert!(!scale_per_factor.is_empty());

    let mut row_scales = vec![R::one(); nr_rows];
    let mut column_scales = vec![R::one(); nr_columns];

    for (factor, (row_changes, column_changes)) in scale_per_factor {
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

    (row_scales, column_scales)
}

#[cfg(test)]
mod test {
    use relp_num::{R8, Rational8};

    use crate::data::linear_algebra::matrix::{ColumnMajor, Order};
    use crate::data::linear_algebra::vector::{DenseVector, Vector};
    use crate::data::linear_algebra::vector::test::TestVector;
    use crate::data::linear_program::elements::{Objective, RangedConstraintRelation, VariableType};
    use crate::data::linear_program::general_form::{GeneralForm, Variable};
    use crate::data::linear_program::general_form::presolve::scale::{RowColumnScaling, Scalable};
    use crate::data::linear_program::general_form::presolve::scale::rational::Factorization;

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
            DenseVector::from_test_data::<u8>(vec![3, 12, 21, 0]),
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

        let factorization = general_form.factorize();

        let expected = Factorization {
            all_factors: vec![2, 3, 7, 11].into_iter().collect(),
            b: vec![Some(vec![(3, 1)]), Some(vec![(2, 2), (3, 1)]), Some(vec![(3, 1), (7, 1)]), None],
            c: vec![Some(vec![(11, 1)]), Some(vec![(2, 2)])],
            A: vec![
                vec![(0, vec![]), (1, vec![(2, 2)]), (2, vec![(7, 1)])],
                vec![(0, vec![(2, 1)]), (1, vec![(2, 1), (3, 1)]), (2, vec![(2, 1), (7, 1)]), (3, vec![(11, 1)])],
            ],
        };

        assert_eq!(factorization, expected);
    }

    #[test]
    fn test_solve() {

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
        assert_eq!(scaling, RowColumnScaling {
            row_scales: vec![R8!(1), R8!(1, 2)],
            column_scales: vec![R8!(1), R8!(1, 2), R8!(1, 7)],
        });
    }
}
