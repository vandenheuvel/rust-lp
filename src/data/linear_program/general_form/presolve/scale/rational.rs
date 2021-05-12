use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};
use std::hash::Hash;
use std::ops::{Add, AddAssign, Div, DivAssign, MulAssign, Sub};

use itertools::izip;
use num::{One, Zero};
use num::traits::Pow;

use crate::data::linear_algebra::SparseTuple;
use crate::data::linear_algebra::traits::{SparseComparator, SparseElement};
use crate::data::linear_program::general_form::GeneralForm;
use crate::data::linear_program::general_form::presolve::scale::{Scalable, ScaleInfo};
use crate::data::number_types::traits::factorization::{NonzeroFactorizable, NonzeroFactorization};

impl<R> Scalable<R> for GeneralForm<R>
where
    R: NonzeroFactorizable + DivAssign<R::Factor> + MulAssign<R::Factor> + MulAssign<R> + One + SparseElement<R> + SparseComparator,
    R::Factor: Hash + Pow<R::Power, Output=R::Factor>,
    R::Power: Sub<Output=R::Power> + Add<Output=R::Power> + Zero + One + AddAssign + Ord,
    for<'r> &'r R: Div<&'r R, Output=R>,
{
    fn scale(&mut self) -> ScaleInfo<R> {
        let factorization = self.factorize();
        let mut row_major_column_indices = vec![Vec::new(); self.nr_constraints()];
        for (j, column) in self.constraints.iter_columns().enumerate() {
            for (data_index, &(i, _)) in column.iter().enumerate() {
                row_major_column_indices[i].push((j, data_index));
            }
        }
        let factorization_per_factor = factorization.split_into_factors();
        let solutions = factorization_per_factor.into_iter().map(|(factor, pf)| {
            (factor, pf.solve(&row_major_column_indices))
        }).collect::<Vec<_>>();

        let mut row_scales = vec![R::one(); self.nr_constraints()];
        let mut column_scales = vec![R::one(); self.nr_variables()];

        for (factor, (row_changes, column_changes)) in solutions {
            for (i, row_change) in row_changes.into_iter().enumerate() {
                match row_change.cmp(&R::Power::zero()) {
                    Ordering::Less => {row_scales[i] /= factor.clone().pow(row_change)},
                    Ordering::Equal => {}
                    Ordering::Greater => {row_scales[i] *= factor.clone().pow(row_change)},
                }
            }
            for (j, column_change) in column_changes.into_iter().enumerate() {
                match column_change.cmp(&R::Power::zero()) {
                    Ordering::Less => {column_scales[j] /= factor.clone().pow(column_change);}
                    Ordering::Equal => {}
                    Ordering::Greater => {column_scales[j] *= factor.clone().pow(column_change);}
                }
            }
        }

        for (j, column) in self.constraints.data.iter_mut().enumerate() {
            for (i, value) in column {
                // The row_scales value is large, the column_scales value is small, so combine them
                // first, some factors might cancel
                let total = &row_scales[*i] / &column_scales[j];
                *value *= total;
            }
        }

        ScaleInfo { row_scales, column_scales }
    }

    fn scale_back(&mut self, scale_info: ScaleInfo<R>) {
        for (j, column) in self.constraints.data.iter_mut().enumerate() {
            for (i, value) in column {
                // The row_scales value is large, the column_scales value is small
                let total = &scale_info.column_scales[j] / &scale_info.row_scales[*i];
                *value *= total;
            }
        }
    }
}

#[derive(Clone, Debug)]
struct PerFactor<P> {
    b: Vec<P>,
    c: Vec<P>,
    // Column major
    A: Vec<Vec<SparseTuple<P>>>,
    column_info: Vec<ColumnInfo<P>>,
}

#[derive(Clone, Debug)]
struct ColumnInfo<P> {
    /// Signed, needs at most 32 bits for arbitrary precision numbers
    lowest_value: P,
    /// Relates to the number of rows
    lowest_count: u32,
    /// Signed, needs at most 32 bits for arbitrary precision numbers
    highest_value: P,
}

impl<T> PerFactor<T>
where
    T: Sub<Output=T> + Add<Output=T> + AddAssign + PartialEq + Zero + One + Ord + Copy
{
    fn solve(mut self, index: &Vec<Vec<(usize, usize)>>) -> (Vec<T>, Vec<T>) {
        let nr_rows = self.b.len();

        let mut row_gains = index.iter().enumerate().map(|(row_index, cols)| {
            let (lower, upper): (Vec<_>, Vec<_>) = cols.iter().map(|&(col_index, data_index)| {
                let value = &self.A[row_index][data_index].1;
                let info = &self.column_info[col_index];

                let lower_relevance = value == &info.lowest_value && info.lowest_count == 1;
                let upper_relevance = value == &info.highest_value;

                (lower_relevance, upper_relevance)
            }).unzip();
            let lower_score = lower.into_iter().filter(|&v| v).count();
            let upper_score = upper.into_iter().filter(|&v| v).count();

            // TODO: Adapt with b and c

            (lower_score, upper_score)
        }).collect::<Vec<_>>();

        let mut row_increments = vec![T::zero(); nr_rows];
        let mut not_incremented = (0..nr_rows).collect::<HashSet<_>>();

        while let Some(row_index) = {
            let change = |i: usize| row_gains[i].0 - row_gains[i].1;
            let candidate = (0..nr_rows).filter(|&index| change(index) >= 0)
                .max_by_key(|index| {
                    (change(*index), if not_incremented.contains(index) { 0 } else { 1 })
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
                let mut info = &mut self.column_info[column];
                let value = &self.A[column][data_index].1;

                if value == &info.lowest_value {
                    info.lowest_count -= 1;
                    if info.lowest_count == 0 {
                        gains.0 -= 1;

                        // A new lost value, scan for the count
                        info.lowest_value += T::one();
                        for (_, value) in &self.A[column] {
                            if value == &info.lowest_value {
                                info.lowest_count += 1;
                            }
                        }
                    }
                }
                if *value == info.highest_value - T::one() {
                    gains.1 += 1;
                }
                if value == &info.highest_value {
                    info.highest_value += T::one();
                    // No need to adapt gains
                }

                self.A[column][data_index].1 += T::one();
            }

            debug_assert_eq!(gains.0, 0);

            row_increments[row_index] += T::one();
        }

        let column_increments = self.A.into_iter().map(|column| {
            let mut powers = column.into_iter().map(|(_, p)| p).collect::<Vec<_>>();
            powers.sort_unstable();
            let middle = powers.len() / 2;
            debug_assert!(!powers.is_empty(), "At least one element per column (even if zero).");
            powers[middle]
        }).collect();

        (row_increments, column_increments)
    }
}

struct Factorization<R: NonzeroFactorizable> {
    /// No duplicate elements
    all_factors: Vec<R::Factor>,
    /// Zero values are None, others have a factorization that might be empty
    b: Vec<Option<Vec<(R::Factor, R::Power)>>>,
    c: Vec<Option<Vec<(R::Factor, R::Power)>>>,
    /// Column major
    A: Vec<Vec<SparseTuple<Vec<(R::Factor, R::Power)>>>>,
}

impl<R: NonzeroFactorizable<Factor: Hash, Power: Zero + Ord>> Factorization<R> {
    fn split_into_factors(self) -> Vec<(R::Factor, PerFactor<R::Power>)> {
        let number_of_factors = self.all_factors.len();
        let nr_rows = self.b.len();

        let mut A_template = Vec::with_capacity(self.A.len());
        for column in &self.A {
            let mut new_column = Vec::with_capacity(column.len());
            for &(row_index, _) in column {
                new_column.push((row_index, R::Power::zero()));
            }
            A_template.push(new_column);
        }
        let mut As = vec![A_template; number_of_factors];

        let mapping = self.all_factors.iter().enumerate()
            .map(|(index, factor)| (factor.clone(), index))
            .collect::<HashMap<R::Factor, usize>>();

        for (column_index, column) in self.A.into_iter().enumerate() {
            for (row_index, factors) in column {
                for (factor, power) in factors {
                    let column = &mut As[*mapping.get(&factor).unwrap()][column_index];
                    let data_index = column.binary_search_by_key(&row_index, |&(i, _)| i).into_ok_or_err();
                    column[data_index].1 = power;
                }
            }
        }
        let mut column_infos = Vec::with_capacity(number_of_factors);
        for factor_index in 0..number_of_factors {
            let mut infos_for_factor = Vec::new();
            for column in &As[factor_index] {
                let mut info = ColumnInfo {
                    lowest_value: column[0].1, // Some initial value
                    // When the first value is encountered again in the loop this becomes 1
                    lowest_count: 0,
                    highest_value: column[0].1, // Some initial value
                };

                for &(_, power) in column {
                    match power.cmp(&info.lowest_value) {
                        Ordering::Less => {
                            info.lowest_value = power;
                            info.lowest_count = 1;
                        }
                        Ordering::Equal => {
                            info.lowest_count += 1;
                        }
                        Ordering::Greater => {}
                    }

                    if power > info.highest_value {
                        info.highest_value = power;
                    }
                }

                infos_for_factor.push(info);
            }
            column_infos.push(infos_for_factor);
        }

        let mut bs = vec![vec![R::Power::zero(); nr_rows]; number_of_factors];
        for (row, value) in self.b.into_iter().enumerate() {
            if let Some(factors) = value {
                for (factor, power) in factors {
                    bs[*mapping.get(&factor).unwrap()][row] = power;
                }
            }
        }

        let mut cs = vec![vec![R::Power::zero(); nr_rows]; number_of_factors];
        for (column, value) in self.c.into_iter().enumerate() {
            if let Some(factors) = value {
                for (factor, power) in factors {
                    cs[*mapping.get(&factor).unwrap()][column] = power;
                }
            }
        }

        izip!(self.all_factors, bs, cs, As, column_infos).map(|(factor, b, c, A, column_info)| {
            (factor, PerFactor { b, c, A, column_info })
        }).collect()
    }
}

impl<R> GeneralForm<R>
where
    R: NonzeroFactorizable<Factor: Hash> + SparseElement<R> + SparseComparator,
{
    fn factorize(&self) -> Factorization<R> {
        let mut all_factors = HashSet::new();
        let mut all_factors_vec = Vec::new();

        let b = self.b.data.iter()
            .map(|v| {
                if v.is_not_zero() {
                    let NonzeroFactorization { factors, .. } = v.factorize();
                    for (factor, _count) in &factors {
                        if all_factors.insert(factor.clone()) {
                            all_factors_vec.push(factor.clone());
                        }
                    }
                    Some(factors)
                } else {
                    None
                }
            })
            .collect();

        let c = self.variables.iter()
            .map(|variable| {
                let NonzeroFactorization { factors, .. } = variable.cost.factorize();
                for (factor, _count) in &factors {
                    if all_factors.insert(factor.clone()) {
                        all_factors_vec.push(factor.clone());
                    }
                }
                Some(factors)
            })
            .collect();

        let A = self.constraints.iter_columns().map(|column| {
            column.iter().map(|(i, v)| {
                let NonzeroFactorization { factors, .. } = v.factorize();
                for (factor, _count) in &factors {
                    if all_factors.insert(factor.clone()) {
                        all_factors_vec.push(factor.clone());
                    }
                }
                (*i, factors)
            }).collect()
        }).collect();

        Factorization {
            all_factors: all_factors_vec,
            b,
            c,
            A,
        }
    }
}

#[cfg(test)]
mod test {
    use crate::data::linear_algebra::matrix::{ColumnMajor, Order, Sparse as SparseMatrix};
    use crate::data::linear_algebra::vector::{DenseVector, Vector};
    use crate::data::linear_program::elements::{Objective, RangedConstraintRelation, VariableType};
    use crate::data::linear_program::general_form::{GeneralForm, Variable};

    #[test]
    fn test_factorize() {
        // let general_form = GeneralForm::new(
        //     Objective::Minimize,
        //     ColumnMajor::from_test_data(&[
        //         vec![]
        //     ], 2),
        //     vec![RangedConstraintRelation::Equal],
        //     DenseVector::new(vec![]),
        //     vec![Variable {
        //         variable_type: VariableType::Continuous,
        //         cost: (),
        //         lower_bound: None,
        //         upper_bound: None,
        //         shift: (),
        //         flipped: false
        //     }],
        //     vec!["x".to_string()],
        //     0,
        // );

        // TODO
    }
}
