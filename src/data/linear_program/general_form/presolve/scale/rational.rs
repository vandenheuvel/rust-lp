use crate::data::linear_program::general_form::presolve::scale::Scalable;
use crate::data::linear_program::general_form::GeneralForm;
use crate::data::number_types::traits::factorization::{NonzeroFactorizable, NonzeroFactorization};
use std::collections::{HashSet, HashMap};
use crate::data::linear_algebra::SparseTuple;
use num::Zero;
use crate::data::linear_algebra::traits::{SparseElement, SparseComparator};
use std::hash::Hash;
use std::cmp::Ordering;
use itertools::izip;

impl<R> Scalable for GeneralForm<R>
where
    R: NonzeroFactorizable<Factor: Hash, Power: Zero + Ord> + SparseElement<R> + SparseComparator,
{
    fn scale(&mut self) {
        let factorization = self.factorize();
        let mut row_major_column_indices = vec![Vec::new(); self.nr_constraints()];
        for (j, column) in self.constraints.iter_columns().enumerate() {
            for (data_index, &(i, _)) in column.iter().enumerate() {
                row_major_column_indices[i].push((j, data_index));
            }
        }
        let factorization_per_factor = factorization.split_into_factors();


        unimplemented!();
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
    lowest_value: P,
    lowest_count: i16,
    highest_value: P,
    highest_count: i16,
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
                    lowest_value: column[0].1.clone(),
                    lowest_count: 0,
                    highest_value: column[0].1.clone(),
                    highest_count: 0,
                };

                for (_, power) in column {
                    match power.cmp(&info.lowest_value) {
                        Ordering::Less => {
                            info.lowest_value = power.clone();
                            info.lowest_count = 1;
                        }
                        Ordering::Equal => {
                            info.lowest_count += 1;
                        }
                        Ordering::Greater => {}
                    }

                    match power.cmp(&info.highest_value) {
                        Ordering::Less => {}
                        Ordering::Equal => {
                            info.highest_count += 1;
                        }
                        Ordering::Greater => {
                            info.highest_value = power.clone();
                            info.highest_count = 1;
                        }
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
