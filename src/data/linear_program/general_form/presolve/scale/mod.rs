use num_traits::One;

use crate::data::linear_algebra::traits::{SparseComparator, SparseElement};
use crate::data::linear_program::general_form::GeneralForm;
use crate::data::linear_algebra::vector::{SparseVector, Vector};
use std::ops::DivAssign;

mod rational;

pub trait Scalable<T> {
    // TODO(ARCHITECTURE): Make the return type more flexible.
    fn scale(&mut self) -> Scaling<T>;
    fn scale_back(&mut self, scale_info: Scaling<T>);
}

impl<T: One + SparseElement<T> + SparseComparator + Clone> Scalable<T> for GeneralForm<T> {
    #[must_use]
    default fn scale(&mut self) -> Scaling<T> {
        // TODO(CORRECTNESS): Are these the right sizes?
        Scaling {
            c_factor: T::one(),
            constraint_row_factors: vec![T::one(); self.nr_active_constraints()],
            b_factor: T::one(),
            constraint_column_factors: vec![T::one(); self.nr_active_variables()]
        }
        // TODO(LOGGING): Log that no scaling is done.
        // TODO(FLOAT): Write the scaling implementation here.
    }

    default fn scale_back(&mut self, _scale_info: Scaling<T>) {
        // TODO(LOGGING): Log that no scaling is done.
        // TODO(FLOAT): Write the scaling implementation here.
    }
}

#[derive(Eq, PartialEq, Debug)]
pub struct Scaling<R> {
    c_factor: R,
    constraint_row_factors: Vec<R>,
    b_factor: R,
    constraint_column_factors: Vec<R>,
}

impl<F> Scaling<F> {
    pub fn scale_back<S>(&self, vector: &mut SparseVector<S, S>)
    where
        for<'r> S: DivAssign<&'r F>,
        S: SparseElement<S> + SparseComparator,
    {
        for (j, value) in vector.iter_values_mut() {
            *value /= &self.constraint_column_factors[*j];
        }
    }
}
