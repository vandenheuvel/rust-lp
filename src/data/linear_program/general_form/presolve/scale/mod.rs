use num_traits::One;

use crate::data::linear_algebra::traits::{SparseComparator, SparseElement};
use crate::data::linear_program::general_form::GeneralForm;

mod rational;

pub trait Scalable<T> {
    // TODO(ARCHITECTURE): Make the return type more flexible.
    fn scale(&mut self) -> Scaling<T>;
    fn scale_back(&mut self, scale_info: Scaling<T>);
}

#[derive(Eq, PartialEq, Debug)]
pub struct Scaling<R> {
    c_factor: R,
    constraint_row_factors: Vec<R>,
    b_factor: R,
    constraint_column_factors: Vec<R>,
}

impl<T: One + SparseElement<T> + SparseComparator + Clone> Scalable<T> for GeneralForm<T> {
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
