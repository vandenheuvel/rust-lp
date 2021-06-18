use num_traits::One;

use crate::data::linear_algebra::traits::{SparseComparator, SparseElement};
use crate::data::linear_program::general_form::GeneralForm;

mod rational;

pub trait Scalable<T> {
    // TODO(ARCHITECTURE): Make the return type more flexible.
    fn scale(&mut self) -> RowColumnScaling<T>;
    fn scale_back(&mut self, scale_info: RowColumnScaling<T>);
}

#[derive(Eq, PartialEq, Debug)]
pub struct RowColumnScaling<R> {
    row_scales: Vec<R>,
    column_scales: Vec<R>,
}

impl<T: One + SparseElement<T> + SparseComparator + Clone> Scalable<T> for GeneralForm<T> {
    default fn scale(&mut self) -> RowColumnScaling<T> {
        // TODO(CORRECTNESS): Are these the right sizes?
        RowColumnScaling {
            row_scales: vec![T::one(); self.nr_active_constraints()],
            column_scales: vec![T::one(); self.nr_active_variables()],
        }
        // TODO(LOGGING): Log that no scaling is done.
        // TODO(FLOAT): Write the scaling implementation here.
    }

    default fn scale_back(&mut self, _scale_info: RowColumnScaling<T>) {
        // TODO(LOGGING): Log that no scaling is done.
        // TODO(FLOAT): Write the scaling implementation here.
    }
}
