use crate::data::linear_program::general_form::GeneralForm;

mod rational;

pub trait Scalable {
    fn scale(&mut self);
}

impl<T> Scalable for GeneralForm<T> {
    default fn scale(&mut self) {
        // TODO(LOGGING): Log that no scaling is done.
        // TODO(FLOAT): Write the scaling implementation here.
    }
}
