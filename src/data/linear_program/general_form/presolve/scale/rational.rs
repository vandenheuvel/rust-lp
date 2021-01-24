use crate::data::linear_program::general_form::presolve::scale::Scalable;
use crate::data::linear_program::general_form::GeneralForm;
use crate::data::number_types::rational::Rational;
use crate::data::number_types::nonzero::Nonzero;
use crate::data::number_types::traits::factorization::NonzeroFactorizable;

impl<R: Rational + NonzeroFactorizable> Scalable for GeneralForm<R> {
    fn scale(&mut self) {
        let factorization = self.factorize();

        unimplemented!();
    }
}

struct Factorization<R: NonzeroFactorizable> {
    /// Zero values are None, others have a factorization that might be empty
    b: Vec<Option<Vec<(R::Factor, R::Power)>>>,
}

impl<R: Rational + NonzeroFactorizable> GeneralForm<R> {
    fn factorize(&self) -> Factorization<R> {
        unimplemented!();
        // let b = self.b.data.iter()
        //     .map(|v| {
        //         if v.is_not_zero() {
        //             Some(v.factorize())
        //         } else {
        //             None
        //         }
        //     })
        //     .collect();
        //
        // Factorization {
        //     b,
        // }
    }
}
