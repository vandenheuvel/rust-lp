use std::num::NonZeroI8;
use std::ops::{Sub};

use crate::algorithm::utilities::merge_sparse_indices;
use crate::data::number_types::nonzero::{NonzeroSign, NonzeroSigned, Nonzero};
use crate::data::number_types::rational::Rational;
use crate::data::number_types::traits::factorization::{NonzeroFactorizable, NonzeroFactorization};

impl<R, Factor: Nonzero + Ord, Power: Sub<Output=Power> + Eq + Nonzero> NonzeroFactorizable for R
where
    R: Rational<
        Numerator: NonzeroFactorizable<Factor=Factor, Power=Power>,
        Denominator: NonzeroFactorizable<Factor=Factor, Power=Power>,
    > + NonzeroSigned,
{
    type Factor = Factor;
    type Power = Power;

    fn factorize(&self) -> NonzeroFactorization<Self::Factor, Self::Power> {
        let NonzeroFactorization { sign: n_sign, factors: n_factors } = self.numerator().factorize();
        let NonzeroFactorization { sign: d_sign, factors: d_factors } = self.denominator().factorize();

        let sign = n_sign * d_sign;
        let factors = merge_sparse_indices(n_factors.into_iter(), d_factors.into_iter(), Sub::sub);

        NonzeroFactorization { sign, factors }
    }
}
