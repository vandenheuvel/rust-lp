use std::num::NonZeroI8;
use std::ops::{Sub};

use crate::algorithm::utilities::merge_sparse_indices;
use crate::data::number_types::nonzero::{NonzeroSign, NonzeroSigned, Nonzero};
use crate::data::number_types::traits::factorization::{NonzeroFactorizable, NonzeroFactorization};
use crate::data::number_types::rational::Ratio;

impl<Numerator, Denominator, Factor, Power> NonzeroFactorizable for Ratio<Numerator, Denominator>
where
    Numerator: NonzeroFactorizable<Factor=Factor, Power=Power>,
    Denominator: NonzeroFactorizable<Factor=Factor, Power=Power>,
    Factor: Ord + Nonzero,
    Power: Sub<Output=Power> + Eq + Nonzero,
    Ratio<Numerator, Denominator>: Nonzero,
{
    type Factor = Factor;
    type Power = Power;

    fn factorize(&self) -> NonzeroFactorization<Self::Factor, Self::Power> {
        let NonzeroFactorization { sign: n_sign, factors: n_factors } = self.numerator.factorize();
        let NonzeroFactorization { sign: d_sign, factors: d_factors } = self.denominator.factorize();

        let sign = n_sign * d_sign;
        let factors = merge_sparse_indices(n_factors.into_iter(), d_factors.into_iter(), Sub::sub);

        NonzeroFactorization { sign, factors }
    }
}

#[cfg(test)]
mod test {
    use crate::data::number_types::rational::Ratio;
    use crate::data::number_types::traits::factorization::{NonzeroFactorizable, NonzeroFactorization};
    use crate::data::number_types::nonzero::sign::Sign;

    #[test]
    fn test_factorize() {
        let ratio = Ratio { numerator: 1, denominator: 2 };
        let expected = NonzeroFactorization { sign: Sign::Positive, factors: vec![(2, -1)] };
        assert_eq!(ratio.factorize(), expected);
    }
}
