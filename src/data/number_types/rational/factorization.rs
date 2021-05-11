use std::ops::{Sub, Neg};

use crate::algorithm::utilities::merge_sparse_indices;
use crate::data::number_types::nonzero::{Nonzero};
use crate::data::number_types::traits::factorization::{NonzeroFactorizable, NonzeroFactorization};
use crate::data::number_types::rational::Ratio;
use std::convert::identity;

impl<Numerator, Denominator, Factor, Power> NonzeroFactorizable for Ratio<Numerator, Denominator>
where
    Numerator: NonzeroFactorizable<Factor=Factor, Power=Power>,
    Denominator: NonzeroFactorizable<Factor=Factor, Power=Power>,
    Factor: Ord + Nonzero,
    Power: Neg<Output=Power> + Sub<Output=Power> + Eq + Nonzero,
    Ratio<Numerator, Denominator>: Nonzero,
{
    type Factor = Factor;
    type Power = Power;

    fn factorize(&self) -> NonzeroFactorization<Self::Factor, Self::Power> {
        debug_assert!(self.numerator.is_not_zero() && self.denominator.is_not_zero());

        let NonzeroFactorization { sign: n_sign, factors: n_factors } = self.numerator.factorize();
        let NonzeroFactorization { sign: d_sign, factors: d_factors } = self.denominator.factorize();

        let sign = n_sign * d_sign;
        let factors = merge_sparse_indices(
            n_factors.into_iter(), d_factors.into_iter(),
            Sub::sub, identity, |x| -x,
        );

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
        let ratio: Ratio<u64, u64> = Ratio { numerator: 1, denominator: 2 };
        let expected = NonzeroFactorization { sign: Sign::Positive, factors: vec![(2, -1)] };
        assert_eq!(ratio.factorize(), expected);

        let ratio: Ratio<u64, u64> = Ratio { numerator: 161, denominator: 3 };
        let expected = NonzeroFactorization { sign: Sign::Positive, factors: vec![(3, -1), (7, 1), (23, 1)] };
        assert_eq!(ratio.factorize(), expected);

        let ratio: Ratio<u64, u64> = Ratio { numerator: 2, denominator: 1 };
        let expected = NonzeroFactorization { sign: Sign::Positive, factors: vec![(2, 1)] };
        assert_eq!(ratio.factorize(), expected);

        let ratio: Ratio<u64, u64> = Ratio { numerator: 1, denominator: 1 };
        let expected = NonzeroFactorization { sign: Sign::Positive, factors: vec![] };
        assert_eq!(ratio.factorize(), expected);
    }
}
