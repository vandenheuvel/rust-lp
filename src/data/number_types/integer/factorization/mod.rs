//! # Factorizing integers
use gcd::Gcd;
use num::Integer;

use crate::data::number_types::integer::factorization::prime::Prime;
use crate::data::number_types::integer::factorization::prime::primes::SMALL_ODD_PRIMES;
use crate::data::number_types::nonzero::NonzeroSigned;
use crate::data::number_types::nonzero::sign::Sign;
use crate::data::number_types::traits::factorization::{NonzeroFactorizable, NonzeroFactorization};

pub mod prime;

impl NonzeroFactorizable for i32 {
    type Factor = u32;
    // TODO(CORRECTNESS): Consider making this type unsigned
    type Power = i8;

    fn factorize(&self) -> NonzeroFactorization<Self::Factor, Self::Power> {
        let sign = NonzeroSigned::signum(self);
        let NonzeroFactorization {
            factors, ..
        } = (self.abs() as u32).factorize();
        NonzeroFactorization { sign, factors }
    }
}

impl NonzeroFactorizable for u32 {
    type Factor = u32;
    // TODO(CORRECTNESS): Consider making this type unsigned
    type Power = i8;

    fn factorize(&self) -> NonzeroFactorization<Self::Factor, Self::Power> {
        debug_assert_ne!(*self, 0);

        let mut x = self.clone();
        let mut unsorted_factors = Vec::with_capacity(32);

        while x.is_even() {
            x /= 2;
            unsorted_factors.push(2);
        }

        'odd_trial_division: {
            // smallest
            for divisor in SMALL_ODD_PRIMES {
                while x % divisor as u32 == 0 {
                    x /= divisor as u32;
                    unsorted_factors.push(divisor as u32);
                }

                if x == 1 {
                    break 'odd_trial_division;
                }
            }
            // small
            let first = *SMALL_ODD_PRIMES.last().unwrap() as u32 + 2;
            let last = 2_u32.pow(32 / 2);
            let mut sqrt = ((x as f64).sqrt() + 2_f64) as u32;
            for divisor in (first..last).step_by(2) {
                while x % divisor == 0 {
                    x /= divisor;
                    sqrt = ((x as f64).sqrt() + 2_f64) as u32;
                    unsorted_factors.push(divisor);
                }

                if x == 1 {
                    break 'odd_trial_division;
                } else if divisor > sqrt {
                    unsorted_factors.push(x);
                    break 'odd_trial_division;
                }
            }
        }

        // Aggregate the factors
        let mut factors = Vec::with_capacity(16);
        for factor in unsorted_factors {
            if let Some((existing_factor, count)) = factors.last_mut() {
                if *existing_factor == factor {
                    *count += 1;
                    continue;
                }
            }

            factors.push((factor, 1));
        }

        NonzeroFactorization { sign: Sign::Positive, factors }
    }
}

impl NonzeroFactorizable for i64 {
    type Factor = u64;
    // TODO(CORRECTNESS): Consider making this type unsigned
    type Power = i8;

    fn factorize(&self) -> NonzeroFactorization<Self::Factor, Self::Power> {
        let sign = NonzeroSigned::signum(self);
        let NonzeroFactorization {
            factors, ..
        } = (self.abs() as u64).factorize();
        NonzeroFactorization { sign, factors }
    }
}

impl NonzeroFactorizable for u64 {
    type Factor = u64;
    // TODO(CORRECTNESS): Consider making this type unsigned
    type Power = i8;

    fn factorize(&self) -> NonzeroFactorization<Self::Factor, Self::Power> {
        debug_assert_ne!(*self, 0);

        let mut x = self.clone();
        let mut unsorted_factors = Vec::with_capacity(64);

        // Trial division
        // 2
        while x.is_even() {
            x /= 2;
            unsorted_factors.push(2);
        }
        'odd: {
            // smallest
            for divisor in SMALL_ODD_PRIMES {
                while x % divisor as u64 == 0 {
                    x /= divisor as u64;
                    unsorted_factors.push(divisor as u64);
                }

                if x == 1 {
                    break 'odd;
                }
            }
            // small
            let first = *SMALL_ODD_PRIMES.last().unwrap() as u64 + 2;
            let last = 500_000;
            // TODO(PERFORMANCE): How far should this loop go?
            let mut sqrt = ((x as f64).sqrt() + 2_f64) as u64;
            for divisor in (first..last).step_by(2) {
                while x % divisor == 0 {
                    x /= divisor;
                    sqrt = ((x as f64).sqrt() + 2_f64) as u64;
                    unsorted_factors.push(divisor);
                }

                if x == 1 {
                    break 'odd;
                } else if divisor > sqrt {
                    // Checked enough, must be prime
                    unsorted_factors.push(x);
                    break 'odd;
                }
            }

            // Prime test and Pollard's rho
            rho_loop(x, &mut unsorted_factors);
        }

        // Sort and aggregate the factors
        unsorted_factors.sort();
        let mut factors = Vec::with_capacity(16);
        for factor in unsorted_factors {
            if let Some((existing_factor, count)) = factors.last_mut() {
                if *existing_factor == factor {
                    *count += 1;
                    continue;
                }
            }

            factors.push((factor, 1));
        }

        NonzeroFactorization { sign: Sign::Positive, factors }
    }
}

// TODO(PERFORMANCE): How long should this be tried?
const RHO_BASE_LIMIT: u64 = 20;

fn rho_loop(mut x: u64, factors: &mut Vec<u64>) {
    debug_assert_ne!(x, 0);

    let mut e = 2;
    while x > 1 {
        if x.is_prime() || e == RHO_BASE_LIMIT {
            factors.push(x);
            break;
        }

        match rho(x, e) {
            None | Some(1) => {
                // TODO(PERFORMANCE): Odd values only?
                e += 1;
            }
            Some(factor) => {
                if factor.is_prime() {
                    factors.push(factor);
                } else {
                    // TODO(PERFORMANCE): Should the `e` variable be passed to the inner call?
                    rho_loop(factor, factors);
                }
                x /= factor;
            }
        }
    }
}

/// Pollard's rho function generates a divisor.
///
/// Up to minor adaptions, this code is from the reikna repository developed by Phillip Heikoop.
pub fn rho(value: u64, entropy: u64) -> Option<u64> {
    debug_assert_ne!(value, 0);
    debug_assert_ne!(value, 1);
    debug_assert_ne!(value, 2);

    let entropy = entropy.wrapping_mul(value);
    let c = entropy & 0xff;
    let u = entropy & 0x7f;

    let mut r: u64 = 1;
    let mut q: u64 = 1;
    let mut y: u64 = entropy & 0xf;

    let mut factor = 1;

    let mut y_old = 0;
    let mut x = 0;

    let f = |x: u64| (x.wrapping_mul(x) + c) % value;

    while factor == 1 {
        x = y;

        for _ in 0..r {
            y = f(y);
        }

        let mut k = 0;
        while k < r && factor == 1 {
            y_old = y;

            for _ in 0..u64::min(u, r - k) {
                y = f(y);

                if x > y {
                    q = q.wrapping_mul(x - y) % value;
                } else {
                    q = q.wrapping_mul(y - x) % value;
                }
            }

            factor = Gcd::gcd(q, value);
            k += u;
        }

        r *= 2;
    }


    while factor == value || factor <= 1 {
        y_old = f(y_old);

        if x > y_old {
            factor = Gcd::gcd(x - y_old, value);
        } else if x < y_old {
            factor = Gcd::gcd(y_old - x, value);
        } else {
            // the algorithm has failed for this entropy,
            // return the factor as-is
            return None;
        }
    }

    Some(factor)
}

#[cfg(test)]
mod test {
    use crate::data::number_types::nonzero::sign::Sign;
    use crate::data::number_types::traits::factorization::{NonzeroFactorizable, NonzeroFactorization};

    #[test]
    fn test_factorize_one() {
        assert_eq!(1_u64.factorize(), NonzeroFactorization { sign: Sign::Positive, factors: vec![]});
    }

    #[test]
    fn test_factorize_two_powers() {
        assert_eq!(
            2_u64.factorize(),
            NonzeroFactorization { sign: Sign::Positive, factors: vec![(2, 1)] },
        );
        assert_eq!(
            4_u64.factorize(),
            NonzeroFactorization { sign: Sign::Positive, factors: vec![(2, 2)] },
        );
        assert_eq!(
            2_u64.pow(16).factorize(),
            NonzeroFactorization { sign: Sign::Positive, factors: vec![(2, 16)] },
        );
    }

    #[test]
    fn test_factorize_three_powers() {
        assert_eq!(
            3_u64.factorize(),
            NonzeroFactorization { sign: Sign::Positive, factors: vec![(3, 1)] },
        );
        assert_eq!(
            9_u64.factorize(),
            NonzeroFactorization { sign: Sign::Positive, factors: vec![(3, 2)] },
        );
        assert_eq!(
            3_u64.pow(16).factorize(),
            NonzeroFactorization { sign: Sign::Positive, factors: vec![(3, 16)] },
        );
    }

    #[test]
    fn test_factorize_mixed() {
        assert_eq!(
            6_u64.factorize(),
            NonzeroFactorization { sign: Sign::Positive, factors: vec![(2, 1), (3, 1)] },
        );
        assert_eq!(
            36_u64.factorize(),
            NonzeroFactorization { sign: Sign::Positive, factors: vec![(2, 2), (3, 2)] },
        );
        assert_eq!(
            6_u64.pow(16).factorize(),
            NonzeroFactorization { sign: Sign::Positive, factors: vec![(2, 16), (3, 16)] },
        );
    }

    #[test]
    fn test_factorize_large_prime() {
        let prime = 2_u64.pow(36) - 5;
        assert_eq!(
            prime.factorize(),
            NonzeroFactorization { sign: Sign::Positive, factors: vec![(prime, 1)] },
        );
        let prime = 2_u64.pow(60) - 93;
        assert_eq!(
            prime.factorize(),
            NonzeroFactorization { sign: Sign::Positive, factors: vec![(prime, 1)] },
        );
        let prime = 2_u64.pow(53) - 111;
        assert_eq!(
            prime.factorize(),
            NonzeroFactorization { sign: Sign::Positive, factors: vec![(prime, 1)] },
        );
    }

    #[test]
    fn test_factorize_large_composite_u64() {
        let composite = 2_u64.pow(36) - 7;
        assert_eq!(
            composite.factorize(),
            NonzeroFactorization { sign: Sign::Positive, factors: vec![(3, 1), (22906492243, 1)] },
        );
        let composite = 2_u64.pow(60) - 95;
        assert_eq!(
            composite.factorize(),
            NonzeroFactorization { sign: Sign::Positive, factors: vec![(3457, 1), (6203, 1), (53764867411, 1)] },
        );
        let composite = 2_u64.pow(53) - 113;
        assert_eq!(
            composite.factorize(),
            NonzeroFactorization { sign: Sign::Positive, factors: vec![(3, 2), (43, 1), (642739, 1), (36211303, 1)] },
        );
    }

    #[test]
    fn test_factorize_large_composite_i64() {
        let composite = 2_i64.pow(36) - 7;
        assert_eq!(
            composite.factorize(),
            NonzeroFactorization { sign: Sign::Positive, factors: vec![(3, 1), (22906492243, 1)] },
        );
        let composite = 2_i64.pow(60) - 95;
        assert_eq!(
            composite.factorize(),
            NonzeroFactorization { sign: Sign::Positive, factors: vec![(3457, 1), (6203, 1), (53764867411, 1)] },
        );
        let composite = 2_i64.pow(53) - 113;
        assert_eq!(
            composite.factorize(),
            NonzeroFactorization { sign: Sign::Positive, factors: vec![(3, 2), (43, 1), (642739, 1), (36211303, 1)] },
        );
    }
}
