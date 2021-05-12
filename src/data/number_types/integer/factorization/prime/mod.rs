//! # Primes
//!
//! uses the method by Michal Forišek  and Jakub Jančina from their 2015 paper "fast primality
//! testing for integers that fit into a machine word". it divides all values a machine word can
//! take over a number of buckets pseudo-randomly through hashing. for each bucket, they provide a
//! precomputed parameter that can be used with a standard check to see whether the number is a
//! probable prime. their computational work ensures that for the values in the bucket and the
//! provided parameter, the probable primes correspond exactly to the true primes.
use num::Integer;

use crate::data::number_types::integer::factorization::prime::buckets::{BUCKETS_32, BUCKETS_64};
use crate::data::number_types::nonzero::Nonzero;

mod buckets;
pub mod primes;

pub trait Prime: Nonzero {
    fn is_prime(&self) -> bool;
}

impl Prime for u64 {
    fn is_prime(&self) -> bool {
        debug_assert!(self.is_not_zero());

        if *self == 2 || *self == 3 || *self == 5 || *self == 7 {
            return true;
        }

        if *self % 2 == 0 || *self % 3 == 0 || *self % 5 == 0 || *self % 7 == 0 {
            return false;
        }

        if *self < 121 {
            return *self > 1;
        }

        if !self.is_probable_prime(2) {
            return false;
        }

        let hash = hash_function_64(*self);
        let basis = BUCKETS_64[hash];

        return self.is_probable_prime((basis & 4095) as u64) && self.is_probable_prime((basis >> 12) as u64);
    }
}

fn hash_function_64(mut hash: u64) -> usize {
    unsafe {
        hash = u64::unchecked_mul((hash >> 32) ^ hash,0x45d9f3b3335b369);
        hash = u64::unchecked_mul((hash >> 32) ^ hash, 0x3335b36945d9f3b);
    }
    hash = (hash >> 32) ^ hash;
    hash &= 0b0000_0000_0000_0000_0011_1111_1111_1111;

    return hash as usize
}


// given 0 <= a,b,n < 2^64, computes (a*b)%n without overflow
fn safe_mul(a: u64, b: u64, n: u64) -> u64 {
    (((a as u128) * b as u128) % n as u128) as u64
}

// given 0 <= a,b,n < 2^64, computes (a^b)%n without overflow
fn safe_exp(mut base: u64, exponent: u64, n: u64) -> u64 {
    let mut result: u64 = 1;
    let mut power: u64 = 1;

    for _ in 0..64 {
        if power > exponent {
            break;
        }

        if exponent & power > 0 {
            result = safe_mul(result, base, n);
        }
        base = safe_mul(base, base, n);
        power <<= 1;
    }

    return result;
}

impl Prime for u32 {
    fn is_prime(&self) -> bool {
        debug_assert!(self.is_not_zero());

        if *self == 2 || *self == 3 || *self == 5 || *self == 7 {
            return true;
        }
        if *self % 2 == 0 || *self % 3 == 0 || *self % 5 == 0 || *self % 7 == 0 {
            return false;
        }
        if *self < 121 {
            return *self > 1;
        }

        let basis_index = hash_function_32(*self as u64);

        return self.is_probable_prime(BUCKETS_32[basis_index] as u32);
    }
}

trait ProbablePrime {
    fn is_probable_prime(&self, base: Self) -> bool;
}

impl ProbablePrime for u32 {
    fn is_probable_prime(&self, base: u32) -> bool {
        debug_assert!(*self > 2);

        let mut d = self - 1;
        let mut s = 0;

        while d.is_even() {
            s += 1;
            d /= 2;
        }

        let mut current: u64 = 1;
        let mut power_to_do = d;
        let mut base_power_d = base;

        while power_to_do > 0 {
            if power_to_do.is_odd() {
                current = (current * base_power_d as u64) % *self as u64;
            }

            base_power_d = ((base_power_d as u64 * base_power_d as u64) % *self as u64) as u32;
            power_to_do >>= 1;
        }

        // base_power_d is now actually equal to base ** d

        if current == 1 {
            return true;
        }

        for _r in 0..s {
            if current == *self as u64 - 1 {
                return true;
            }
            current = (current * current) % *self as u64;
        }

        return false;
    }
}

impl ProbablePrime for u64 {
    // given 2 <= n,a < 2^64, a prime, check whether n is a-SPRP
    fn is_probable_prime(&self, base: Self) -> bool {
        if *self == base {
            return true;
        }

        if *self % base == 0 {
            return false;
        }

        let mut d: u64 = *self - 1;
        let mut s: i32 = 0;
        while d % 2 == 0 {
            s += 1;
            d /= 2;
        }

        let mut current: u64 = safe_exp(base, d, *self);
        if current == 1 {
            return true;
        }

        for _ in 0..s {
            if current == *self - 1 {
                return true;
            }

            current = safe_mul(current, current, *self);
        }

        return false;
    }
}

fn hash_function_32(mut hash: u64) -> usize {
    for _ in 0..2 {
        hash = (hash >> 16) ^ hash;
        unsafe {
            hash = hash.unchecked_mul(0x45d9f3b);
        }
    }
    hash = ((hash >> 16) ^ hash) & 0b0000_0000_0000_0000_0000_0000_1111_1111;

    return hash as usize
}

#[cfg(test)]
mod test {
    use crate::data::number_types::integer::factorization::prime::Prime;

    #[test]
    fn test_is_prime_32() {
        assert!(2_u32.is_prime());
        assert!(11_u32.is_prime());
        assert!(173_u32.is_prime());
        assert!(7741_u32.is_prime());
        assert!(((2_u64.pow(32) - 65) as u32).is_prime());
        assert!(((2_u64.pow(32) - 17) as u32).is_prime());
        assert!(((2_u64.pow(32) - 5) as u32).is_prime());

        assert!(!4_u32.is_prime());
        assert!(!27_u32.is_prime());
        assert!(!1819682_u32.is_prime());
        assert!(!68484865_u32.is_prime());
        assert!(!((2_u64.pow(32) - 64) as u32).is_prime());
        assert!(!((2_u64.pow(32) - 18) as u32).is_prime());
        assert!(!((2_u64.pow(32) - 4) as u32).is_prime());
        assert!(!((2_u64.pow(32) - 6) as u32).is_prime());
    }

    #[test]
    fn test_is_prime_64() {
        assert!(2_u64.is_prime());
        assert!(11_u64.is_prime());
        assert!(173_u64.is_prime());
        assert!(7741_u64.is_prime());
        assert!(((2_u128.pow(64) - 179) as u64).is_prime());
        assert!(((2_u128.pow(64) - 95) as u64).is_prime());
        assert!(((2_u128.pow(64) - 83) as u64).is_prime());
        assert!(((2_u128.pow(64) - 59) as u64).is_prime());

        assert!(!4_u64.is_prime());
        assert!(!27_u64.is_prime());
        assert!(!1819682_u64.is_prime());
        assert!(!68484865_u64.is_prime());
        assert!(!((2_u128.pow(64) - 82) as u64).is_prime());
        assert!(!((2_u128.pow(64) - 3) as u64).is_prime());
        assert!(!((2_u128.pow(64) - 1) as u64).is_prime());
    }

    #[test]
    fn test_count_primes_32() {
        assert_eq!((1_u32..1_000_000_u32).filter(Prime::is_prime).count(), 78498);
        assert_eq!((1_u32..10_000_000_u32).filter(Prime::is_prime).count(), 664579);
    }

    #[test]
    fn test_count_primes64() {
        assert_eq!((1_u64..1_000_000_u64).filter(Prime::is_prime).count(), 78498);
        assert_eq!((1_u64..10_000_000_u64).filter(Prime::is_prime).count(), 664579);
    }
}
