use std::num::{NonZeroI8, NonZeroU8};
use crate::data::number_types::traits::factorization::{NonzeroFactorizable, NonzeroFactorization};
use crate::data::number_types::nonzero::{NonzeroSign, NonzeroSigned};
use crate::data::number_types::integer::factorization::buckets::BUCKETS_32;


mod buckets;


const MILLER_RABIN_BASES_32: [u32; 3] = [2, 7, 61];
const MILLER_RABIN_BASES_64: [u64; 7] = [2, 325, 9375, 28178, 450775, 9780504, 1795265022];


// const SMALL_PRIMES: [u8; 64] = [
//     2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
//     73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
//     179, 181, 191, 193, 197, 199, 211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
//     283, 293, 307, 311,
// ];

impl NonzeroFactorizable for u64 {
    type Factor = u8;
    // TODO(CORRECTNESS): Consider making this type unsigned
    type Power = i8;

    fn factorize(&self) -> NonzeroFactorization<Self::Factor, Self::Power> {
        unimplemented!();
        // let mut x = self.clone();
        //
        // let sign = NonzeroSigned::signum(&x);
        //
        // let mut factors = Vec::new();
        //
        // // TODO: Use an external library instead
        // // Trial division
        // for divisor in SMALL_PRIMES {
        //     divide_out_factor(&mut x, divisor, &mut factors);
        //     if divisor == 1 {
        //         break;
        //     }
        // }
        // for divisor in ((SMALL_PRIMES[63] + 2)..(2_u8.pow(8 / 2))).step_by(2) {
        //     divide_out_factor(&mut x, divisor, &mut factors);
        // }
        //
        // NonzeroFactorization { sign, factors }
    }
}

// fn divide_out_factor(x: &mut i8, divisor: u8, factors: &mut [(u8, u8)]) {
//     let mut division_count = 0;
//     while x % divisor == 0 {
//         *x /= divisor;
//         division_count += 1;
//     }
//
//     if division_count > 0 {
//         factors.push((divisor.clone(), division_count));
//     }
// }

fn is_sprp(n: u32, mut a: u32) -> bool {
    let mut d = n - 1;
    let mut s = 0;
    while d & 1 == 0 {
        s += 1;
        d >>= 1;
    }

    let mut current: u64 = 1;
    let mut pw = d;

    while pw > 0 {
        if pw & 1 > 0 {
            current = (current * a as u64) % n as u64;
        }

        a = ((a as u64 * a as u64) % n as u64) as u32;
        pw >>= 1;
    }

    if current == 1 {
        return true;
    }

    for _ in 0..s {
        if current == n as u64 - 1 {
            return true;
        }
        current = (current * current) % n as u64;
    }

    return false;
}

pub fn is_prime(x: u32) -> bool {
    if x == 2 || x == 3 || x == 5 || x == 7 {
        return true;
    }
    if x % 2 == 0 || x % 3 == 0 || x % 5 == 0 || x % 7 == 0 {
        return false;
    }
    if x < 121 {
        return x > 1;
    }

    let basis_index = hash_function(x as u64);
    let basis = BUCKETS_32[basis_index];

    return is_sprp(x, basis as u32);
}

fn hash_function(mut hash: u64) -> usize {
    for _ in 0..2 {
        hash = ((hash >> 16) ^ hash);
        unsafe {
            hash = hash.unchecked_mul(0x45d9f3b);
        }
    }
    hash = ((hash >> 16) ^ hash) & 255;

    return hash as usize
}

macro_rules! impl_factorization_non_zero {
    ($t:ident, $f:ident) => {
        impl NonzeroFactorizable for $t {
            type Factor = $f;
            // TODO(CORRECTNESS): Consider making this type unsigned
            type Power = i8;

            fn factorize(&self) -> NonzeroFactorization<Self::Factor, Self::Power> {
                let sign = self.signum();

                let mut factors = Vec::new();


                NonzeroFactorization { sign,factors }
            }
        }
    }
}

// impl_factorization_non_zero!(NonZeroI8, NonZeroU8);
impl_factorization_non_zero!(i8, u8);
impl_factorization_non_zero!(u8, u8);
impl_factorization_non_zero!(i16, u16);
impl_factorization_non_zero!(u16, u16);
impl_factorization_non_zero!(i32, u32);
impl_factorization_non_zero!(u32, u32);
impl_factorization_non_zero!(i64, u64);
// impl_factorization_non_zero!(u64, u64);
impl_factorization_non_zero!(i128, u128);
impl_factorization_non_zero!(u128, u128);

#[cfg(test)]
mod test {
    use crate::data::number_types::traits::factorization::{NonzeroFactorizable, NonzeroFactorization};
    use crate::data::number_types::nonzero::sign::Sign;
    use crate::data::number_types::integer::factorization::is_prime;

    #[test]
    fn test_factorize() {
        assert_eq!(2_u8.factorize(), NonzeroFactorization { sign: Sign::Positive, factors: vec![(2, 1)] });
    }

    #[test]
    fn test_is_prime() {
        assert_eq!(is_prime(2), true);
        assert_eq!(is_prime(11), true);
        assert_eq!(is_prime(173), true);
        assert_eq!(is_prime(7741), true);
        assert_eq!(is_prime(7741), true);
        assert_eq!(is_prime((2_u64.pow(32) - 65) as u32), true);
        assert_eq!(is_prime((2_u64.pow(32) - 17) as u32), true);
        assert_eq!(is_prime((2_u64.pow(32) - 5) as u32), true);

        assert_eq!(is_prime(4), false);
        assert_eq!(is_prime(27), false);
        assert_eq!(is_prime(1819682), false);
        assert_eq!(is_prime(68484865), false);
        assert_eq!(is_prime((2_u64.pow(32) - 64) as u32), false);
        assert_eq!(is_prime((2_u64.pow(32) - 18) as u32), false);
        assert_eq!(is_prime((2_u64.pow(32) - 4) as u32), false);
        assert_eq!(is_prime((2_u64.pow(32) - 6) as u32), false);
    }

    #[test]
    fn test_count_primes() {
        assert_eq!((1..1_000_000).filter(|&n| is_prime(n)).count(), 78498);
    }
}
