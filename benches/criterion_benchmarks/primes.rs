use criterion::{black_box, Criterion, criterion_group, criterion_main};

use rust_lp::data::number_types::integer::factorization::prime::Prime;

pub fn is_prime_small_positive(c: &mut Criterion) {
    c.bench_function("is_prime with small prime", |b| b.iter(|| {
        Prime::is_prime(black_box(&31_u64))
    }));
}

pub fn is_prime_small_negative(c: &mut Criterion) {
    c.bench_function("is_prime with small nonprime", |b| b.iter(|| {
        Prime::is_prime(black_box(&30_u64))
    }));
}

pub fn is_prime_large_positive(c: &mut Criterion) {
    c.bench_function("is_prime with large prime", |b| b.iter(|| {
        Prime::is_prime(black_box(&(2_u64.pow(32) - 5)))
    }));
}

pub fn is_prime_large_negative(c: &mut Criterion) {
    c.bench_function("is_prime with large nonprime", |b| b.iter(|| {
        Prime::is_prime(black_box(&(2_u64.pow(32) - 21)))
    }));
}

criterion_group!(primes,
    is_prime_small_positive,
    is_prime_small_negative,
    is_prime_large_positive,
    is_prime_large_negative,
);
