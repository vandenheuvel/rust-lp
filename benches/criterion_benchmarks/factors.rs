use criterion::{black_box, Criterion, criterion_group};

use rust_lp::data::number_types::traits::factorization::{NonZeroFactorizable, NonZeroFactorization};

pub fn factor_small(c: &mut Criterion) {
    c.bench_function("factorize a small number", |b| b.iter(|| {
        NonZeroFactorizable::factorize(black_box(&31_u64))
    }));
}

criterion_group!(factors,
    factor_small,
);
