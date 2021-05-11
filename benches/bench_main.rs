use criterion::criterion_main;

mod criterion_benchmarks;

criterion_main! {
    criterion_benchmarks::primes::primes,
    criterion_benchmarks::factors::factors,
}
