[![Build Status](https://travis-ci.com/vandenheuvel/rust-lp.svg?branch=master)](https://travis-ci.com/vandenheuvel/rust-lp) [![codecov](https://codecov.io/gh/vandenheuvel/rust-lp/branch/master/graph/badge.svg)](https://codecov.io/gh/vandenheuvel/rust-lp)

# rust-lp
A linear program solver written in Rust.

This project is a work in progress. Integer variables are not supported.

## Priorities

- Correctness and safety 
- Performance
- Maintainability and flexibility

## Roadmap

- [x] Simplex algorithm
- [x] Lazy column generation
- [x] Arbitrary ordered fields
- [x] Arbitrary precision
- [x] Presolving framework
- [x] Reading MPS files
- [ ] Smart matrix inverse maintenance
- [ ] Floating point support
- [ ] A convenient API
- [ ] Dual algorithm
- [ ] Integration into branch and bound frameworks

## Usage

Notice the "A convenient API" box above is still unchecked...
```rust
extern crate rust_lp as lp;
extern crate num;
```
<details><summary>...</summary>
<p>
    
```rust
use std::env;
use std::path::Path;

use lp::algorithm::simplex::solve_relaxation;
use lp::io::import;
use lp::algorithm::simplex::logic::OptimizationResult;

use num_traits::rational::Ratio;
use num_traits::BigInt;
use lp::data::linear_program::general_form::GeneralForm;
use std::convert::TryInto;
use lp::algorithm::simplex::matrix_provider::MatrixProvider;
use lp::algorithm::simplex::strategy::pivot_rule::{FirstProfitable, FirstProfitableWithMemory};
```
</p>
</details>

```rust
fn main() {
    type T = Ratio<BigInt>;

    let name = env::args().nth(1).unwrap();
    let pathname = format!("./test_files/{}.mps", name);
    let path = Path::new(&pathname);
    let result = import(path).unwrap();

    let mut general_form: GeneralForm<T, T> = result.try_into().ok().unwrap();
    let matrix_data = general_form.derive_matrix_data().ok().unwrap();

    let result = solve_relaxation::<_, _, _, FirstProfitableWithMemory, FirstProfitableWithMemory>(&matrix_data);

    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let reconstructed = matrix_data.reconstruct_solution(vector);
            let solution = general_form.compute_full_solution_with_reduced_solution(reconstructed);
            println!("{:?}", solution.objective_value);
        },
        _ => panic!(),
    }
}
```
