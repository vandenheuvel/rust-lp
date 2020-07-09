use std::convert::TryInto;

use num_traits::{BigInt, ToPrimitive};
use num_traits::rational::Ratio;

use rust_lp::algorithm::two_phase::matrix_provider::MatrixProvider;
use rust_lp::algorithm::two_phase::tableau::inverse_maintenance::carry::Carry;
use rust_lp::io::import;

use super::get_test_file_path;

/// Testing a problem, comparing only the objective value.
fn test(file_name: &str, objective: f64) {
    type T = Ratio<BigInt>;

    let path = get_test_file_path(file_name);
    let result = import::<T>(&path).unwrap();

    let mut general = result.try_into().ok().unwrap();
    let data = general.derive_matrix_data().ok().unwrap();
    let result = data.solve_relaxation::<Carry<_>>();

    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let reconstructed = data.reconstruct_solution(vector);
            let solution = general.compute_full_solution_with_reduced_solution(reconstructed);
            assert_eq!(solution.objective_value.floor().to_integer().to_i64().unwrap(), objective.floor() as i64);
        },
        _ => assert!(false),
    }
}

#[test]
fn test_50v() {
    test("50v-10", 2879.065687f64);
}  // GLPK

#[test]
#[ignore = "Incorrectly determined infeasible"]
fn test_30n() {
    test("30n20b8", 43.33557298f64);
}  // GLPK

#[test]
#[ignore = "Too computationally intensive"]
fn test_acc() {
    test("acc-tight4", 0f64);
}  // GLPK
