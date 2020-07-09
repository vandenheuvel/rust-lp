use std::convert::TryInto;

use num_traits::{BigInt, FromPrimitive};
use num_traits::rational::Ratio;

use rust_lp::{BR, R128};
use rust_lp::algorithm::two_phase::matrix_provider::MatrixProvider;
use rust_lp::algorithm::two_phase::tableau::inverse_maintenance::carry::Carry;
use rust_lp::data::linear_program::solution::Solution;
use rust_lp::data::number_types::traits::{OrderedField};
use rust_lp::io::import;

use super::get_test_file_path;

fn solve<T: 'static + OrderedField + FromPrimitive>(
    file_name: &str,
) -> Solution<T>
where
    T: OrderedFieldRef<T>,
{
    let path = get_test_file_path(file_name);
    let mps = import::<T>(&path).unwrap();

    let mut general = mps.try_into().ok().unwrap();
    let data = general.derive_matrix_data().ok().unwrap();
    let result = data.solve_relaxation::<Carry<_>>();

    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let reconstructed = data.reconstruct_solution(vector);
            general.compute_full_solution_with_reduced_solution(reconstructed)
        },
        _ => panic!(),
    }
}

#[test]
#[ignore = "Has a row name consisting of only numbers, the parser doesn't support that."]
fn test_25FV47() {
    type T = Ratio<i128>;

    let result = solve::<T>("25FV47");
    assert_eq!(result.objective_value, R128!(5.5018459e+03)); // Gurobi
}

#[test]
fn test_80BAU3B() {
    type T = Ratio<BigInt>;

    let result = solve::<T>("80BAU3B");
    assert_eq!(result.objective_value, BR!(9.872241924e+05)); // Gurobi
}

#[test]
#[ignore = "Too computationally expensive."]
fn test_ADLITTLE() {
    type T = Ratio<BigInt>;

    let result = solve::<T>("ADLITTLE");
    assert_eq!(result.objective_value, BR!(2.254949632e+05)); // Gurobi
}
