//! # The Mixed Integer Programming Library tests
//!
//! This module integration-tests the project on a number of problem from the
//! [MIPLIB](http://miplib.zib.de/) project.
//!
//! ## Note
//!
//! The tests in this module are only ran when the `miplib_tests` feature is enabled as these tests
//! take a long time (minutes to hours) to run.
use std::convert::TryInto;
use std::path::{Path, PathBuf};

use num::{BigInt, ToPrimitive};
use num::rational::Ratio;

use rust_lp::algorithm::{OptimizationResult, SolveRelaxation};
use rust_lp::algorithm::two_phase::matrix_provider::MatrixProvider;
use rust_lp::algorithm::two_phase::tableau::inverse_maintenance::carry::Carry;
use rust_lp::io::import;

/// # Generation and execution
#[allow(missing_docs)]
mod test;

/// Relative path of the folder where the mps files are stored.
///
/// The path is relative to the project root folder.
fn problem_file_directory() -> PathBuf {
    Path::new(file!()).parent().unwrap().join("problem_files")
}

/// Compute the path of the problem file, based on the problem name.
///
/// # Arguments
///
/// * `name`: Problem name without extension.
///
/// # Return value
///
/// File path relative to the project root folder.
fn get_test_file_path(name: &str) -> PathBuf {
    problem_file_directory().join(name).with_extension("mps")
}


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
