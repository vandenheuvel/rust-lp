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

use relp::algorithm::{OptimizationResult, SolveRelaxation};
use relp::algorithm::two_phase::matrix_provider::MatrixProvider;
use relp::algorithm::two_phase::tableau::inverse_maintenance::carry::basis_inverse_rows::BasisInverseRows;
use relp::algorithm::two_phase::tableau::inverse_maintenance::carry::Carry;
use relp_num::RationalBig;
use relp_num::Abs;
use relp::io::import;
use relp_num::RB;

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
fn test(file_name: &str, objective: f64, tolerance: f64) {
    type T = RationalBig;
    type S = RationalBig;

    let path = get_test_file_path(file_name);
    let result = import::<T>(&path).unwrap();

    let mut general = result.try_into().ok().unwrap();
    general.presolve().unwrap();
    general.standardize();
    let data = general.derive_matrix_data();
    let result = data.solve_relaxation::<Carry<S, BasisInverseRows<_>>>();

    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let reconstructed = data.reconstruct_solution(vector);
            let solution = general.compute_full_solution_with_reduced_solution(reconstructed);

            // TODO(CORRECTNESS): Make this a relative error check.
            assert!((solution.objective_value - RB!(objective)).abs() < RB!(tolerance));
        },
        _ => assert!(false),
    }
}
