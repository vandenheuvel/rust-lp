//! # Netlib
//!
//! Hosted [here](http://www.numerical.rl.ac.uk/cute/netlib.html).
use std::convert::TryInto;
use std::fs::File;
use std::io::Read;
use std::path::{Path, PathBuf};

use relp::algorithm::{OptimizationResult, SolveRelaxation};
use relp::algorithm::two_phase::matrix_provider::MatrixProvider;
use relp::algorithm::two_phase::tableau::inverse_maintenance::carry::basis_inverse_rows::BasisInverseRows;
use relp::algorithm::two_phase::tableau::inverse_maintenance::carry::Carry;
use relp::data::linear_program::general_form::{GeneralForm, Scalable};
use relp::data::linear_program::solution::Solution;
use relp::io::error::Import;
use relp::io::mps::parse_fixed;
use relp_num::{RationalBig, NonZeroFactorizable};
use std::ops::{MulAssign, DivAssign, Mul, Neg, Sub, Add, AddAssign, Div};
use num_traits::{One, Zero};
use std::hash::Hash;

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
    problem_file_directory().join(name).with_extension("SIF")
}

type T = RationalBig;
type S = RationalBig;

fn solve(file_name: &str) -> Solution<S> {
    let file_path = get_test_file_path(file_name);

    let mut program = String::new();
    File::open(&file_path)
        .map_err(Import::IO).unwrap()
        .read_to_string(&mut program)
        .map_err(Import::IO).unwrap();
    let mps = parse_fixed(&program).unwrap();

    let mut general: GeneralForm<T> = mps.try_into().unwrap();
    general.presolve().unwrap();
    let constraint_type_counts = general.standardize();
    let scaling = Scalable::<S>::scale(&mut general);
    f::<S>();

    println!("{:?}", scaling);
    let data = general.derive_matrix_data(constraint_type_counts);
    let result = data.solve_relaxation::<Carry<S, BasisInverseRows<_>>>();

    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let reconstructed = data.reconstruct_solution(vector);
            general.compute_full_solution_with_reduced_solution(reconstructed)
        },
        _ => panic!(),
    }
}
use relp::data::linear_algebra::traits::SparseElement;
use relp::data::linear_algebra::traits::SparseComparator;
fn f<R>()
where
    for<'r> R:
        NonZeroFactorizable +
        MulAssign<R::Factor> +
        DivAssign<R::Factor> +
        MulAssign<&'r R::Factor> +
        DivAssign<&'r R::Factor> +
        MulAssign<R> +
        DivAssign<R> +
        MulAssign<&'r R> +
        DivAssign<&'r R> +
        Mul<&'r R, Output=R> +
        One +
        SparseElement<R> +
        SparseComparator +
    ,
    for<'r> &'r R: Mul<&'r R, Output=R>,
    R::Factor: Hash,
    R::Power: Neg<Output=R::Power> + Sub<Output=R::Power> + Add<Output=R::Power> + Zero + One + AddAssign + Ord,
    for<'r> &'r R: Div<&'r R, Output=R>,
{
}
