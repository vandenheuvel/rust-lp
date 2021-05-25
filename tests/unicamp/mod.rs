//! # Problems from Unicamp
//!
//! Hosted [here](https://www.cenapad.unicamp.br/parque/manuais/OSL/oslweb/features/feat24DT.htm).
use std::convert::TryInto;
use std::ops::{Add, AddAssign, Sub};
use std::path::{Path, PathBuf};

use num_traits::{One, Zero};

use relp::algorithm::{OptimizationResult, SolveRelaxation};
use relp::algorithm::two_phase::matrix_provider::MatrixProvider;
use relp::algorithm::two_phase::tableau::inverse_maintenance::carry::Carry;
use relp::algorithm::two_phase::tableau::inverse_maintenance::carry::lower_upper::LUDecomposition;
use relp::algorithm::two_phase::tableau::inverse_maintenance::ops as im_ops;
use relp::algorithm::two_phase::tableau::kind::artificial::Cost;
use relp::data::linear_algebra::traits::Element;
use relp::data::linear_program::elements::LinearProgramType;
use relp::data::linear_program::solution::Solution;
use relp::data::number_types::rational::Rational64;
use relp::data::number_types::traits::{OrderedField, OrderedFieldRef};
use relp::io::import;

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

fn solve<
    IMT: im_ops::Field + im_ops::FieldHR + im_ops::Column<GFT> + AddAssign<GFT> + PartialEq<GFT> + Ord,
    GFT: 'static + From<Rational64> + Zero + One + Ord + Element + OrderedField,
>(file_name: &str) -> Solution<IMT>
where
    for<'r> &'r GFT: OrderedFieldRef<GFT>,
    for<'r> &'r IMT: Add<&'r GFT, Output=IMT> + Sub<&'r IMT, Output=IMT>,
    for<'r> IMT: im_ops::Cost<Option<&'r GFT>> + im_ops::Cost<Cost>,
{
    let path = get_test_file_path(file_name);
    let mps = import::<GFT>(&path).unwrap();

    let mut general = mps.try_into().unwrap();
    let data = match general.derive_matrix_data() {
        Ok(data) => data,
        Err(LinearProgramType::FiniteOptimum(Solution {
                                                 objective_value, solution_values,
                                             })) => {
            return Solution {
                objective_value: objective_value.into(),
                solution_values: solution_values.into_iter()
                    .map(|(name, value)| (name, value.into()))
                    .collect(),
            }
        },
        _ => panic!(),
    };
    let result = data.solve_relaxation::<Carry<IMT, LUDecomposition<_>>>();

    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let reconstructed = data.reconstruct_solution(vector);
            general.compute_full_solution_with_reduced_solution::<IMT>(reconstructed)
        },
        _ => panic!(),
    }
}
