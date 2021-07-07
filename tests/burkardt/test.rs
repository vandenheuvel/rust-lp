use std::convert::TryInto;
use std::ops::{Add, Neg};

use num_traits::{One, Zero};

use relp::algorithm::{OptimizationResult, SolveRelaxation};
use relp::algorithm::two_phase::matrix_provider::MatrixProvider;
use relp::algorithm::two_phase::tableau::inverse_maintenance::carry::basis_inverse_rows::BasisInverseRows;
use relp::algorithm::two_phase::tableau::inverse_maintenance::carry::Carry;
use relp::algorithm::two_phase::tableau::inverse_maintenance::carry::lower_upper::LUDecomposition;
use relp::data::linear_algebra::traits::Element;
use relp::data::linear_program::general_form::GeneralForm;
use relp::data::linear_program::solution::Solution;
use relp::io::import;
use relp_num::{Rational64, RationalBig, NonZero};
use relp_num::OrderedField;
use relp_num::RB;

use super::get_test_file_path;
use std::str::FromStr;

fn to_general_form<T: From<Rational64> + Zero + NonZero + One + Ord + Element + Neg<Output=T>>(
    file_name: &str,
) -> GeneralForm<T>
where
    for<'r> T: Add<&'r T, Output=T>,
    for<'r> &'r T: Neg<Output=T>,
{
    let path = get_test_file_path(file_name);
    let result = import(&path).unwrap();

    result.try_into().unwrap()
}

#[test]
fn adlittle() {
    type T = RationalBig;
    type S = RationalBig;

    let path = get_test_file_path("adlittle");
    let result = import::<T>(&path).unwrap();

    let mut general = result.try_into().unwrap();
    general.presolve().unwrap();
    let constraint_type_counts = general.standardize();
    let data = general.derive_matrix_data(constraint_type_counts);
    let result = data.solve_relaxation::<Carry<S, BasisInverseRows<_>>>();

    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let reconstructed = data.reconstruct_solution(vector);
            let solution = general.compute_full_solution_with_reduced_solution(reconstructed);
            assert_eq!(solution.objective_value, S::from_str("24975305659811992079614961229/120651674036153428931840").unwrap());
        },
        _ => assert!(false),
    }
}

#[test]
fn afiro() {
    type T = Rational64;
    type S = RationalBig;

    let mut general = to_general_form::<T>("afiro");
    general.presolve().unwrap();
    let constraint_type_counts = general.standardize();
    let data = general.derive_matrix_data(constraint_type_counts);
    let result = data.solve_relaxation::<Carry<_, LUDecomposition<S>>>();

    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let reconstructed = data.reconstruct_solution(vector);
            let solution = general.compute_full_solution_with_reduced_solution(reconstructed);
            let expected = Solution::new(
                RB!(-406_659, 875),  // GLPK
                vec![
                    ("X01".to_string(), RB!(80, 1)),
                    ("X02".to_string(), RB!(51, 2)),
                    ("X03".to_string(), RB!(109, 2)),
                    ("X04".to_string(), RB!(424, 5)),
                    ("X06".to_string(), RB!(255, 14)),
                    ("X07".to_string(), S::zero()),
                    ("X08".to_string(), S::zero()),
                    ("X09".to_string(), S::zero()),
                    ("X10".to_string(), S::zero()),
                    ("X11".to_string(), S::zero()),
                    ("X12".to_string(), S::zero()),
                    ("X13".to_string(), S::zero()),
                    ("X14".to_string(), RB!(255, 14)),
                    ("X15".to_string(), S::zero()),
                    ("X16".to_string(), RB!(999, 1)),
                    ("X22".to_string(), RB!(500, 1)),
                    ("X23".to_string(), RB!(11898, 25)),
                    ("X24".to_string(), RB!(602, 25)),
                    ("X25".to_string(), S::zero()),
                    ("X26".to_string(), RB!(215, 1)),
                    ("X28".to_string(), S::zero()),
                    ("X29".to_string(), S::zero()),
                    ("X30".to_string(), S::zero()),
                    ("X31".to_string(), S::zero()),
                    ("X32".to_string(), S::zero()),
                    ("X33".to_string(), S::zero()),
                    ("X34".to_string(), S::zero()),
                    ("X35".to_string(), S::zero()),
                    ("X36".to_string(), RB!(11898, 35)),
                    ("X37".to_string(), RB!(11898, 35)),
                    ("X38".to_string(), S::zero()),
                    ("X39".to_string(), S::zero()),
                ],
            );

            assert!(expected.is_probably_equal_to(&solution, 0.1_f64));
        },
        _ => assert!(false),
    }
}

#[test]
#[ignore = "Not yet implemented: The same range value occurring twice for a single row while being equal should be accepted."]
fn empstest() {
    type T = Rational64;

    let path = get_test_file_path("empstest");
    import::<T>(&path).unwrap();
}

#[test]
fn maros() {
    type T = Rational64;
    type S = RationalBig;

    let mut general = to_general_form::<T>("maros");
    general.presolve().unwrap();
    let constraint_type_counts = general.standardize();
    let data = general.derive_matrix_data(constraint_type_counts);
    let result = data.solve_relaxation::<Carry<S, LUDecomposition<_>>>();

    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let reconstructed = data.reconstruct_solution(vector);
            let solution = general.compute_full_solution_with_reduced_solution(reconstructed);
            assert_eq!(solution, Solution::new(
                S::new(385, 3).unwrap(),  // GLPK
                vec![
                    ("VOL1".to_string(), S::new(10, 3).unwrap()),
                    ("VOL2".to_string(), S::new(40, 3).unwrap()),
                    ("VOL3".to_string(), S::new(20, 1).unwrap()),
                    ("VOL4".to_string(), S::new(0, 1).unwrap()),
                ],
            ));
        },
        _ => assert!(false),
    }
}

#[test]
fn nazareth() {
    type T = Rational64;
    type S = RationalBig;

    let mut general = to_general_form::<T>("nazareth");
    general.presolve().unwrap();
    let constraint_type_counts = general.standardize();
    let data = general.derive_matrix_data(constraint_type_counts);
    let result = data.solve_relaxation::<Carry<S, LUDecomposition<_>>>();
    assert_eq!(result, OptimizationResult::Unbounded);  // GLPK
}

#[test]
fn testprob() {
    type T = Rational64;
    type S = RationalBig;

    let mut general = to_general_form::<T>("testprob");
    general.presolve().unwrap();
    let constraint_type_counts = general.standardize();
    let data = general.derive_matrix_data(constraint_type_counts);
    let result = data.solve_relaxation::<Carry<S, LUDecomposition<_>>>();

    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let reconstructed = data.reconstruct_solution(vector);
            let solution = general.compute_full_solution_with_reduced_solution(reconstructed);
            assert_eq!(solution, Solution::new(
                S::new(54, 1).unwrap(),  // GLPK
                vec![
                    ("X1".to_string(), S::new(4, 1).unwrap()),
                    ("X2".to_string(), S::new(-1, 1).unwrap()),
                    ("X3".to_string(), S::new(6, 1).unwrap()),
                ],
            ));
        },
        _ => assert!(false),
    }
}
