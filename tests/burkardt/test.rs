use std::convert::TryInto;

use num::{FromPrimitive, One, Zero};

use rust_lp::algorithm::{OptimizationResult, SolveRelaxation};
use rust_lp::algorithm::two_phase::matrix_provider::MatrixProvider;
use rust_lp::algorithm::two_phase::tableau::inverse_maintenance::carry::Carry;
use rust_lp::data::linear_algebra::traits::Element;
use rust_lp::data::linear_program::general_form::GeneralForm;
use rust_lp::data::linear_program::solution::Solution;
use rust_lp::data::number_types::rational::{Rational64, RationalBig};
use rust_lp::data::number_types::traits::OrderedField;
use rust_lp::io::import;
use rust_lp::RB;

use super::get_test_file_path;

fn to_general_form<T: From<Rational64> + Zero + One + Ord + Element>(
    file_name: &str,
) -> GeneralForm<T> {
    let path = get_test_file_path(file_name);
    let result = import(&path).unwrap();

    result.try_into().ok().unwrap()
}

#[test]
#[ignore = "Takes too long, fill in the right values"]
fn adlittle() {
    type T = RationalBig;

    let path = get_test_file_path("adlittle");
    let result = import::<T>(&path).unwrap();

    let mut general = result.try_into().ok().unwrap();
    let data = general.derive_matrix_data().ok().unwrap();
    let result = data.solve_relaxation::<Carry<_>>();

    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let reconstructed = data.reconstruct_solution(vector);
            let solution = general.compute_full_solution_with_reduced_solution(reconstructed);
            let expected = Solution::new(
                RB!(207003, 1),  // GLPK
                vec![
                    ("...100".to_string(), RB!(0, 1)),
                ],
            );

            assert!(expected.is_probably_equal_to(&solution, 0.1_f64));
        },
        _ => assert!(false),
    }
}

#[test]
#[ignore = "Fill in the accurate values"]
fn afiro_big_int() {
    type T = RationalBig;

    let mut general = to_general_form::<T>("afiro");
    let data = general.derive_matrix_data().ok().unwrap();
    let result = data.solve_relaxation::<Carry<_>>();

    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let reconstructed = data.reconstruct_solution(vector);
            let solution = general.compute_full_solution_with_reduced_solution(reconstructed);
            let expected = Solution::new(
                RB!(-406659, 875),  // GLPK
                vec![
                    ("X01".to_string(), RB!(80, 1)),
                    ("X02".to_string(), RB!(51, 2)),
                    ("X03".to_string(), RB!(109, 2)),
                    ("X04".to_string(), RB!(424, 5)),
                    ("X06".to_string(), RB!(255, 14)),
                    ("X07".to_string(), T::zero()),
                    ("X08".to_string(), T::zero()),
                    ("X09".to_string(), T::zero()),
                    ("X10".to_string(), T::zero()),
                    ("X11".to_string(), T::zero()),
                    ("X12".to_string(), T::zero()),
                    ("X13".to_string(), T::zero()),
                    ("X14".to_string(), RB!(255, 14)),
                    ("X15".to_string(), T::zero()),
                    ("X16".to_string(), RB!(999, 1)),
                    ("X22".to_string(), RB!(500, 1)),
                    ("X23".to_string(), RB!(11898, 25)),
                    ("X24".to_string(), RB!(602, 25)),
                    ("X25".to_string(), T::zero()),
                    ("X26".to_string(), RB!(215, 1)),
                    ("X28".to_string(), T::zero()),
                    ("X29".to_string(), T::zero()),
                    ("X30".to_string(), T::zero()),
                    ("X31".to_string(), T::zero()),
                    ("X32".to_string(), T::zero()),
                    ("X33".to_string(), T::zero()),
                    ("X34".to_string(), T::zero()),
                    ("X35".to_string(), T::zero()),
                    ("X36".to_string(), RB!(11898, 35)),
                    ("X37".to_string(), RB!(11898, 35)),
                    ("X38".to_string(), T::zero()),
                    ("X39".to_string(), T::zero()),
                ],
            );

            assert!(expected.is_probably_equal_to(&solution, 0.1_f64));
        },
        _ => assert!(false),
    }
}

/// TODO: Why does the bigint give a different solution?
#[test]
#[ignore = "Different values than the big integer variant of this test."]
fn afiro() {
    type T = RationalBig;

    let mut general = to_general_form::<T>("afiro");
    let data = general.derive_matrix_data().ok().unwrap();
    let result = data.solve_relaxation::<Carry<_>>();

    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let reconstructed = data.reconstruct_solution(vector);
            let solution = general.compute_full_solution_with_reduced_solution(reconstructed);
            let expected = Solution::new(
                T::new(-406659, 875),  // GLPK
                vec![
                    ("X01".to_string(), T::new(80, 1)),
                    ("X02".to_string(), T::new(51, 2)),
                    ("X03".to_string(), T::new(109, 2)),
                    ("X04".to_string(), T::new(424, 5)),
                    ("X06".to_string(), T::new(255, 14)),
                    ("X07".to_string(), T::zero()),
                    ("X08".to_string(), T::zero()),
                    ("X09".to_string(), T::zero()),
                    ("X10".to_string(), T::zero()),
                    ("X11".to_string(), T::zero()),
                    ("X12".to_string(), T::zero()),
                    ("X13".to_string(), T::zero()),
                    ("X14".to_string(), T::new(255, 14)),
                    ("X15".to_string(), T::zero()),
                    ("X16".to_string(), T::new(9999999999999999, 1)),
                    ("X22".to_string(), T::new(500, 1)),
                    ("X23".to_string(), T::new(11898, 25)),
                    ("X24".to_string(), T::new(602, 25)),
                    ("X25".to_string(), T::zero()),
                    ("X26".to_string(), T::new(215, 1)),
                    ("X28".to_string(), T::zero()),
                    ("X29".to_string(), T::zero()),
                    ("X30".to_string(), T::zero()),
                    ("X31".to_string(), T::zero()),
                    ("X32".to_string(), T::zero()),
                    ("X33".to_string(), T::zero()),
                    ("X34".to_string(), T::zero()),
                    ("X35".to_string(), T::zero()),
                    ("X36".to_string(), T::new(11898, 35)),
                    ("X37".to_string(), T::new(11898, 35)),
                    ("X38".to_string(), T::zero()),
                    ("X39".to_string(), T::zero()),
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
fn maros_bigint() {
    type T = RationalBig;

    let mut general = to_general_form::<T>("maros");
    let data = general.derive_matrix_data().ok().unwrap();
    let result = data.solve_relaxation::<Carry<_>>();

    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let reconstructed = data.reconstruct_solution(vector);
            let solution = general.compute_full_solution_with_reduced_solution(reconstructed);
            assert_eq!(solution, Solution::new(
                RB!(385, 3),  // GLPK
                vec![
                    ("VOL1".to_string(), RB!(10, 3)),
                    ("VOL2".to_string(), RB!(40, 3)),
                    ("VOL3".to_string(), RB!(20, 1)),
                    ("VOL4".to_string(), RB!(0, 1)),
                ],
            ));
        },
        _ => assert!(false),
    }
}

#[test]
fn maros() {
    type T = Rational64;

    let mut general = to_general_form::<T>("maros");
    let data = general.derive_matrix_data().ok().unwrap();
    let result = data.solve_relaxation::<Carry<_>>();

    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let reconstructed = data.reconstruct_solution(vector);
            let solution = general.compute_full_solution_with_reduced_solution(reconstructed);
            assert_eq!(solution, Solution::new(
                T::new(385, 3),  // GLPK
                vec![
                    ("VOL1".to_string(), T::new(10, 3)),
                    ("VOL2".to_string(), T::new(40, 3)),
                    ("VOL3".to_string(), T::new(20, 1)),
                    ("VOL4".to_string(), T::new(0, 1)),
                ],
            ));
        },
        _ => assert!(false),
    }
}

#[test]
fn nazareth_bigint() {
    type T = RationalBig;

    let mut general = to_general_form::<T>("nazareth");
    let data = general.derive_matrix_data().ok().unwrap();
    let result = data.solve_relaxation::<Carry<T>>();
    assert_eq!(result, OptimizationResult::Unbounded);  // GLPK
}

#[test]
fn nazareth() {
    type T = Rational64;

    let mut general = to_general_form::<T>("nazareth");
    let data = general.derive_matrix_data().ok().unwrap();
    let result = data.solve_relaxation::<Carry<T>>();
    assert_eq!(result, OptimizationResult::Unbounded);  // GLPK
}

#[test]
fn testprob_bigint() {
    type T = RationalBig;

    let mut general = to_general_form::<T>("testprob");
    let data = general.derive_matrix_data().ok().unwrap();
    let result = data.solve_relaxation::<Carry<_>>();

    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let reconstructed = data.reconstruct_solution(vector);
            let solution = general.compute_full_solution_with_reduced_solution(reconstructed);
            assert_eq!(solution, Solution::new(
                RB!(54),  // GLPK
                vec![
                    ("X1".to_string(), RB!(4)),
                    ("X2".to_string(), RB!(-1)),
                    ("X3".to_string(), RB!(6)),
                ],
            ));
        },
        _ => assert!(false),
    }
}

#[test]
fn testprob() {
    type T = Rational64;

    let mut general = to_general_form::<T>("testprob");
    let data = general.derive_matrix_data().ok().unwrap();
    let result = data.solve_relaxation::<Carry<_>>();

    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let reconstructed = data.reconstruct_solution(vector);
            let solution = general.compute_full_solution_with_reduced_solution(reconstructed);
            assert_eq!(solution, Solution::new(
                T::new(54, 1),  // GLPK
                vec![
                    ("X1".to_string(), T::new(4, 1)),
                    ("X2".to_string(), T::new(-1, 1)),
                    ("X3".to_string(), T::new(6, 1)),
                ],
            ));
        },
        _ => assert!(false),
    }
}
