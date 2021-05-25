use relp_num::R32;
use relp_num::Rational32;

use crate::data::linear_algebra::matrix::{ColumnMajor, Order};
use crate::data::linear_algebra::vector::DenseVector;
use crate::data::linear_algebra::vector::test::TestVector;
use crate::data::linear_program::elements::{LinearProgramType, Objective, VariableType};
use crate::data::linear_program::general_form::{GeneralForm, RangedConstraintRelation, Variable};
use crate::data::linear_program::solution::Solution;

type T = Rational32;

/// MIN Z = 211x1 + 223x2 + 227x3 - 229x4 + 233x5 + 0x6
/// subject to
/// 2x1 = 101
/// 3x1 + 5x2 <= 103
/// 7x1 + 11x2 + 13x3 >= 107
/// x2 >= -97/10
/// 17x1 + 19x2 + 23x3 + 29x5 + 31x6 = 109
/// x4 <= 131
/// x5 >= -30736/1885
/// x5 <= 123
/// x6 >= 5
/// and x1,x2,x3,x4,x5,x6 unrestricted in sign
#[test]
fn presolve() {
    let data = vec![
        // Column 3 should be removed because empty
        vec![2, 0, 0, 0, 0, 0], // Should be removed because simple bound
        vec![3, 5, 0, 0, 0, 0], // Should be removed because simple bound after removal of the row above
        vec![7, 11, 13, 0, 0, 0], // Should be removed because of fixed variable after the removal of above two
        vec![17, 19, 23, 0, 29, 31], // Removal by variable bounds
    ];
    let constraints = ColumnMajor::from_test_data::<T, T, _>(&data, 6);
    let b = DenseVector::from_test_data(vec![
        101,
        103,
        107,
        109,
    ]);
    let constraint_types = vec![
        RangedConstraintRelation::Equal,
        RangedConstraintRelation::Less,
        RangedConstraintRelation::Greater,
        RangedConstraintRelation::Equal,
    ];
    let column_info = vec![
        Variable {
            variable_type: VariableType::Continuous,
            cost: R32!(211),
            lower_bound: None,
            upper_bound: None,
            shift: R32!(0),
            flipped: false
        }, Variable {
            variable_type: VariableType::Continuous,
            cost: R32!(223),
            lower_bound: Some((R32!(103) - R32!(101) / R32!(2) * R32!(3)) / R32!(5)),
            upper_bound: None,
            shift: R32!(0),
            flipped: false
        }, Variable {
            variable_type: VariableType::Continuous,
            cost: R32!(227),
            lower_bound: None,
            upper_bound: None,
            shift: R32!(0),
            flipped: false
        }, Variable {
            variable_type: VariableType::Continuous,
            cost: R32!(-229),
            lower_bound: None,
            upper_bound: Some(R32!(131)),
            shift: R32!(0),
            flipped: false
        }, Variable {
            variable_type: VariableType::Continuous,
            cost: R32!(233),
            lower_bound: Some(R32!(-30736, 65 * 29)),
            upper_bound: Some(R32!(123)),
            shift: R32!(0),
            flipped: false
        }, Variable {
            variable_type: VariableType::Continuous,
            cost: R32!(0),
            lower_bound: Some(R32!(5)),
            upper_bound: None,
            shift: R32!(0),
            flipped: false
        },
    ];
    let variable_names = vec![
        "XONE".to_string(),
        "XTWO".to_string(),
        "XTHREE".to_string(),
        "XFOUR".to_string(),
        "XFIVE".to_string(),
        "XSIX".to_string(),
    ];
    let mut initial = GeneralForm::new(
        Objective::Minimize,
        constraints,
        constraint_types,
        b,
        column_info,
        variable_names,
        R32!(1),
    );
    debug_assert_eq!(
        initial.presolve(),
        Err(LinearProgramType::FiniteOptimum(Solution::new(
            R32!(1)
                + R32!(211 * 101, 2)
                + R32!(223 * -97, 10)
                + R32!(227 * -699, 65)
                + R32!(-229 * 131)
                + R32!(233 * -30736, 1885),
            vec![
                ("XONE".to_string(), R32!(101, 2)),
                ("XTWO".to_string(), (R32!(103) - R32!(101) / R32!(2) * R32!(3)) / R32!(5)),
                ("XTHREE".to_string(), (R32!(-3601, 5) + R32!(29 * 30736, 1885)) / R32!(23)),
                ("XFOUR".to_string(), R32!(131)),
                ("XFIVE".to_string(), R32!(-30736, 65 * 29)),
                ("XSIX".to_string(), R32!(5)),
            ],
        ))),
    );
}
