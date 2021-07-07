use relp_num::{Rational16, Rational8, NonZeroFactorizable};
use relp_num::{R8, RB};

use crate::data::linear_algebra::matrix::{ColumnMajor, Order};
use crate::data::linear_algebra::vector::DenseVector;
use crate::data::linear_algebra::vector::test::TestVector;
use crate::data::linear_algebra::vector::Vector;
use crate::data::linear_program::elements::{Objective, RangedConstraintRelation, VariableType};
use crate::data::linear_program::general_form::{GeneralForm, Scalable, Scaling, Variable};
use crate::data::linear_program::general_form::presolve::scale::rational::{ColumnInfo, GeneralFormFactorization};
use std::ops::{MulAssign, DivAssign, Mul, Neg, Sub, Add, AddAssign, Div};
use num_traits::{One, Zero};
use std::hash::Hash;
use crate::data::linear_algebra::traits::{SparseElement, SparseComparator};

#[test]
fn test_scale_cost() {
    let mut general_form: GeneralForm<Rational8> = GeneralForm::new(
        Objective::Minimize,
        ColumnMajor::from_test_data::<_, _, u8>(&[
            vec![1, 2],
        ], 2),
        vec![
            RangedConstraintRelation::Equal,
        ],
        DenseVector::from_test_data::<u8>(vec![3]),
        vec![
            Variable {
                variable_type: VariableType::Continuous,
                cost: R8!(13 * 2),
                lower_bound: Some(R8!(3)),
                upper_bound: Some(R8!(5)),
                shift: R8!(0),
                flipped: false
            },
            Variable {
                variable_type: VariableType::Continuous,
                cost: R8!(13),
                lower_bound: Some(R8!(7)),
                upper_bound: Some(R8!(11)),
                shift: R8!(0),
                flipped: false
            },
        ],
        vec!["x".to_string(), "y".to_string()],
        R8!(16),
    );
    let scaling = general_form.scale();
    let expected_scaling = Scaling {
        cost_factor: R8!(1, 13),
        constraint_row_factors: vec![R8!(1)],
        constraint_column_factors: vec![R8!(1), R8!(1)],
    };
    assert_eq!(scaling, expected_scaling);
}

#[test]
fn test_scale() {
    let mut general_form: GeneralForm<Rational8> = GeneralForm::new(
        Objective::Minimize,
        ColumnMajor::from_test_data::<_, _, u8>(&[
            vec![11, 2],
            vec![4, 6],
            vec![7, 14],
            vec![0, 11],
        ], 2),
        vec![
            RangedConstraintRelation::Equal,
            RangedConstraintRelation::Less,
            RangedConstraintRelation::Greater,
            RangedConstraintRelation::Equal,
        ],
        DenseVector::from_test_data::<u8>(vec![3, 0, 21, 11]),
        vec![
            Variable {
                variable_type: VariableType::Continuous,
                cost: R8!(4),
                lower_bound: Some(R8!(0)),
                upper_bound: Some(R8!(6)),
                shift: R8!(0),
                flipped: false
            },
            Variable {
                variable_type: VariableType::Continuous,
                cost: R8!(11),
                lower_bound: Some(R8!(1)),
                upper_bound: Some(R8!(2)),
                shift: R8!(0),
                flipped: false
            },
        ],
        vec!["x".to_string(), "y".to_string()],
        R8!(16),
    );

    let mut factorization = general_form.factorize();

    let expected = GeneralFormFactorization {
        all_factors: vec![2, 3, 7, 11],
        b: vec![Some(vec![(3, 1)]), None, Some(vec![(3, 1), (7, 1)]), Some(vec![(11, 1)])],
        c: vec![Some(vec![(2, 2)]), Some(vec![(11, 1)])],
        bounds: vec![(None, Some(vec![(2, 1), (3, 1)])), (Some(vec![]), Some(vec![(2, 1)]))],
        A: vec![
            vec![(0, vec![(11, 1)]), (1, vec![(2, 2)]), (2, vec![(7, 1)])],
            vec![(0, vec![(2, 1)]), (1, vec![(2, 1), (3, 1)]), (2, vec![(2, 1), (7, 1)]), (3, vec![(11, 1)])],
        ],
    };
    assert_eq!(factorization, expected);

    let index = vec![
        vec![(0, 0), (1, 0)],
        vec![(0, 1), (1, 1)],
        vec![(0, 2), (1, 2)],
        vec![(1, 3)],
    ];
    let (column_info, gains) = factorization.solve_single_setup(&index);
    let expected_column_info = (
        Some(ColumnInfo {
            lowest_value: 0,
            lowest_count: 5,
            highest_value: 1,
        }),
        vec![
            ColumnInfo {
                lowest_value: 0,
                lowest_count: 4,
                highest_value: 1,
            },
            ColumnInfo {
                lowest_value: 0,
                lowest_count: 5,
                highest_value: 1,
            },
        ],
    );
    assert_eq!(column_info, expected_column_info);
    let expected_gains = (
        (0, 1),
        vec![
            (0, 4),
            (0, 0),
            (0, 0),
            (0, 1 * 2 + 1 * 4),
        ],
    );
    assert_eq!(gains, expected_gains);

    let row_increments = factorization.solve_single_rows(gains, column_info, &index);
    let expected_row_increments = (0, vec![0, 1, 1, 0]);
    assert_eq!(row_increments, expected_row_increments);

    let column_changes = factorization.solve_single_columns(&row_increments);
    let expected_column_changes = vec![1, 1];
    assert_eq!(column_changes, expected_column_changes);

    let factor = factorization.remove_factor_info();
    let expected_factor = 11;
    assert_eq!(factor, expected_factor);
    let expected_factorization = GeneralFormFactorization {
        all_factors: vec![2, 3, 7],
        b: vec![Some(vec![(3, 1)]), None, Some(vec![(3, 1), (7, 1)]), Some(vec![])],
        c: vec![Some(vec![(2, 2)]), Some(vec![])],
        bounds: vec![(None, Some(vec![(2, 1), (3, 1)])), (Some(vec![]), Some(vec![(2, 1)]))],
        A: vec![
            vec![(0, vec![]), (1, vec![(2, 2)]), (2, vec![(7, 1)])],
            vec![(0, vec![(2, 1)]), (1, vec![(2, 1), (3, 1)]), (2, vec![(2, 1), (7, 1)]), (3, vec![])],
        ],
    };
    assert_eq!(factorization, expected_factorization);

    let scaling = general_form.scale();
    let expected_scaling = Scaling {
        cost_factor: R8!(1),
        constraint_row_factors: vec![R8!(1), R8!(1, 2), R8!(1, 7), R8!(1, 11)],
        constraint_column_factors: vec![R8!(1), R8!(1)],
    };
    assert_eq!(scaling, expected_scaling);
    let expected_general_form = GeneralForm::new(
        Objective::Minimize,
        ColumnMajor::from_test_data::<_, _, u8>(&[
            vec![11, 2],
            vec![2, 3],
            vec![1, 2],
            vec![0, 1],
        ], 2),
        vec![
            RangedConstraintRelation::Equal,
            RangedConstraintRelation::Less,
            RangedConstraintRelation::Greater,
            RangedConstraintRelation::Equal,
        ],
        DenseVector::new(vec![R8!(1), R8!(0), R8!(1), R8!(1, 3)], 4),
        vec![
            Variable {
                variable_type: VariableType::Continuous,
                cost: R8!(4),
                lower_bound: Some(R8!(0)),
                upper_bound: Some(R8!(6)),
                shift: R8!(0),
                flipped: false,
            },
            Variable {
                variable_type: VariableType::Continuous,
                cost: R8!(11),
                lower_bound: Some(R8!(1)),
                upper_bound: Some(R8!(2)),
                shift: R8!(0),
                flipped: false,
            },
        ],
        vec!["x".to_string(), "y".to_string()],
        R8!(16),
    );
    assert_eq!(general_form, expected_general_form);
}

#[test]
fn test_solve_single_setup_without_b() {
    let factorization = GeneralFormFactorization::<Rational16> {
        all_factors: vec![2, 3, 7, 11],
        b: vec![None, None, None, None],
        c: vec![Some(vec![(11, 1)]), Some(vec![(2, 2)])],
        bounds: vec![(None, None), (None, None)],
        A: vec![
            vec![(0, vec![]), (1, vec![(2, 2)]), (2, vec![(7, 1)])],
            vec![(0, vec![(2, 1)]), (1, vec![(2, 1), (3, 1)]), (2, vec![(2, 1), (7, 1)]), (3, vec![(11, 1)])],
        ],
    };

    let index = vec![
        vec![(0, 0), (1, 0)],
        vec![(0, 1), (1, 1)],
        vec![(0, 2), (1, 2)],
        vec![(1, 3)],
    ];
    let result = factorization.solve_single_setup(&index);
    // For factor `11`
    let expected = (
        (
            None,
            vec![
                ColumnInfo {
                    lowest_value: 0,
                    lowest_count: 3,
                    highest_value: 1,
                },
                ColumnInfo {
                    lowest_value: 0,
                    lowest_count: 4,
                    highest_value: 1,
                },
            ],
        ),
        (
            (0, 1),
            vec![
                (0, 0),
                (0, 0),
                (0, 0),
                (0, 4),
            ],
        )
    );
    assert_eq!(result, expected);
}