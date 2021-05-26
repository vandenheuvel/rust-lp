use std::collections::HashMap;

use relp_num::R32;

use crate::data::linear_algebra::matrix::{ColumnMajor, Order};
use crate::data::linear_algebra::vector::{DenseVector, Vector};
use crate::data::linear_program::elements::{Objective, VariableType};
use crate::data::linear_program::elements::BoundDirection;
use crate::data::linear_program::elements::LinearProgramType;
use crate::data::linear_program::general_form::{GeneralForm, RangedConstraintRelation, Variable};
use crate::data::linear_program::general_form::presolve::updates::Changes;
use crate::data::linear_program::general_form::RemovedVariable;

/// A simple problem that should not be presolved (further).
///
/// Many test problems that follow are variants of the problem below; they add a redundancy that
/// should be presolved.
#[test]
fn no_changes() {
    let initial = GeneralForm::new(
        Objective::Maximize,
        ColumnMajor::from_test_data(&[
                vec![R32!(1), R32!(1)],
                vec![R32!(4, 5), R32!(3, 2)],
                vec![R32!(3, 2), R32!(4, 5)],
            ],
            2,
        ),
        vec![
            RangedConstraintRelation::Less,
            RangedConstraintRelation::Less,
            RangedConstraintRelation::Less,
        ],
        DenseVector::new(vec![R32!(1), R32!(6, 5), R32!(6, 5)], 3),
        vec![Variable {
            variable_type: VariableType::Continuous,
            cost: R32!(1),
            lower_bound: Some(R32!(0)),
            upper_bound: Some(R32!(3, 4)),
            shift: R32!(0),
            flipped: false,
        }; 2],
        vec!["x".to_string(); 2],
        R32!(0),
    );
    assert_eq!(
        initial.compute_presolve_changes(),
        Ok(Changes {
            b: HashMap::default(),
            constraints: vec![],
            fixed_cost: R32!(0),
            bounds: HashMap::default(),
            removed_variables: vec![],
            constraints_marked_removed: vec![],
        }),
    );
}

/// Rules that get applied when indexes etc. are initialized.
mod initialization {
    use super::*;

    mod empty_variable {
        use super::*;

        /// Remove an empty column with a cost value.
        #[test]
        fn empty_column() {
            let initial = GeneralForm::new(
                Objective::Maximize,
                ColumnMajor::from_test_data(&[
                    vec![R32!(1), R32!(1), R32!(0)],
                    vec![R32!(4, 5), R32!(3, 2), R32!(0)],
                    vec![R32!(3, 2), R32!(4, 5), R32!(0)],
                ],
                                            3,
                ),
                vec![
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                ],
                DenseVector::new(vec![R32!(1), R32!(6, 5), R32!(6, 5)], 3),
                vec![Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(1),
                    lower_bound: Some(R32!(0)),
                    upper_bound: Some(R32!(3, 4)),
                    shift: R32!(0),
                    flipped: false,
                }; 3],
                vec!["x".to_string(); 3],
                R32!(0),
            );
            assert_eq!(
                initial.compute_presolve_changes(),
                Ok(Changes {
                    b: HashMap::default(),
                    constraints: vec![],
                    fixed_cost: R32!(3, 4),
                    bounds: HashMap::default(),
                    removed_variables: vec![(2, RemovedVariable::Solved(R32!(3, 4)))],
                    constraints_marked_removed: vec![],
                }),
            );
        }

        /// Remove an empty column that doesn't appear in the objective function.
        #[test]
        fn empty_column_no_cost() {
            let mut initial = GeneralForm::new(
                Objective::Maximize,
                ColumnMajor::from_test_data(&[
                    vec![R32!(1), R32!(1), R32!(0)],
                    vec![R32!(4, 5), R32!(3, 2), R32!(0)],
                    vec![R32!(3, 2), R32!(4, 5), R32!(0)],
                ],
                                            3,
                ),
                vec![
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                ],
                DenseVector::new(vec![R32!(1), R32!(6, 5), R32!(6, 5)], 3),
                vec![Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(1),
                    lower_bound: Some(R32!(0)),
                    upper_bound: Some(R32!(3, 4)),
                    shift: R32!(0),
                    flipped: false,
                }; 3],
                vec!["x".to_string(); 3],
                R32!(0),
            );
            initial.variables[2].cost = R32!(0);

            assert_eq!(
                initial.compute_presolve_changes(),
                Ok(Changes {
                    b: HashMap::default(),
                    constraints: vec![],
                    fixed_cost: R32!(0),
                    bounds: HashMap::default(),
                    removed_variables: vec![(2, RemovedVariable::Solved(R32!(3, 4)))],
                    constraints_marked_removed: vec![],
                }),
            );
        }

        /// Remove an empty column that doesn't appear in the objective function without any
        /// bounds.
        #[test]
        fn empty_column_no_cost_no_bound() {
            let mut initial = GeneralForm::new(
                Objective::Maximize,
                ColumnMajor::from_test_data(&[
                    vec![R32!(1), R32!(1), R32!(0)],
                    vec![R32!(4, 5), R32!(3, 2), R32!(0)],
                    vec![R32!(3, 2), R32!(4, 5), R32!(0)],
                ],
                                            3,
                ),
                vec![
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                ],
                DenseVector::new(vec![R32!(1), R32!(6, 5), R32!(6, 5)], 3),
                vec![Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(1),
                    lower_bound: None,
                    upper_bound: None,
                    shift: R32!(0),
                    flipped: false,
                }; 3],
                vec!["x".to_string(); 3],
                R32!(0),
            );
            initial.variables[2].cost = R32!(0);

            assert_eq!(
                initial.compute_presolve_changes(),
                Ok(Changes {
                    b: HashMap::default(),
                    constraints: vec![],
                    fixed_cost: R32!(0),
                    bounds: HashMap::default(),
                    removed_variables: vec![(2, RemovedVariable::Solved(R32!(0)))],
                    constraints_marked_removed: vec![],
                }),
            );
        }

    }


    /// Removal of empty constraint.
    mod empty_constraint {
        use super::*;

        /// Rhs is nonzero and constraint type is equality.
        #[test]
        fn infeasible_equality_bound() {
            let initial = GeneralForm::new(
                Objective::Maximize,
                ColumnMajor::from_test_data(&[
                    vec![R32!(1), R32!(1)],
                    vec![R32!(4, 5), R32!(3, 2)],
                    vec![R32!(3, 2), R32!(4, 5)],
                    vec![R32!(0), R32!(0)],
                ],
                                            2,
                ),
                vec![
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Equal,
                ],
                DenseVector::new(vec![R32!(1), R32!(6, 5), R32!(6, 5), R32!(123)], 4),
                vec![Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(1),
                    lower_bound: Some(R32!(0)),
                    upper_bound: Some(R32!(3, 4)),
                    shift: R32!(0),
                    flipped: false,
                }; 2],
                vec!["x".to_string(); 2],
                R32!(0),
            );
            assert_eq!(
                initial.compute_presolve_changes(),
                Err(LinearProgramType::Infeasible),
            );
        }

        /// Rhs is nonzero and constraint type is in the wrong direction.
        #[test]
        fn infeasible_inequality_bound() {
            let initial = GeneralForm::new(
                Objective::Maximize,
                ColumnMajor::from_test_data(&[
                    vec![R32!(1), R32!(1)],
                    vec![R32!(4, 5), R32!(3, 2)],
                    vec![R32!(3, 2), R32!(4, 5)],
                    vec![R32!(0), R32!(0)],
                ],
                                            2,
                ),
                vec![
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                ],
                DenseVector::new(vec![R32!(1), R32!(6, 5), R32!(6, 5), R32!(-123)], 4),
                vec![Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(1),
                    lower_bound: None,
                    upper_bound: Some(R32!(3, 4)),
                    shift: R32!(0),
                    flipped: false,
                }; 2],
                vec!["x".to_string(); 2],
                R32!(0),
            );
            assert_eq!(
                initial.compute_presolve_changes(),
                Err(LinearProgramType::Infeasible),
            );
        }

        /// Rhs is zero, constraint can be removed.
        #[test]
        fn feasible_multiple() {
            let initial = GeneralForm::new(
                Objective::Maximize,
                ColumnMajor::from_test_data(&[
                    vec![R32!(1), R32!(1)],
                    vec![R32!(4, 5), R32!(3, 2)],
                    vec![R32!(0), R32!(0)],
                    vec![R32!(3, 2), R32!(4, 5)],
                    vec![R32!(0), R32!(0)],
                    vec![R32!(0), R32!(0)],
                ],
                                            2,
                ),
                vec![
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Equal,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Greater,
                    RangedConstraintRelation::Less,
                ],
                DenseVector::new(vec![
                    R32!(1),
                    R32!(6, 5),
                    R32!(0),
                    R32!(6, 5),
                    R32!(-1234),
                    R32!(0),
                ], 6),
                vec![Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(1),
                    lower_bound: Some(R32!(0)),
                    upper_bound: None,
                    shift: R32!(0),
                    flipped: false,
                }; 2],
                vec!["x".to_string(); 2],
                R32!(0),
            );
            assert_eq!(
                initial.compute_presolve_changes(),
                Ok(Changes {
                    b: HashMap::default(),
                    constraints: vec![],
                    fixed_cost: R32!(0),
                    bounds: HashMap::default(),
                    removed_variables: vec![],
                    constraints_marked_removed: vec![2, 4, 5],
                }),
            );
        }

        /// Multiple can be removed, but one can't.
        #[test]
        fn infeasible_multiple() {
            let initial = GeneralForm::new(
                Objective::Maximize,
                ColumnMajor::from_test_data(&[
                    vec![R32!(1), R32!(1)],
                    vec![R32!(4, 5), R32!(3, 2)],
                    vec![R32!(0), R32!(0)],
                    vec![R32!(3, 2), R32!(4, 5)],
                    vec![R32!(0), R32!(0)],
                    vec![R32!(0), R32!(0)],
                    vec![R32!(0), R32!(0)],
                ],
                                            2,
                ),
                vec![
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Equal,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Greater,
                    RangedConstraintRelation::Less,
                ],
                DenseVector::new(vec![
                    R32!(1),
                    R32!(6, 5),
                    R32!(0),
                    R32!(6, 5),
                    R32!(-5678),
                    R32!(-1234),
                    R32!(0),
                ], 7),
                vec![Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(1),
                    lower_bound: Some(R32!(0)),
                    upper_bound: None,
                    shift: R32!(0),
                    flipped: false,
                }; 2],
                vec!["x".to_string(); 2],
                R32!(0),
            );
            assert_eq!(
                initial.compute_presolve_changes(),
                Err(LinearProgramType::Infeasible),
            );
        }
    }
}

/// Removal of fixed variable.
mod fixed_variable {
    use super::*;

    /// Shifting bound
    #[test]
    fn shift_bounds() {
        let initial = {
            let mut initial = GeneralForm::new(
                Objective::Maximize,
                    ColumnMajor::from_test_data(&[
                        vec![R32!(1), R32!(1), R32!(1, 1000)],
                        vec![R32!(4, 5), R32!(3, 2), R32!(0)],
                        vec![R32!(3, 2), R32!(4, 5), R32!(0)],
                    ],
                    3,
                ),
                vec![
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                ],
                DenseVector::new(vec![R32!(1), R32!(6, 5), R32!(6, 5)], 3),
                vec![Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(1),
                    lower_bound: Some(R32!(0)),
                    upper_bound: Some(R32!(3, 4)),
                    shift: R32!(0),
                    flipped: false,
                }; 3],
                vec!["x".to_string(); 3],
                R32!(0),
            );

            initial.variables[2].lower_bound = Some(R32!(1));
            initial.variables[2].upper_bound = Some(R32!(1));
            initial.variables[2].cost = R32!(0);

            initial
        };
        assert_eq!(
            initial.compute_presolve_changes(),
            Ok(Changes {
                b: {
                    let mut set = HashMap::new();
                    set.insert(0, R32!(999, 1000));
                    set
                },
                constraints: vec![],
                fixed_cost: R32!(0),
                bounds: HashMap::default(),
                removed_variables: vec![(2, RemovedVariable::Solved(R32!(1)))],
                constraints_marked_removed: vec![],
            }),
        );
    }

    /// Shifting global fixed costs
    #[test]
    fn shift_costs() {
        let initial = {
            let mut initial = GeneralForm::new(
                Objective::Maximize,
                ColumnMajor::from_test_data(&[
                        vec![R32!(1), R32!(1), R32!(0)],
                        vec![R32!(4, 5), R32!(3, 2), R32!(0)],
                        vec![R32!(3, 2), R32!(4, 5), R32!(0)],
                    ],
                    3,
                ),
                vec![
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                ],
                DenseVector::new(vec![R32!(1), R32!(6, 5), R32!(6, 5)], 3),
                vec![Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(1),
                    lower_bound: Some(R32!(0)),
                    upper_bound: Some(R32!(3, 4)),
                    shift: R32!(0),
                    flipped: false,
                }; 3],
                vec!["x".to_string(); 3],
                R32!(0),
            );

            initial.variables[2].lower_bound = Some(R32!(3));
            initial.variables[2].upper_bound = Some(R32!(3));
            initial.variables[2].cost = R32!(10);

            initial
        };
        assert_eq!(
            initial.compute_presolve_changes(),
            Ok(Changes {
                b: HashMap::default(),
                constraints: vec![],
                fixed_cost: R32!(30),
                bounds: HashMap::default(),
                removed_variables: vec![(2, RemovedVariable::Solved(R32!(3)))],
                constraints_marked_removed: vec![],
            }),
        );
    }

    /// Triggering an empty row to be presolved.
    ///
    /// Middle variable removed, first constraint removed.
    #[test]
    fn trigger_empty_row() {
        let initial = {
            let mut initial = GeneralForm::new(
                Objective::Maximize,
                ColumnMajor::from_test_data(&[
                        vec![R32!(0), R32!(-51), R32!(0)],
                        vec![R32!(1), R32!(1, 100), R32!(1)],
                        vec![R32!(4, 5), R32!(0), R32!(3, 2)],
                        vec![R32!(3, 2), R32!(0), R32!(4, 5)],
                    ],
                    3,
                ),
                vec![
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                ],
                DenseVector::new(vec![R32!(-20), R32!(1), R32!(6, 5), R32!(6, 5)], 4),
                vec![Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(1),
                    lower_bound: Some(R32!(0)),
                    upper_bound: Some(R32!(3, 4)),
                    shift: R32!(0),
                    flipped: false,
                }; 3],
                vec!["x".to_string(); 3],
                R32!(0),
            );

            initial.variables[1].lower_bound = Some(R32!(1, 2));
            initial.variables[1].upper_bound = Some(R32!(1, 2));
            initial.variables[1].cost = R32!(15);

            initial
        };
        assert_eq!(
            initial.compute_presolve_changes(),
            Ok(Changes {
                b: {
                    let mut set = HashMap::new();
                    set.insert(1, R32!(199, 200));
                    set
                },
                constraints: vec![],
                fixed_cost: R32!(15, 2),
                bounds: HashMap::default(),
                removed_variables: vec![(1, RemovedVariable::Solved(R32!(1, 2)))],
                constraints_marked_removed: vec![0],
            }),
        );
    }


    /// Triggering a bound row to be removed.
    ///
    /// Middle variable removed, first constraint read as bound.
    #[test]
    fn trigger_bound_rule() {
        let initial = {
            let mut initial = GeneralForm::new(
                Objective::Maximize,
                ColumnMajor::from_test_data(&[
                        vec![R32!(0), R32!(1, 2), R32!(1)],
                        vec![R32!(1), R32!(0), R32!(1)],
                        vec![R32!(4, 5), R32!(1, 20), R32!(3, 2)],
                        vec![R32!(3, 2), R32!(0), R32!(4, 5)],
                    ],
                    3,
                ),
                vec![
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                ],
                DenseVector::new(vec![R32!(1), R32!(1), R32!(6, 5), R32!(6, 5)], 4),
                vec![Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(-12),
                    lower_bound: Some(R32!(0)),
                    upper_bound: None,
                    shift: R32!(0),
                    flipped: false,
                }; 3],
                vec!["x".to_string(); 3],
                R32!(0),
            );

            initial.variables[1].lower_bound = Some(R32!(1, 2));
            initial.variables[1].upper_bound = Some(R32!(1, 2));
            initial.variables[1].cost = R32!(-11);

            initial
        };
        assert_eq!(
            initial.compute_presolve_changes(),
            Ok(Changes {
                b: {
                    let mut set = HashMap::new();
                    set.insert(2, R32!(47, 40));
                    set
                },
                constraints: vec![],
                fixed_cost: R32!(-11, 2),
                bounds: {
                    let mut set = HashMap::new();
                    set.insert((2, BoundDirection::Upper), R32!(3, 4));
                    set
                },
                removed_variables: vec![(1, RemovedVariable::Solved(R32!(1, 2)))],
                constraints_marked_removed: vec![0],
            }),
        );
    }
}

/// Removal of bound constraint.
mod bound_constraint {
    use super::*;

    /// Read inequality bound (i.e. one side).
    #[test]
    fn feasible_inequality() {
        let initial = GeneralForm::new(
            Objective::Maximize,
            ColumnMajor::from_test_data(&[
                vec![R32!(1), R32!(1)],
                vec![R32!(4, 5), R32!(3, 2)],
                vec![R32!(0), R32!(1)],
                vec![R32!(3, 2), R32!(4, 5)],
            ],
                                        2,
            ),
            vec![
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Less,
            ],
            DenseVector::new(vec![R32!(1), R32!(6, 5), R32!(3, 4), R32!(6, 5)], 4),
            vec![Variable {
                variable_type: VariableType::Continuous,
                cost: R32!(-2),
                lower_bound: None,
                upper_bound: None,
                shift: R32!(14),
                flipped: false,
            }; 2],
            vec!["x".to_string(); 2],
            R32!(12),
        );
        assert_eq!(
            initial.compute_presolve_changes(),
            Ok(Changes {
                b: HashMap::default(),
                constraints: vec![],
                fixed_cost: R32!(0),
                bounds: {
                    let mut map = HashMap::new();
                    map.insert((1, BoundDirection::Upper), R32!(3, 4));
                    map
                },
                removed_variables: vec![],
                constraints_marked_removed: vec![2],
            }),
        );
    }

    /// Read inequality bound (i.e. one side).
    #[test]
    fn feasible_inequality_negative_sign() {
        let initial = GeneralForm::new(
            Objective::Maximize,
            ColumnMajor::from_test_data(&[
                    vec![R32!(1), R32!(1)],
                    vec![R32!(4, 5), R32!(3, 2)],
                    vec![R32!(0), R32!(-1)],
                    vec![R32!(3, 2), R32!(4, 5)],
                ],
                2,
            ),
            vec![
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Less,
            ],
            DenseVector::new(vec![R32!(1), R32!(6, 5), -R32!(3, 4), R32!(6, 5)], 4),
            vec![Variable {
                variable_type: VariableType::Continuous,
                cost: R32!(-2),
                lower_bound: None,
                upper_bound: None,
                shift: R32!(14),
                flipped: false,
            }; 2],
            vec!["x".to_string(); 2],
            R32!(12),
        );
        assert_eq!(
            initial.compute_presolve_changes(),
            Ok(Changes {
                b: HashMap::default(),
                constraints: vec![],
                fixed_cost: R32!(0),
                bounds: {
                    let mut map = HashMap::new();
                    map.insert((0, BoundDirection::Upper), R32!(3, 32));
                    map.insert((1, BoundDirection::Lower), R32!(3, 4));
                    map
                },
                removed_variables: vec![],
                constraints_marked_removed: vec![2],
            }),
        );
    }

    /// Read equality bound, not requiring substitution.
    #[test]
    fn feasible_equality_unsolved_no_substitution() {
        let initial = GeneralForm::new(
            Objective::Maximize,
            ColumnMajor::from_test_data(&[
                    vec![R32!(1),    R32!(1),    R32!(0)],
                    vec![R32!(4, 5), R32!(3, 2), R32!(0)],
                    vec![R32!(0),    R32!(0),    R32!(1)],
                    vec![R32!(3, 2), R32!(4, 5), R32!(0)],
                ],
                3,
            ),
            vec![
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Equal,
                RangedConstraintRelation::Less,
            ],
            DenseVector::new(vec![R32!(1), R32!(6, 5), R32!(10), R32!(6, 5)], 4),
            vec![Variable {
                variable_type: VariableType::Continuous,
                cost: R32!(2),
                lower_bound: None,
                upper_bound: None,
                shift: R32!(14),
                flipped: false,
            }; 3],
            vec!["x".to_string(); 3],
            R32!(12),
        );
        assert_eq!(
            initial.compute_presolve_changes(),
            Ok(Changes {
                b: HashMap::default(),
                constraints: vec![],
                fixed_cost: R32!(20),
                bounds: HashMap::default(),
                removed_variables: vec![(2, RemovedVariable::Solved(R32!(10)))],
                constraints_marked_removed: vec![2],
            }),
        );
    }

    /// Read equality bound, requiring substitution.
    #[test]
    fn feasible_equality_unsolved_substitution() {
        let initial = GeneralForm::new(
            Objective::Maximize,
            ColumnMajor::from_test_data(&[
                    vec![R32!(1),    R32!(1),    R32!(0)     ],
                    vec![R32!(4, 5), R32!(3, 2), R32!(0)     ],
                    vec![R32!(0),    R32!(0),    R32!(-1)    ],
                    vec![R32!(3, 2), R32!(4, 5), R32!(1, 100)],
                ],
                3,
            ),
            vec![
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Equal,
                RangedConstraintRelation::Less,
            ],
            DenseVector::new(vec![R32!(1), R32!(6, 5), -R32!(3), R32!(6, 5)], 4),
            vec![Variable {
                variable_type: VariableType::Continuous,
                cost: R32!(2),
                lower_bound: None,
                upper_bound: None,
                shift: R32!(14),
                flipped: false,
            }; 3],
            vec!["x".to_string(); 3],
            R32!(12),
        );
        assert_eq!(
            initial.compute_presolve_changes(),
            Ok(Changes {
                b: {
                    let mut map = HashMap::new();
                    map.insert(3,R32!(117, 100));
                    map
                },
                constraints: vec![],
                fixed_cost: R32!(6),
                bounds: HashMap::default(),
                removed_variables: vec![(2, RemovedVariable::Solved(R32!(3)))],
                constraints_marked_removed: vec![2],
            }),
        );
    }

    /// Read inequality bound (i.e. one side).
    #[test]
    fn feasible_equality_solved() {
        let initial = GeneralForm::new(
            Objective::Maximize,
            ColumnMajor::from_test_data(&[
                    vec![R32!(1), R32!(1)],
                    vec![R32!(4, 5), R32!(3, 2)],
                    vec![R32!(1), R32!(0)],
                    vec![R32!(3, 2), R32!(4, 5)],
                ],
                2,
            ),
            vec![
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Equal,
                RangedConstraintRelation::Less,
            ],
            DenseVector::new(vec![R32!(1), R32!(6, 5), R32!(3, 4), R32!(6, 5)], 4),
            vec![Variable {
                variable_type: VariableType::Continuous,
                cost: R32!(2),
                lower_bound: None,
                upper_bound: None,
                shift: R32!(14),
                flipped: false,
            }; 2],
            vec!["x".to_string(); 2],
            R32!(12),
        );
        assert_eq!(
            initial.compute_presolve_changes(),
            Ok(Changes {
                b: HashMap::default(),
                constraints: vec![],
                fixed_cost: (R32!(3, 4) + R32!(3, 32)) * R32!(2),
                bounds: HashMap::default(),
                removed_variables: vec![
                    (0, RemovedVariable::Solved(R32!(3, 4))),
                    (1, RemovedVariable::Solved(R32!(3, 32))),
                ],
                constraints_marked_removed: vec![0, 1, 2, 3],
            }),
        );
    }

    /// Read bound leading to variable solved to an arbitrary value.
    #[test]
    fn feasible_equality_unsolved_independent_column() {
        let initial = {
            let mut initial = GeneralForm::new(
                Objective::Minimize,
                ColumnMajor::from_test_data(&[
                    vec![R32!(1),    R32!(1),    R32!(0) ],
                    vec![R32!(4, 5), R32!(3, 2), R32!(0) ],
                    vec![R32!(0),    R32!(0),    R32!(1)],
                    vec![R32!(3, 2), R32!(4, 5), R32!(0) ],
                ],
                                            3,
                ),
                vec![
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Greater,
                    RangedConstraintRelation::Less,
                ],
                DenseVector::new(vec![R32!(1), R32!(6, 5), -R32!(3), R32!(6, 5)], 4),
                vec![Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(-3),
                    lower_bound: None,
                    upper_bound: None,
                    shift: R32!(-1),
                    flipped: true,
                }; 3],
                vec!["x".to_string(); 3],
                R32!(-12),
            );
            initial.variables[2].cost = R32!(0);
            initial
        };
        assert_eq!(
            initial.compute_presolve_changes(),
            Ok(Changes {
                b: HashMap::default(),
                constraints: vec![],
                fixed_cost: R32!(0),
                bounds: HashMap::default(),
                removed_variables: vec![(2, RemovedVariable::Solved(R32!(-3)))],
                constraints_marked_removed: vec![2],
            }),
        );
    }

    /// Read bound leading to variable solved as a function of others.
    #[test]
    fn feasible_equality_unsolved_slack() {
        let initial = {
            let mut initial = GeneralForm::new(
                Objective::Minimize,
                ColumnMajor::from_test_data(&[
                        vec![R32!(1),    R32!(1),    R32!(1)],
                        vec![R32!(4, 5), R32!(3, 2), R32!(0)],
                        vec![R32!(0),    R32!(0),    R32!(1)],
                        vec![R32!(3, 2), R32!(4, 5), R32!(0)],
                    ],
                    3,
                ),
                vec![
                    RangedConstraintRelation::Equal,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Greater,
                    RangedConstraintRelation::Less,
                ],
                DenseVector::new(vec![R32!(2), R32!(6, 5), R32!(1), R32!(6, 5)], 4),
                vec![Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(-3),
                    lower_bound: None,
                    upper_bound: None,
                    shift: R32!(-1),
                    flipped: true,
                }; 3],
                vec!["x".to_string(); 3],
                R32!(-12),
            );
            initial.variables[2].cost = R32!(0);
            initial
        };
        assert_eq!(
            initial.compute_presolve_changes(),
            Ok(Changes {
                b: {
                    let mut map = HashMap::new();
                    map.insert(0, R32!(1));
                    map
                },
                constraints: vec![(0, RangedConstraintRelation::Less)],
                fixed_cost: R32!(0),
                bounds: HashMap::default(),
                removed_variables: vec![(2, RemovedVariable::FunctionOfOthers {
                    constant: R32!(2),
                    coefficients: vec![(0, R32!(1)), (1, R32!(1))],
                })],
                constraints_marked_removed: vec![2],
            }),
        );
    }
}

/// Removal of a slack variable.
mod slack_variable {
    use super::*;

    /// Removing a simple slack variable.
    #[test]
    fn remove_slack() {
        let initial = {
            let mut initial = GeneralForm::new(
                Objective::Maximize,
                ColumnMajor::from_test_data(&[
                        vec![R32!(1), R32!(1), R32!(1)],
                        vec![R32!(4, 5), R32!(3, 2), R32!(0)],
                        vec![R32!(3, 2), R32!(4, 5), R32!(0)],
                    ],
                    3,
                ),
                vec![
                    RangedConstraintRelation::Equal,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                ],
                DenseVector::new(vec![R32!(1), R32!(6, 5), R32!(6, 5)], 3),
                vec![Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(1),
                    lower_bound: Some(R32!(0)),
                    upper_bound: Some(R32!(3, 4)),
                    shift: R32!(0),
                    flipped: false,
                }; 3],
                vec!["x".to_string(); 3],
                R32!(0),
            );
            initial.variables[2].cost = R32!(0);
            initial.variables[2].upper_bound = None;
            initial
        };
        assert_eq!(
            initial.compute_presolve_changes(),
            Ok(Changes {
                b: HashMap::default(),
                constraints: vec![(0, RangedConstraintRelation::Less)],
                fixed_cost: R32!(0),
                bounds: HashMap::default(),
                removed_variables: vec![(2, RemovedVariable::FunctionOfOthers {
                    constant: R32!(1),
                    coefficients: vec![(0, R32!(1)), (1, R32!(1))]
                })],
                constraints_marked_removed: vec![],
            }),
        );
    }

    /// Don't remove a range slack variable.
    #[test]
    fn remove_range_slack() {
        let initial = {
            let mut initial = GeneralForm::new(
                Objective::Maximize,
                ColumnMajor::from_test_data(&[
                    vec![R32!(1), R32!(1), R32!(1)],
                    vec![R32!(4, 5), R32!(3, 2), R32!(0)],
                    vec![R32!(3, 2), R32!(4, 5), R32!(0)],
                ],
                    3,
                ),
                vec![
                    RangedConstraintRelation::Equal,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                ],
                DenseVector::new(vec![R32!(1), R32!(6, 5), R32!(6, 5)], 3),
                vec![Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(1),
                    lower_bound: Some(R32!(0)),
                    upper_bound: Some(R32!(3, 4)),
                    shift: R32!(0),
                    flipped: false,
                }; 3],
                vec!["x".to_string(); 3],
                R32!(0),
            );
            initial.variables[2].cost = R32!(0);
            initial
        };
        assert_eq!(
            initial.compute_presolve_changes(),
            Ok(Changes {
                b: HashMap::default(),
                constraints: vec![(0, RangedConstraintRelation::Range(R32!(3, 4)))],
                fixed_cost: R32!(0),
                bounds: HashMap::default(),
                removed_variables: vec![(2, RemovedVariable::FunctionOfOthers {
                    constant: R32!(1),
                    coefficients: vec![(0, R32!(1)), (1, R32!(1))],
                })],
                constraints_marked_removed: vec![],
            }),
        );
    }

    /// Remove the entire constraint.
    #[test]
    fn free_variable() {
        let initial = {
            let mut initial = GeneralForm::new(
                Objective::Maximize,
                ColumnMajor::from_test_data(&[
                        vec![R32!(1), R32!(1), R32!(1)],
                        vec![R32!(4, 5), R32!(3, 2), R32!(0)],
                        vec![R32!(3, 2), R32!(4, 5), R32!(0)],
                    ],
                    3,
                ),
                vec![
                    RangedConstraintRelation::Equal,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                ],
                DenseVector::new(vec![R32!(1), R32!(6, 5), R32!(6, 5)], 3),
                vec![Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(1),
                    lower_bound: Some(R32!(0)),
                    upper_bound: Some(R32!(3, 4)),
                    shift: R32!(0),
                    flipped: false,
                }; 3],
                vec!["x".to_string(); 3],
                R32!(0),
            );
            initial.variables[2].cost = R32!(0);
            initial.variables[2].lower_bound = None;
            initial.variables[2].upper_bound = None;
            initial
        };
        assert_eq!(
            initial.compute_presolve_changes(),
            Ok(Changes {
                b: HashMap::default(),
                constraints: vec![],
                fixed_cost: R32!(0),
                bounds: HashMap::default(),
                removed_variables: vec![(2, RemovedVariable::FunctionOfOthers {
                    constant: R32!(1),
                    coefficients: vec![(0, R32!(1)), (1, R32!(1))],
                })],
                constraints_marked_removed: vec![0],
            }),
        );
    }

    /// Derive variable bounds that solve the problem.
    #[test]
    fn range_slack_for_range_bound() {
        let initial = {
            let mut initial = GeneralForm::new(
                Objective::Maximize,
                ColumnMajor::from_test_data(&[
                        vec![R32!(1), R32!(1), R32!(1)],
                        vec![R32!(4, 5), R32!(3, 2), R32!(0)],
                        vec![R32!(3, 2), R32!(4, 5), R32!(0)],
                    ],
                    3,
                ),
                vec![
                    RangedConstraintRelation::Range(R32!(1, 2)),
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                ],
                DenseVector::new(vec![R32!(1), R32!(6, 5), R32!(6, 5)], 3),
                vec![Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(14),
                    lower_bound: None,
                    upper_bound: None,
                    shift: R32!(0),
                    flipped: false,
                }; 3],
                vec!["x".to_string(); 3],
                R32!(0),
            );
            initial.variables[2].lower_bound = Some(R32!(0));
            initial.variables[2].upper_bound = Some(R32!(1, 2));
            initial.variables[2].cost = R32!(0);
            initial
        };

        assert_eq!(
            initial.compute_presolve_changes(),
            Ok(Changes {
                b: HashMap::default(),
                constraints: vec![(0, RangedConstraintRelation::Range(R32!(1)))],
                fixed_cost: R32!(0),
                bounds: HashMap::default(),
                removed_variables: vec![(2, RemovedVariable::FunctionOfOthers {
                    constant: R32!(1),
                    coefficients: vec![(0, R32!(1)), (1, R32!(1))],
                })],
                constraints_marked_removed: vec![],
            }),
        );
    }

    /// Remove a range slack from an inequality constraint and trigger the variable bound rule.
    #[test]
    fn range_slack_inequality_bound() {
        let mut initial = GeneralForm::new(
            Objective::Maximize,
            ColumnMajor::from_test_data(&[
                    vec![R32!(1), R32!(0), R32!(1)],
                    vec![R32!(1), R32!(1), R32!(0)],
                    vec![R32!(4, 5), R32!(3, 2), R32!(0)],
                    vec![R32!(3, 2), R32!(4, 5), R32!(0)],
                ],
                3,
            ),
            vec![
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Less,
            ],
            DenseVector::new(vec![R32!(1), R32!(1), R32!(6, 5), R32!(6, 5)], 4),
            vec![Variable {
                variable_type: VariableType::Continuous,
                cost: R32!(1),
                lower_bound: None,
                upper_bound: None,
                shift: R32!(0),
                flipped: false,
            }; 3],
            vec!["x".to_string(); 3],
            R32!(0),
        );
        initial.variables[2].lower_bound = Some(R32!(0));
        initial.variables[2].upper_bound = Some(R32!(1));
        initial.variables[2].cost = R32!(0);

        assert_eq!(
            initial.compute_presolve_changes(),
            Ok(Changes {
                b: HashMap::default(),
                constraints: vec![],
                fixed_cost: R32!(0),
                bounds: {
                    let mut map = HashMap::new();
                    map.insert((0, BoundDirection::Upper), R32!(1));
                    map
                },
                removed_variables: vec![(2, RemovedVariable::Solved(R32!(0)))],
                constraints_marked_removed: vec![0],
            }),
        );
    }
}

/// Domain propagation, activity bounds
mod domain_propagation {
    use super::*;

    /// Remove a constraint that is always satisfied due to the variable bounds.
    #[test]
    fn remove_redundant_constraint() {
        let initial = GeneralForm::new(
            Objective::Maximize,
            ColumnMajor::from_test_data(&[
                    vec![R32!(1), R32!(1)],
                    vec![R32!(4, 5), R32!(3, 2)],
                    vec![R32!(3, 2), R32!(4, 5)],
                ],
                2,
            ),
            vec![
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Less,
            ],
            DenseVector::new(vec![R32!(3, 2), R32!(6, 5), R32!(6, 5)], 3),
            vec![Variable {
                variable_type: VariableType::Continuous,
                cost: R32!(1),
                lower_bound: Some(R32!(0)),
                upper_bound: Some(R32!(3, 4)),
                shift: R32!(0),
                flipped: false,
            }; 2],
            vec!["x".to_string(); 2],
            R32!(0),
        );
        assert_eq!(
            initial.compute_presolve_changes(),
            Ok(Changes {
                b: HashMap::default(),
                constraints: vec![],
                fixed_cost: R32!(0),
                bounds: HashMap::default(),
                removed_variables: vec![],
                constraints_marked_removed: vec![0],
            }),
        );
    }

    /// Derive a variable / activity bound that shows the problem to be infeasible.
    #[test]
    fn trigger_fixed_unfeasible() {
        let initial = {
            let mut initial = GeneralForm::new(
                Objective::Maximize,
                ColumnMajor::from_test_data(&[
                        vec![R32!(1), R32!(1)],
                        vec![R32!(4, 5), R32!(3, 2)],
                        vec![R32!(3, 2), R32!(4, 5)],
                    ],
                    2,
                ),
                vec![
                    RangedConstraintRelation::Equal,
                    RangedConstraintRelation::Less,
                    RangedConstraintRelation::Less,
                ],
                DenseVector::new(vec![R32!(1), R32!(6, 5), R32!(6, 5)], 3),
                vec![Variable {
                    variable_type: VariableType::Continuous,
                    cost: R32!(14),
                    lower_bound: None,
                    upper_bound: None,
                    shift: R32!(0),
                    flipped: false,
                }; 2],
                vec!["x".to_string(); 2],
                R32!(0),
            );

            initial.variables[0].lower_bound = Some(R32!(1, 4));
            initial.variables[1].lower_bound = Some(R32!(3, 4));
            initial
        };

        assert_eq!(initial.compute_presolve_changes(), Err(LinearProgramType::Infeasible));
    }

    /// Derive variable bounds that solve the problem.
    #[test]
    fn trigger_substitution_feasible() {
        let initial = GeneralForm::new(
            Objective::Maximize,
            ColumnMajor::from_test_data(&[
                vec![R32!(1), R32!(1)],
                vec![R32!(4, 5), R32!(3, 2)],
                vec![R32!(3, 2), R32!(4, 5)],
            ],
                                        2,
            ),
            vec![
                RangedConstraintRelation::Equal,
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Less,
            ],
            DenseVector::new(vec![R32!(1), R32!(6, 5), R32!(6, 5)], 3),
            vec![Variable {
                variable_type: VariableType::Continuous,
                cost: R32!(14),
                lower_bound: Some(R32!(1, 2)),
                upper_bound: None,
                shift: R32!(0),
                flipped: false,
            }; 2],
            vec!["x".to_string(); 2],
            R32!(0),
        );

        assert_eq!(
            initial.compute_presolve_changes(),
            Ok(Changes {
                b: HashMap::default(),
                constraints: vec![],
                fixed_cost: R32!(14) * (R32!(1, 2) + R32!(1, 2)),
                bounds: HashMap::default(),
                removed_variables: vec![
                    (0, RemovedVariable::Solved(R32!(1, 2))),
                    (1, RemovedVariable::Solved(R32!(1, 2))),
                ],
                constraints_marked_removed: vec![0, 1, 2],
            }),
        );
    }

    /// Derive variable bounds that make a constraint unsatisfiable.
    #[test]
    fn trigger_substitution_infeasible() {
        let initial = GeneralForm::new(
            Objective::Maximize,
            ColumnMajor::from_test_data(&[
                vec![R32!(1), R32!(1)],
                vec![R32!(4, 5), R32!(3, 2)],
                vec![R32!(3, 2), R32!(4, 5)],
            ],
                                        2,
            ),
            vec![
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Less,
            ],
            DenseVector::new(vec![R32!(1), R32!(6, 5), R32!(6, 5)], 3),
            vec![Variable {
                variable_type: VariableType::Continuous,
                cost: R32!(14),
                lower_bound: Some(R32!(2, 3)),
                upper_bound: None,
                shift: R32!(0),
                flipped: false,
            }; 2],
            vec!["x".to_string(); 2],
            R32!(0),
        );

        assert_eq!(
            initial.compute_presolve_changes(),
            Err(LinearProgramType::Infeasible),
        );
    }

    /// Remove all constraints, because they are all satisfied by activity bounds.
    #[test]
    fn all_constraints_redundant() {
        let initial = GeneralForm::new(
            Objective::Maximize,
            ColumnMajor::from_test_data(&[
                vec![R32!(1), R32!(1)],
                vec![R32!(4, 5), R32!(3, 2)],
                vec![R32!(3, 2), R32!(4, 5)],
            ],
                                        2,
            ),
            vec![
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Less,
                RangedConstraintRelation::Less,
            ],
            DenseVector::new(vec![R32!(3, 2), R32!(6, 5), R32!(6, 5)], 3),
            vec![Variable {
                variable_type: VariableType::Continuous,
                cost: R32!(1),
                lower_bound: Some(R32!(0)),
                upper_bound: Some(R32!(1, 4)),
                shift: R32!(0),
                flipped: false,
            }; 2],
            vec!["x".to_string(); 2],
            R32!(0),
        );
        assert_eq!(
            initial.compute_presolve_changes(),
            Ok(Changes {
                b: HashMap::default(),
                constraints: vec![],
                fixed_cost: R32!(1, 2),
                bounds: HashMap::default(),
                removed_variables: vec![
                    (0, RemovedVariable::Solved(R32!(1, 4))),
                    (1, RemovedVariable::Solved(R32!(1, 4))),
                ],
                constraints_marked_removed: vec![0, 1, 2],
            }),
        );
    }
}
