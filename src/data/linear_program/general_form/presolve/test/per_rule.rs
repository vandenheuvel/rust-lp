use relp_num::R32;
use relp_num::Rational32;

use crate::data::linear_algebra::matrix::ColumnMajor;
use crate::data::linear_algebra::matrix::Order;
use crate::data::linear_algebra::vector::{DenseVector, Vector};
use crate::data::linear_algebra::vector::test::TestVector;
use crate::data::linear_program::elements::{BoundDirection, LinearProgramType, Objective, VariableType};
use crate::data::linear_program::general_form::{GeneralForm, RangedConstraintRelation, RemovedVariable, Variable};
use crate::data::linear_program::general_form::presolve::Index;

type T = Rational32;

#[test]
fn presolve_fixed_variable_feasible() {
    let initial = GeneralForm::<_>::new(
        Objective::Minimize,
        ColumnMajor::from_test_data(&[vec![1], vec![2]], 1),
        vec![RangedConstraintRelation::Equal, RangedConstraintRelation::Greater],
        DenseVector::from_test_data(vec![1; 2]),
        vec![Variable {
            variable_type: VariableType::Continuous,
            cost: R32!(1),
            lower_bound: Some(R32!(1)),
            upper_bound: Some(R32!(1)),
            shift: R32!(0),
            flipped: false
        }],
        vec!["X".to_string()],
        R32!(7),
    );
    let mut index = Index::new(&initial).unwrap();

    index.presolve_fixed_variable(0).ok().unwrap();
    assert_eq!(index.counters.constraint, vec![0, 0]);
    assert_eq!(index.counters.variable, vec![0]);
    assert_eq!(index.updates.constraints_marked_removed, vec![0, 1]);
    assert_eq!(index.updates.removed_variables, vec![(0, RemovedVariable::Solved(R32!(1)))]);
    assert_eq!(index.updates.b.get(&0), Some(&R32!(0)));
    assert_eq!(index.updates.b.get(&1), Some(&R32!(-1)));
    assert_eq!(index.updates.fixed_cost, R32!(1) * R32!(1));
    assert!(!index.queues.are_empty());
}

#[test]
fn presolve_fixed_variable_infeasible() {
    let initial = GeneralForm::<_>::new(
        Objective::Minimize,
        ColumnMajor::from_test_data(&[vec![1], vec![2]], 1),
        vec![RangedConstraintRelation::Equal; 2],
        DenseVector::from_test_data(vec![1; 2]),
        vec![Variable {
            variable_type: VariableType::Continuous,
            cost: R32!(1),
            lower_bound: Some(R32!(1)),
            upper_bound: Some(R32!(1)),
            shift: R32!(0),
            flipped: false
        }],
        vec!["X".to_string()],
        R32!(7),
    );
    let mut index = Index::new(&initial).unwrap();

    assert_eq!(index.presolve_fixed_variable(0), Err(LinearProgramType::Infeasible));
}

#[test]
fn presolve_simple_bound_constraint() {
    let initial = GeneralForm::<T>::new(
        Objective::Minimize,
        ColumnMajor::from_test_data(&vec![vec![1]; 2], 1),
        vec![RangedConstraintRelation::Equal; 2],
        DenseVector::from_test_data(vec![2; 2]),
        vec![Variable {
            variable_type: VariableType::Continuous,
            cost: R32!(1),
            lower_bound: None,
            upper_bound: None,
            shift: R32!(0),
            flipped: false
        }],
        vec!["X".to_string()],
        R32!(0),
    );
    let mut index = Index::new(&initial).unwrap();

    assert_eq!(index.presolve_bound_constraint(0), Ok(()));
    assert_eq!(index.counters.constraint, vec![0, 1]);
    assert_eq!(index.counters.variable, vec![1]);
    assert_eq!(index.updates.constraints_marked_removed, vec![0]);
    assert_eq!(index.updates.removed_variables, vec![]);
    assert_eq!(index.updates.bounds.get(&(0, BoundDirection::Lower)), Some(&R32!(2)));
    assert_eq!(index.updates.bounds.get(&(0, BoundDirection::Upper)), Some(&R32!(2)));
    assert!(!index.queues.are_empty());
    assert_eq!(initial.b, DenseVector::from_test_data(vec![2; 2]));
}

#[test]
fn presolve_constraint_if_slack_with_suitable_bounds() {
    let create = |constraint_type, lower, upper| {
        let nr_variables = 3;
        GeneralForm::<T>::new(
            Objective::Minimize,
            ColumnMajor::from_test_data(&[vec![2; nr_variables]], nr_variables),
            vec![constraint_type],
            DenseVector::from_test_data(vec![3]),
            vec![Variable {
                variable_type: VariableType::Continuous,
                cost: R32!(0),
                lower_bound: lower,
                upper_bound: upper,
                shift: R32!(0),
                flipped: false,
            }; nr_variables],
            vec!["X".to_string(); nr_variables],
            R32!(0),
        )
    };

    // Read a range variable
    let initial = create(RangedConstraintRelation::Equal, Some(R32!(1)), Some(R32!(2)));
    let mut index = Index::new(&initial).unwrap();
    assert!(index.presolve_slack(0).is_ok());
    assert_eq!(index.counters.constraint, vec![2]);
    assert_eq!(index.counters.variable, vec![0, 1, 1]);
    assert_eq!(index.updates.constraints_marked_removed, vec![]);
    assert_eq!(index.updates.removed_variables, vec![(0, RemovedVariable::FunctionOfOthers {
        constant: R32!(3, 2),
        coefficients: vec![(1, R32!(1)), (2, R32!(1))],
    })]);

    // Change a thing
    let initial = create(RangedConstraintRelation::Equal, Some(R32!(1)), None);
    let mut index = Index::new(&initial).unwrap();
    assert!(index.presolve_slack(0).is_ok());
    assert_eq!(index.counters.constraint, vec![2]);
    assert_eq!(index.counters.variable, vec![0, 1, 1]);
    assert_eq!(index.updates.constraints_marked_removed, vec![]);
    assert_eq!(index.updates.constraint_type(0), &RangedConstraintRelation::<Rational32>::Less);
    assert_eq!(index.updates.removed_variables, vec![(0, RemovedVariable::FunctionOfOthers {
        constant: R32!(3, 2),
        coefficients: vec![(1, R32!(1)), (2, R32!(1))],
    })]);
    assert_eq!(index.updates.b(0), &R32!(3 - 2));

    let initial = create(RangedConstraintRelation::Equal, None, Some(R32!(1)));
    let mut index = Index::new(&initial).unwrap();
    assert!(index.presolve_slack(0).is_ok());
    assert_eq!(index.counters.constraint, vec![2]);
    assert_eq!(index.counters.variable, vec![0, 1, 1]);
    assert_eq!(index.updates.constraints_marked_removed, vec![]);
    assert_eq!(index.updates.constraint_type(0), &RangedConstraintRelation::<Rational32>::Greater);
    assert_eq!(index.updates.removed_variables, vec![(0, RemovedVariable::FunctionOfOthers {
        constant: R32!(3, 2),
        coefficients: vec![(1, R32!(1)), (2, R32!(1))],
    })]);
    assert_eq!(index.updates.b(0), &R32!(3 - 2));

    let initial = create(RangedConstraintRelation::Greater, Some(R32!(1)), None);
    let mut index = Index::new(&initial).unwrap();
    assert!(index.presolve_slack(0).is_ok());
    assert_eq!(index.counters.constraint, vec![0]);
    assert_eq!(index.counters.variable, vec![0; 3]);
    assert_eq!(index.updates.constraints_marked_removed, vec![0]);
    assert_eq!(initial.constraint_types.len(), 1); // Removed after
    assert_eq!(index.updates.removed_variables, vec![
        (1, RemovedVariable::Solved(R32!(1))),
        (2, RemovedVariable::Solved(R32!(1))),
        (0, RemovedVariable::FunctionOfOthers {
            constant: R32!(3, 2),
            coefficients: vec![(1, R32!(1)), (2, R32!(1))],
        }),
    ]);
    assert_eq!(initial.b.len(), 1); // Removed after
}
