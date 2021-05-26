use std::collections::HashSet;

use relp_num::{Rational64, RationalBig};
use relp_num::Field;
use relp_num::NonZero;
use relp_num::RB;

use crate::algorithm::OptimizationResult;
use crate::algorithm::two_phase::{phase_one, phase_two};
use crate::algorithm::two_phase::matrix_provider::matrix_data::MatrixData;
use crate::algorithm::two_phase::phase_one::{Rank, RankedFeasibilityResult};
use crate::algorithm::two_phase::strategy::pivot_rule::FirstProfitable;
use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::{BasisInverse, Carry};
use crate::algorithm::two_phase::tableau::inverse_maintenance::carry::basis_inverse_rows::BasisInverseRows;
use crate::algorithm::two_phase::tableau::inverse_maintenance::InverseMaintener;
use crate::algorithm::two_phase::tableau::kind::artificial::partially::Partially;
use crate::algorithm::two_phase::tableau::kind::non_artificial::NonArtificial;
use crate::algorithm::two_phase::tableau::Tableau;
use crate::data::linear_algebra::matrix::{ColumnMajor, Order, Sparse};
use crate::data::linear_algebra::vector::{DenseVector, SparseVector};
use crate::data::linear_algebra::vector::test::TestVector;
use crate::data::linear_program::elements::VariableType;
use crate::data::linear_program::general_form::Variable;

type T = Rational64;
type S = RationalBig;

#[test]
fn conversion_pipeline() {
    let (constraints, b, variables) = create_matrix_data_data::<T>();
    let matrix_data_form = matrix_data_form(&constraints, &b, &variables);

    // Artificial tableau form
    let artificial_tableau_form_computed = Tableau::<_, Partially<_>>::new(&matrix_data_form);
    assert_eq!(artificial_tableau_form_computed, artificial_tableau_form(&matrix_data_form));

    // Get to a basic feasible solution
    let feasibility_result = phase_one::primal::<_, _, MatrixData<T>, FirstProfitable>(artificial_tableau_form_computed);
    let mut tableau_form_computed = match feasibility_result {
        RankedFeasibilityResult::Feasible {
            rank,
            nr_artificial_variables,
            inverse_maintainer,
            basis,
        } => {
            assert_eq!(rank, Rank::Full);
            Tableau::<_, NonArtificial<_>>::from_artificial(
                inverse_maintainer,
                nr_artificial_variables,
                basis,
                &matrix_data_form,
            )
        },
        _ => panic!(),
    };

    // Non-artificial tableau form
    assert_eq!(tableau_form_computed, tableau_form(&matrix_data_form));

    // Get to a basic feasible solution
    let result = phase_two::primal::<_, _, FirstProfitable>(&mut tableau_form_computed);
    assert_eq!(result, OptimizationResult::FiniteOptimum(SparseVector::from_test_tuples(vec![
        (1, 0.5f64),
        (3, 2.5f64),
        (4, 1.5f64),
    ], 5)));
}

pub fn create_matrix_data_data<'a, T: Field + From<u8> + From<&'a u8> + NonZero>(
) -> (Sparse<T, T, ColumnMajor>, DenseVector<T>, Vec<Variable<T>>) {
    let constraints = ColumnMajor::from_test_data(
        &vec![
            vec![3, 2, 1, 0, 0],
            vec![5, 1, 1, 1, 0],
            vec![2, 5, 1, 0, 1],
        ],
        5,
    );
    
    let b = DenseVector::from_test_data(vec![
        1,
        3,
        4,
    ]);

    let variables = vec![
        Variable {
            cost: T::from(1),
            lower_bound: Some(T::from(0)),
            upper_bound: None,
            shift: T::from(0),
            variable_type: VariableType::Continuous,
            flipped: false,
        }; 5
    ];

    (constraints, b, variables)
}

pub fn matrix_data_form<'a>(
    constraints: &'a Sparse<T, T, ColumnMajor>,
    b: &'a DenseVector<T>,
    variables: &'a Vec<Variable<T>>,
) -> MatrixData<'a, T> {
    MatrixData::new(
        &constraints,
        &b,
        Vec::with_capacity(0),
        3,
        0,
        0,
        0,
        &variables,
    )
}

pub fn artificial_tableau_form<'a>(
    data: &'a MatrixData<'a, T>,
) -> Tableau<Carry<S, BasisInverseRows<S>>, Partially<MatrixData<'a, T>>> {
    let m = 3;
    let artificials = (0..m).collect::<Vec<_>>();
    let carry = {
        let minus_objective = RB!(-8);
        let minus_pi = DenseVector::from_test_data(vec![-1; 3]);
        let b = DenseVector::from_test_data(vec![1, 3, 4]);
        let basis_indices = artificials.clone();
        let basis_inverse_rows = BasisInverseRows::identity(m);
        Carry::new(minus_objective, minus_pi, b, basis_indices, basis_inverse_rows)
    };
    let basis_columns = (0..m).map(|i| carry.basis_column_index_for_row(i)).collect();

    Tableau::<_, Partially<_>>::new_with_basis(
        data,
        carry,
        basis_columns,
        artificials,
    )
}

pub fn tableau_form<'a>(
    data: &'a MatrixData<'a, T>,
) -> Tableau<Carry<S, BasisInverseRows<S>>, NonArtificial<MatrixData<'a, T>>> {
    let carry = {
        let minus_objective = RationalBig::from((-9, 2));
        let minus_pi = DenseVector::from_test_data(vec![(5, 2), (-1, 1), (-1, 1)]);
        let b = DenseVector::from_test_data(vec![
            (1, 2),
            (5, 2),
            (3, 2),
        ]);
        let basis_indices = vec![1, 3, 4];
        let basis_inverse_rows = BasisInverseRows::new(vec![
            SparseVector::from_test_data(vec![(1, 2), (0, 1), (0, 1)]),
            SparseVector::from_test_data(vec![(-1, 2), (1, 1), (0, 1)]),
            SparseVector::from_test_data(vec![(-5, 2), (0, 1), (1, 1)]),
        ]);

        Carry::new(minus_objective, minus_pi, b, basis_indices, basis_inverse_rows)
    };
    let basis_columns = {
        let mut basis_columns = HashSet::new();
        basis_columns.insert(1);
        basis_columns.insert(3);
        basis_columns.insert(4);

        basis_columns
    };

    Tableau::<_, NonArtificial<_>>::new_with_inverse_maintainer(
        data,
        carry,
        basis_columns,
    )
}
