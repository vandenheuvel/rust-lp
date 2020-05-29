//! Simple linear program.
//!
//! From https://en.wikipedia.org/w/index.php?title=MPS_(format)&oldid=941892011
use std::collections::HashSet;
use std::convert::{TryFrom, TryInto};

use num::rational::Ratio;
use num::FromPrimitive;

use crate::R32;
use crate::algorithm::simplex::data::{Artificial, CarryMatrix, NonArtificial, Tableau};
use crate::algorithm::simplex::logic::{artificial_primal, FeasibilityResult, OptimizationResult, primal, Rank};
use crate::algorithm::simplex::matrix_provider::matrix_data::{MatrixData, Variable as MatrixDataVariable};
use crate::algorithm::simplex::matrix_provider::MatrixProvider;
use crate::algorithm::simplex::strategy::pivot_rule::FirstProfitable;
use crate::data::linear_algebra::matrix::{MatrixOrder, RowMajorOrdering};
use crate::data::linear_algebra::vector::{DenseVector, SparseVector};
use crate::data::linear_algebra::vector::test::TestVector;
use crate::data::linear_program::elements::{ConstraintType, Objective, VariableType};
use crate::data::linear_program::general_form::GeneralForm;
use crate::data::linear_program::general_form::Variable as GeneralFormVariable;
use crate::io::mps::{Bound, BoundType, Constraint, Rhs, Variable};
use crate::io::mps::parsing::{into_atom_lines, UnstructuredBound, UnstructuredColumn, UnstructuredMPS, UnstructuredRhs, UnstructuredRow};
use crate::io::mps::structuring::MPS;

type T = Ratio<i32>;

#[test]
fn test_conversion_pipeline() {
    let input = &MPS_LITERAL_STRING;

    // Unstructured MPS
    let lines = into_atom_lines(input);
    let result = UnstructuredMPS::try_from(lines);
    assert!(result.is_ok());
    let unstructured_mps_computed: UnstructuredMPS<T> = result.unwrap();
    assert_eq!(unstructured_mps_computed, unstructured_mps_form());

    // MPS
    let result = unstructured_mps_computed.try_into();
    assert!(result.is_ok());
    let structured_mps_computed: MPS<T> = result.unwrap();
    assert_eq!(structured_mps_computed, mps_form());

    // General form
    let result = structured_mps_computed.try_into();
    assert!(result.is_ok());
    let mut general_form_computed: GeneralForm<_> = result.unwrap();
    assert_eq!(general_form_computed, general_form());

    // Matrix data form
    let result = general_form_computed.derive_matrix_data();
    // General form, canonicalized
    assert_eq!(general_form_computed, general_form_canonicalized());
    // The resulting matrix data form
    assert!(result.is_ok());
    let matrix_data_form_computed = result.unwrap();
    assert_eq!(matrix_data_form_computed, matrix_data_form());

    // Artificial tableau form
    let mut artificial_tableau_form_computed = Tableau::<_, Artificial, _>::new(&matrix_data_form_computed);
    assert_eq!(artificial_tableau_form_computed, artificial_tableau_form(&matrix_data_form()));

    // Get to a basic feasible solution
    let feasibility_result = artificial_primal::<_, _, FirstProfitable>(&mut artificial_tableau_form_computed);
    assert_eq!(feasibility_result, FeasibilityResult::Feasible(Rank::Full));

    // Non-artificial tableau form
    let mut tableau_form_computed = Tableau::<_, NonArtificial, _>::from_artificial(artificial_tableau_form_computed);
    assert_eq!(tableau_form_computed, tableau_form(&matrix_data_form()));

    // Get to a basic feasible solution
    let result = primal::<Ratio<i32>, _, FirstProfitable>(&mut tableau_form_computed);
    assert_eq!(result, OptimizationResult::FiniteOptimum(SparseVector::from_test_tuples(vec![
        (0, 2_f64),
        (5, 2_f64),
        (6, 2_f64),
    ], 7)));
}

const MPS_LITERAL_STRING: &str = "NAME          TESTPROB
ROWS
 N  COST
 L  LIM1
 G  LIM2
 E  MYEQN
COLUMNS
    XONE      COST                 1   LIM1                 1
    XONE      LIM2                 1
    MARKER0   'MARKER'      'INTORG'
    YTWO      COST                 4   LIM1                 1
    YTWO      MYEQN               -1
    MARKER0   'MARKER'      'INTEND'
    ZTHREE    COST                 9   LIM2                 1
    ZTHREE    MYEQN                1
RHS
    RHS1      LIM1                 5   LIM2                10
    RHS1      MYEQN                7
BOUNDS
 UP BND1      XONE                 4
 LO BND1      YTWO                -1
 UP BND1      YTWO                 1
ENDATA";


/// Build the expected `MPS` instance, corresponding to the MPS file string.
pub(super) fn unstructured_mps_form() -> UnstructuredMPS<'static, T> {
    UnstructuredMPS {
        name: "TESTPROB",
        cost_row_name: "COST",
        rows: vec![
            UnstructuredRow {
                name: "LIM1",
                constraint_type: ConstraintType::Less,
            },
            UnstructuredRow {
                name: "LIM2",
                constraint_type: ConstraintType::Greater,
            },
            UnstructuredRow {
                name: "MYEQN",
                constraint_type: ConstraintType::Equal,
            },
        ],
        columns: vec![
            UnstructuredColumn {
                name: "XONE",
                variable_type: VariableType::Continuous,
                row_name: "COST",
                value: R32!(1),
            },
            UnstructuredColumn {
                name: "XONE",
                variable_type: VariableType::Continuous,
                row_name: "LIM1",
                value: R32!(1),
            },
            UnstructuredColumn {
                name: "XONE",
                variable_type: VariableType::Continuous,
                row_name: "LIM2",
                value: R32!(1),
            },
            UnstructuredColumn {
                name: "YTWO",
                variable_type: VariableType::Integer,
                row_name: "COST",
                value: R32!(4),
            },
            UnstructuredColumn {
                name: "YTWO",
                variable_type: VariableType::Integer,
                row_name: "LIM1",
                value: R32!(1),
            },
            UnstructuredColumn {
                name: "YTWO",
                variable_type: VariableType::Integer,
                row_name: "MYEQN",
                value: -R32!(1),
            },
            UnstructuredColumn {
                name: "ZTHREE",
                variable_type: VariableType::Continuous,
                row_name: "COST",
                value: R32!(9),
            },
            UnstructuredColumn {
                name: "ZTHREE",
                variable_type: VariableType::Continuous,
                row_name: "LIM2",
                value: R32!(1),
            },
            UnstructuredColumn {
                name: "ZTHREE",
                variable_type: VariableType::Continuous,
                row_name: "MYEQN",
                value: R32!(1),
            },
        ],
        rhss: vec![
            UnstructuredRhs {
                name: "RHS1",
                row_name: "LIM1",
                value: R32!(5),
            },
            UnstructuredRhs {
                name: "RHS1",
                row_name: "LIM2",
                value: R32!(10),
            },
            UnstructuredRhs {
                name: "RHS1",
                row_name: "MYEQN",
                value: R32!(7),
            },
        ],
        ranges: vec![],
        bounds: vec![
            UnstructuredBound {
                name: "BND1",
                bound_type: BoundType::UpperContinuous(R32!(4)),
                column_name: "XONE",
            },
            UnstructuredBound {
                name: "BND1",
                bound_type: BoundType::LowerContinuous(-R32!(1)),
                column_name: "YTWO",
            },
            UnstructuredBound {
                name: "BND1",
                bound_type: BoundType::UpperContinuous(R32!(1)),
                column_name: "YTWO",
            },
        ],
    }
}

/// Build the expected `MPS` instance, corresponding to the MPS file string.
pub fn mps_form() -> MPS<T> {
    let name = "TESTPROB".to_string();
    let cost_row_name = "COST".to_string();
    let cost_values = vec![(0, R32!(1)), (1, R32!(4)), (2, R32!(9))];
    let row_names = vec!["LIM1", "LIM2", "MYEQN"].into_iter().map(String::from).collect();
    let rows = vec![
        Constraint { name: 0, constraint_type: ConstraintType::Less },
        Constraint { name: 1, constraint_type: ConstraintType::Greater },
        Constraint { name: 2, constraint_type: ConstraintType::Equal },
    ];
    let column_names = vec!["XONE", "YTWO", "ZTHREE"].into_iter().map(String::from).collect();
    let columns = vec![
        Variable {
            name: 0,
            variable_type: VariableType::Continuous,
            values: vec![
                (0, R32!(1)),
                (1, R32!(1)),
            ],
        },
        Variable {
            name: 1,
            variable_type: VariableType::Integer,
            values: vec![
                (0, R32!(1)),
                (2, -R32!(1)),
            ],
        },
        Variable {
            name: 2,
            variable_type: VariableType::Continuous,
            values: vec![
                (1, R32!(1)),
                (2, R32!(1)),
            ],
        },
    ];
    let ranges = vec![];
    let rhss = vec![
        Rhs {
            name: "RHS1".to_string(),
            values: vec![(0, R32!(5)), (1, R32!(10)), (2, R32!(7))],
        }
    ];
    let bounds = vec![
        Bound {
            name: "BND1".to_string(),
            values: vec![
                (BoundType::UpperContinuous(R32!(4)), 0),
                (BoundType::LowerContinuous(-R32!(1)), 1),
                (BoundType::UpperContinuous(R32!(1)), 1),
            ],
        }
    ];

    MPS::new(
        name,
        cost_row_name,
        cost_values,
        row_names,
        rows,
        column_names,
        columns,
        rhss,
        ranges,
        bounds,
    )
}

/// The linear program in expected `GeneralForm`.
pub fn general_form() -> GeneralForm<T> {
    let data = vec![
        vec![1_f64, 1_f64, 0_f64],
        vec![1_f64, 0_f64, 1_f64],
        vec![0_f64, -1_f64, 1_f64],
    ];
    let rows = RowMajorOrdering::from_test_data(&data, 3);

    let constraints = vec![
        ConstraintType::Less,
        ConstraintType::Greater,
        ConstraintType::Equal,
    ];

    let b = DenseVector::from_test_data(vec![
        5_f64,
        10_f64,
        7_f64,
    ]);

    let variables = vec![
        GeneralFormVariable {
            variable_type: VariableType::Continuous,
            cost: R32!(1),
            lower_bound: Some(R32!(0)),
            upper_bound: Some(R32!(4)),
            shift: R32!(0),
            flipped: false
        },
        GeneralFormVariable {
            variable_type: VariableType::Integer,
            cost: R32!(4),
            lower_bound: Some(R32!(-1)),
            upper_bound: Some(R32!(1)),
            shift: R32!(0),
            flipped: false
        },
        GeneralFormVariable {
            variable_type: VariableType::Continuous,
            cost: R32!(9),
            lower_bound: Some(R32!(0)),
            upper_bound: None,
            shift: R32!(0),
            flipped: false
        },
    ];
    let variable_names = vec!["XONE".to_string(), "YTWO".to_string(), "ZTHREE".to_string()];

    GeneralForm::new(
        Objective::Minimize,
        rows,
        constraints,
        b,
        variables,
        variable_names,
        R32!(0),
    )
}

pub fn general_form_canonicalized() -> GeneralForm<T> {
    let data = vec![
        vec![1_f64, 0_f64, 1_f64],
        vec![0_f64, -1_f64, 1_f64],
    ];
    let constraints = RowMajorOrdering::from_test_data(&data, 3);

    let constraint_types = vec![
        ConstraintType::Greater,
        ConstraintType::Equal,
    ];

    let b = DenseVector::from_test_data(vec![
        2_f64,
        0_f64,
    ]);

    let variables = vec![
        GeneralFormVariable {
            variable_type: VariableType::Continuous,
            cost: R32!(1),
            lower_bound: Some(R32!(0)),
            upper_bound: Some(R32!(2)),
            shift: R32!(-2),
            flipped: false
        },
        GeneralFormVariable {
            variable_type: VariableType::Integer,
            cost: R32!(4),
            lower_bound: Some(R32!(0)),
            upper_bound: Some(R32!(2)),
            shift: R32!(1),
            flipped: false
        },
        GeneralFormVariable {
            variable_type: VariableType::Continuous,
            cost: R32!(9),
            lower_bound: Some(R32!(0)),
            upper_bound: Some(R32!(2)),
            shift: R32!(-6),
            flipped: false
        },
    ];
    let variable_names = vec!["XONE".to_string(), "YTWO".to_string(), "ZTHREE".to_string()];

    GeneralForm::new(
        Objective::Minimize,
        constraints,
        constraint_types,
        b,
        variables,
        variable_names,
        -R32!(-2 * 1) - R32!(1 * 4) - R32!(-6 * 9),
    )
}

pub fn matrix_data_form() -> MatrixData<T> {
    let equal = vec![vec![(1, R32!(-1)), (2, R32!(1))]];
    let upper_bounded = vec![];
    let lower_bounded = vec![vec![(0, R32!(1)), (2, R32!(1))]];
    let b = DenseVector::from_test_data(vec![
        0_f64,
        2_f64,
    ]);
    let variables = vec![
        MatrixDataVariable {
            cost: R32!(1),
            upper_bound: Some(R32!(2)),
            variable_type: VariableType::Continuous,
        },
        MatrixDataVariable {
            cost: R32!(4),
            upper_bound: Some(R32!(2)),
            variable_type: VariableType::Integer,
        },
        MatrixDataVariable {
            cost: R32!(9),
            upper_bound: Some(R32!(2)),
            variable_type: VariableType::Continuous,
        },
    ];
    MatrixData::new(
        equal,
        upper_bounded,
        lower_bounded,
        b,
        variables,
        Vec::with_capacity(0),
    )
}

pub fn artificial_tableau_form<MP: MatrixProvider<T>>(provider: &MP) -> Tableau<T, Artificial, MP> {
    let m = 5;
    let carry = {
        let minus_objective = R32!(-8);
        let minus_pi = DenseVector::from_test_data(vec![-1_f64; m]);
        let b = DenseVector::from_test_data(vec![0_f64, 2_f64, 2_f64, 2_f64, 2_f64]);
        let basis_inverse_rows = (0..m)
            .map(|i| SparseVector::standard_basis_vector(i, m))
            .collect();
        CarryMatrix::create(minus_objective, minus_pi, b, basis_inverse_rows)
    };
    let basis_indices = (0..m).collect();
    let basis_columns = (0..m).collect();

    Tableau::new_with_basis(
        &provider,
        carry,
        basis_indices,
        basis_columns,
    )
}

pub fn tableau_form<MP: MatrixProvider<T>>(provider: &MP) -> Tableau<T, NonArtificial, MP> {
    let carry = {
        let minus_objective = R32!(-28);
        let minus_pi = DenseVector::from_test_data(vec![-9_f64, 0_f64, -1_f64, -13_f64, 0_f64]);
        let b = DenseVector::from_test_data(vec![2_f64, 2_f64, 2_f64, 2_f64, 0_f64]);
        let basis_inverse_rows = vec![
            SparseVector::from_test_data(vec![1_f64, 0_f64, 0_f64, 1_f64, 0_f64]),
            SparseVector::from_test_data(vec![0_f64, 0_f64, 1_f64, 0_f64, 0_f64]),
            SparseVector::from_test_data(vec![1_f64, -1_f64, 1_f64, 1_f64, 0_f64]),
            SparseVector::from_test_data(vec![0_f64, 0_f64, 0_f64, 1_f64, 0_f64]),
            SparseVector::from_test_data(vec![-1_f64, 0_f64, 0_f64, -1_f64, 1_f64]),
        ];
        CarryMatrix::<T>::create(minus_objective, minus_pi, b, basis_inverse_rows)
    };
    let basis_indices = vec![2, 0, 3, 1, 6];
    let mut basis_columns = HashSet::new();
    basis_columns.insert(2);
    basis_columns.insert(0);
    basis_columns.insert(3);
    basis_columns.insert(1);
    basis_columns.insert(6);

    Tableau::new_with_basis(
        provider,
        carry,
        basis_indices,
        basis_columns,
    )
}
