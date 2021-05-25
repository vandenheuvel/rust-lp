//! # Reducing column ranges of a sparse integer matrix
//!
//! We would like to reduce the sum of all column ranges of a sparse integer matrix. More
//! specifically, this concerns a sparse matrix in the sense that each `(row, column)` combination
//! may or may not contain a value, but if it contains a value, it can be zero.
//!
//! This problem formulation introduces a non-negative variable for each row, which will be
//! subtracted from each value in the corresponding row. Then, a variable for the minimum of each
//! column, as well as a variable for the maximum of each column.
//!
//! Two constraints are introduced for each `(row, column)` combination that contains a value (that
//! value might be zero):
//!
//! * `column minimum + row subtraction amount <= matrix value`
//! * `column maximum + row subtraction amount >= matrix value`
//!
//! None of the "subtraction amount" variables are in the initial basis, but all of the column
//! minimums and column maximums are, as well as the slack variable for each of the constraints that
//! don't already contain a pivot from the column minimum or maximum.
use num_traits::Zero;
use relp_num::RationalBig;
use relp_num::RB;

use relp::algorithm::OptimizationResult;
use relp::algorithm::two_phase::matrix_provider::matrix_data::MatrixData;
use relp::algorithm::two_phase::matrix_provider::MatrixProvider;
use relp::algorithm::two_phase::phase_two;
use relp::algorithm::two_phase::strategy::pivot_rule::FirstProfitable;
use relp::algorithm::two_phase::tableau::inverse_maintenance::carry::basis_inverse_rows::BasisInverseRows;
use relp::algorithm::two_phase::tableau::inverse_maintenance::carry::Carry;
use relp::algorithm::two_phase::tableau::inverse_maintenance::InverseMaintener;
use relp::algorithm::two_phase::tableau::kind::non_artificial::NonArtificial;
use relp::algorithm::two_phase::tableau::Tableau;
use relp::data::linear_algebra::matrix::{ColumnMajor, Order};
use relp::data::linear_algebra::vector::{DenseVector, Vector};
use relp::data::linear_program::elements::VariableType;
use relp::data::linear_program::general_form::Variable;

fn main() {
    // Needs at least one `Some` in each row.
    // Values should be large enough, such that the minimal value in the matrix will always be
    // non-negative.
    let input_matrix = [
        [Some(3), Some(3), Some(3)],
        [None,    Some(3), Some(3)],
        [Some(1), Some(2), Some(3)],
    ];

    let m = input_matrix.len();
    let n = input_matrix[0].len();

    let column_extreme = |f: fn(_) -> Option<i32>| (0..n)
        .map(|j| f((0..m).filter_map(move |i| input_matrix[i][j])).unwrap())
        .collect::<Vec<_>>();
    let column_min = column_extreme(Iterator::min);
    let column_max = column_extreme(Iterator::max);

    // Variables
    let subtraction_amount = (0..m).map(|_| Variable {
            variable_type: VariableType::Continuous,
            cost: RB!(0),
            lower_bound: Some(RB!(0)),
            upper_bound: None,
            shift: RB!(0),
            flipped: false,
        })
        .collect::<Vec<_>>();
    let column_minimum = (0..n).map(|_| Variable {
            variable_type: VariableType::Continuous,
            // We subtract this variable from a corresponding variable to compute a column range.
            cost: RB!(-1),
            lower_bound: Some(RB!(0)),
            upper_bound: None,
            shift: RB!(0),
            flipped: false,
        })
        .collect::<Vec<_>>();
    let column_maximum = (0..n).map(|_| Variable {
            variable_type: VariableType::Continuous,
            cost: RB!(1),
            lower_bound: Some(RB!(0)),
            upper_bound: None,
            shift: RB!(0),
            flipped: false,
        })
        .collect::<Vec<_>>();
    // Note the order of the variables.
    let variables = [subtraction_amount, column_minimum, column_maximum].concat();

    // Constraint matrix
    // We collect the constraints row major.
    let mut row_major_constraints = Vec::new();
    let mut b = Vec::new();
    // We'll be choosing a basis column for each row as we construct them.
    let mut basis_columns = Vec::new();

    // Lowest value in each column
    let mut nr_upper_bounded_constraints = 0;
    // Only the constraint containing the variable for the lowest value of a column has an initial
    // basis column in that same row. For other constraints containing that variable, we use a slack
    // as a basis column.
    let mut had_a_min = vec![false; n];
    // We walk the entire matrix and check whether a `column_maximum_j + x_i >= f(i, j)` constraint
    // should be added.
    for j in 0..n {
        for i in 0..m {
            if let Some(value) = input_matrix[i][j] {
                // There is a value there, so we add a constraint.
                row_major_constraints.push(vec![
                    (i, RB!(1)),     // The "subtraction amount" variable index
                    (m + j, RB!(1)), // The "column minimum" variable index (there are `m` of the
                                     // "subtraction amount" variables)
                ]);
                b.push(RB!(value));
                basis_columns.push(if value == column_min[j] {
                    // Values that are one of the extreme values in a column "quality" to have their
                    // constraint contain the "column minimum" variable pivot.
                    if had_a_min[j] {
                        // We already have a basis column for the "column minimum" variable, choose
                        // the constraint slack instead.
                        m + 2 * n + nr_upper_bounded_constraints
                    } else {
                        had_a_min[j] = true;
                        // We don't yet have a basis column for the "column minimum" variable, so we
                        // choose its column as the basis column for this constraint.
                        m + j
                    }
                } else {
                    m + 2 * n + nr_upper_bounded_constraints
                });
                nr_upper_bounded_constraints += 1;
            }
        }
    }

    // Largest value in each column
    let mut nr_lower_bounded_constraints = 0;
    let mut had_a_max = vec![false; n];
    for j in 0..n {
        for i in 0..m {
            if let Some(value) = input_matrix[i][j] {
                row_major_constraints.push(vec![(i, RB!(1)), (m + n + j, RB!(1))]);
                b.push(RB!(value));
                basis_columns.push(if value == column_max[j] {
                    if had_a_max[j] {
                        m + 2 * n + nr_upper_bounded_constraints + nr_lower_bounded_constraints
                    } else {
                        had_a_max[j] = true;
                        m + n + j
                    }
                } else {
                    m + 2 * n + nr_upper_bounded_constraints + nr_lower_bounded_constraints
                });
                nr_lower_bounded_constraints += 1;
            }
        }
    }

    // Transposing constraints to column major
    let mut constraints = vec![vec![]; m + 2 * n];
    for (row_index, row) in row_major_constraints.into_iter().enumerate() {
        for (column_index, value) in row {
            constraints[column_index].push((row_index, value));
        }
    }
    let constraints = ColumnMajor::new(constraints, nr_upper_bounded_constraints + nr_lower_bounded_constraints, m + 2 * n);
    let b = DenseVector::new(b, nr_upper_bounded_constraints + nr_lower_bounded_constraints);

    // Create the datastructure that will serve as a `MatrixProvider`
    let matrix = MatrixData::new(
        &constraints,
        &b,
        Vec::with_capacity(0),
        0, 0, nr_upper_bounded_constraints, nr_lower_bounded_constraints,
        &variables,
    );

    // We will maintain a basis inverse explicitly
    type IM = Carry<RationalBig, BasisInverseRows<RationalBig>>;
    let inverse_maintainer = IM::from_basis(&basis_columns, &matrix);

    // We create a tableau using the constructed matrix. The basis is initialized using the
    // specified columns.
    let mut tableau = Tableau::<_, NonArtificial<_>>::new_with_inverse_maintainer(
        &matrix, inverse_maintainer, basis_columns.into_iter().collect(),
    );

    // We apply primal simplex to improve the solution.
    let result = phase_two::primal::<_, _, FirstProfitable>(&mut tableau);
    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let solution = matrix.reconstruct_solution(vector);
            let shifts = (0..m)
                .map(|i| match solution.get(i) {
                    None => Zero::zero(),
                    Some(value) => value.clone(),
                })
                .collect::<Vec<_>>();
            println!("{:?}", shifts);
        },
        _ => panic!("We started with a feasible solution, and has at least value 0."),
    }
}
