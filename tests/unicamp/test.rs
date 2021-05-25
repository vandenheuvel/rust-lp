use relp::data::linear_program::solution::Solution;
use relp_num::{Rational64, RationalBig};
use relp_num::RB;

use crate::unicamp::solve;

#[test]
fn model_data_1() {
    let result = solve::<RationalBig, Rational64>("model_data_1");
    assert!(result.is_probably_equal_to(&Solution::new(
        RB!(123, 38),  // GLPK
        vec![
            ("COL01".to_string(), RB!(5, 2)),
            ("COL02".to_string(), RB!(0)),
            ("COL03".to_string(), RB!(0)),
            ("COL04".to_string(), RB!(9, 14)),
            ("COL05".to_string(), RB!(1, 2)),
            ("COL06".to_string(), RB!(4)),
            ("COL07".to_string(), RB!(0)),
            ("COL08".to_string(), RB!(5, 19)),
        ],
    ), 0.5));
}

#[test]
#[ignore = "In this implementation, at least one RHS is needed."]
fn model_data_2() {
    let result = solve::<RationalBig, Rational64>("model_data_2");
    assert_eq!(result, Solution::new(  // GLPK
        RB!(0),
        vec![
            ("DCOL1".to_string(), RB!(0)),
        ],
    ));
}

#[test]
fn model_data_3_1() {
    let result = solve::<RationalBig, Rational64>("model_data_3_1");
    assert_eq!(result, Solution::new(  // GLPK
        RB!(70),
        vec![
            ("SUP1".to_string(), RB!(200, 3)),
            ("SUP2".to_string(), RB!(100, 3)),
            ("SUP3".to_string(), RB!(100)),
        ],
    ));
}

#[test]
fn model_data_3_2() {
    let result = solve::<RationalBig, Rational64>("model_data_3_2");
    assert_eq!(result, Solution::new(  // GLPK
        RB!(180),
        vec![
            ("SUP1".to_string(), RB!(25)),
            ("SUP2".to_string(), RB!(75)),
        ],
    ));
}

#[test]
fn model_data_3_3() {
    let result = solve::<RationalBig, Rational64>("model_data_3_3");
    assert_eq!(result, Solution::new(  // GLPK
        RB!(245),
        vec![
            ("SUP1".to_string(), RB!(100)),
            ("SUP2".to_string(), RB!(150)),
        ],
    ));
}

#[test]
fn model_data_3_4() {
    let result = solve::<RationalBig, Rational64>("model_data_3_4");
    assert_eq!(result, Solution::new(  // GLPK
        RB!(2250),
        vec![
            ("RAW1".to_string(), RB!(5)),
            ("RAW2".to_string(), RB!(3)),
            ("RAW3".to_string(), RB!(4)),
            ("PRODUCT".to_string(), RB!(500)),
        ],
    ));
}

#[test]
fn model_data_4() {
    let result = solve::<RationalBig, Rational64>("model_data_4");
    assert_eq!(result, Solution::new(  // GLPK
        RB!(7, 1),
        vec![
            ("COL01".to_string(), RB!(1)),
            ("COL02".to_string(), RB!(2)),
            ("COL03".to_string(), RB!(2)),
        ],
    ));
}

#[test]
#[ignore = "This problem type is not supported."]
fn model_data_5() {
    let result = solve::<RationalBig, Rational64>("model_data_5");
    assert_eq!(result, Solution::new(  // GLPK
        RB!(332, 1),
        vec![
            ("COL01".to_string(), RB!(1, 1)),
            ("COL02".to_string(), RB!(2, 1)),
            ("COL03".to_string(), RB!(2, 1)),
        ],
    ));
}

#[test]
fn model_data_6() {
    let result = solve::<RationalBig, Rational64>("model_data_6");
    assert!(result.is_probably_equal_to(&Solution::new(  // GLPK
        RB!(28, 1),
        vec![
            ("X00".to_string(), RB!(0)),
            ("X01".to_string(), RB!(1)),
            ("X02".to_string(), RB!(1)),
            ("X03".to_string(), RB!(0)),
            ("X04".to_string(), RB!(0)),
            ("X05".to_string(), RB!(0)),
            ("X06".to_string(), RB!(0)),
            ("X07".to_string(), RB!(0)),
            ("X10".to_string(), RB!(1)),
            ("X11".to_string(), RB!(0)),
            ("X12".to_string(), RB!(0)),
            ("X13".to_string(), RB!(2)),
            ("X14".to_string(), RB!(0)),
            ("X15".to_string(), RB!(0)),
            ("X16".to_string(), RB!(0)),
            ("X17".to_string(), RB!(0)),
            ("X20".to_string(), RB!(1)),
            ("X21".to_string(), RB!(0)),
            ("X22".to_string(), RB!(0)),
            ("X23".to_string(), RB!(3)),
            ("X24".to_string(), RB!(0)),
            ("X25".to_string(), RB!(0)),
            ("X26".to_string(), RB!(0)),
            ("X27".to_string(), RB!(0)),
        ],
    ), 0.5));
}

#[test]
#[ignore = "Identical to model_data_1."]
fn model_data_7() {
    let _result = solve::<RationalBig, Rational64>("model_data_7");
}

#[test]
#[ignore = "Unsupported modification of to model_data_7."]
fn model_data_8() {
    let _result = solve::<RationalBig, Rational64>("model_data_8");
}

#[test]
#[ignore = "Unnamed problem files are not supported."]
fn model_data_9() {
    let result = solve::<RationalBig, Rational64>("model_data_9");
    let expected = Solution::new(  // GLPK
        RB!(-100, 1),
        vec![
            ("C0000001".to_string(), RB!(0)),
            ("C0000002".to_string(), RB!(1)),
            ("C0000003".to_string(), RB!(1)),
            ("C0000004".to_string(), RB!(0)),
        ],
    );
    assert!(result.is_probably_equal_to(&expected, 0.5));
}
