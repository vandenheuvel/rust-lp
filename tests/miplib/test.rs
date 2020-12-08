use super::test;

#[test]
fn test_50v() {
    test("50v-10", 2879.065687f64, 1e-3);
}  // GLPK

#[test]
#[ignore = "Incorrectly determined infeasible"]
fn test_30n() {
    test("30n20b8", 43.33557298f64, 1e-3);
}  // GLPK

#[test]
#[ignore = "Too computationally expensive"]
fn test_acc() {
    test("acc-tight4", 0f64, 1e-3);
}  // GLPK
