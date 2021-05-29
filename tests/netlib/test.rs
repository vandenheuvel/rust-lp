use relp_num::Abs;
use relp_num::RB;

use crate::netlib::solve;

#[test]
fn test_ADLITTLE() {
    let result = solve("ADLITTLE");
    let expected = 2.254949632e+05; // Gurobi
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-3)); // Gurobi
}

#[test]
fn test_AFIRO() {
    let result = solve("AFIRO");
    let expected = -464.75314; // Coin LP 1.17.6
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-3));
}

#[test]
fn test_SC50A() {
    let result = solve("SC50A");
    let expected = -6.457507706e+01; // GLPK 4.65
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-5));
}

#[test]
fn test_SC50B() {
    let result = solve("SC50B");
    let expected = -70; // GLPK 4.65
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-10));
}

#[test]
fn test_KB2() {
    let result = solve("KB2");
    let expected = -1.749900130e+03; // GLPK 4.65
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-3));
}

#[test]
fn test_SC105() {
    let result = solve("SC105");
    let expected = -5.220206121e+01; // GLPK 4.65
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-3));
}

#[test]
fn test_STOCFOR1() {
    let result = solve("STOCFOR1");
    let expected = -4.113197622e+04; // GLPK 4.65
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-3));
}

#[test]
fn test_BLEND() {
    let result = solve("BLEND");
    let expected = -30.81215; // GLPK 4.65
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-3));
}

#[test]
fn test_SCAGR7() {
    let result = solve("SCAGR7");
    let expected = -2.331389824e+06; // GLPK 4.65
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-1));
}

#[test]
fn test_SC205() {
    let result = solve("SC205");
    let expected = -5.220206121e+01; // GLPK 4.65
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-5));
}

#[test]
fn test_SHARE2B() {
    let result = solve("SHARE2B");
    let expected = -4.157322407e+02; // GLPK 4.65
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-5));
}

#[test]
fn test_RECIPELP() {
    let result = solve("RECIPELP");
    let expected = -0.266616e3; // Koch - The final Netlib-LP results
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-2));
}

#[test]
fn test_LOTFI() {
    let result = solve("LOTFI");
    let expected = -0.2526470606188e2; // Koch - The final Netlib-LP results
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-6));
}

#[test]
fn test_VTP_BASE() {
    let result = solve("VTP-BASE");
    let expected = 0.1298314624613613657395984384889e6; // Koch - The final Netlib-LP results
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-2));
}

#[test]
fn test_SHARE1B() {
    let result = solve("SHARE1B");
    let expected = -0.7658931857918568112797274346007e5; // Koch - The final Netlib-LP results
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-3));
}

#[test]
fn test_BOEING2() {
    let result = solve("BOEING2");
    let expected = -0.31501872801520287870462195913263e3; // Koch - The final Netlib-LP results
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-3));
}

#[test]
fn test_BORE3D() {
    let result = solve("BORE3D");
    let expected = 0.13730803942084927215581987251301e4; // Koch - The final Netlib-LP results
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-2));
}

#[test]
#[ignore = "Incorrect optimal value."]
fn test_SCORPION() {
    let result = solve("SCORPION");
    let expected = 0.18781248227381066296479411763586e4; // Koch - The final Netlib-LP results
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-2));
}

#[test]
#[ignore = "Too computationally intensive."]
fn test_GREENBEA() {
    let result = solve("GREENBEA");
    let expected = -0.72555248129845987457557870574845e8; // Koch - The final Netlib-LP results
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e0));
}

#[test]
#[ignore = "Too computationally intensive."]
fn test_GREENBEB() {
    let result = solve("GREENBEB");
    let expected = -0.43022602612065867539213672544432e7; // Koch - The final Netlib-LP results
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e1));
}

#[test]
#[ignore = "Too computationally intensive."]
fn test_25FV47() {
    let result = solve("25FV47");
    let expected = 5.5018459e+03; // Gurobi
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-5));
}

#[test]
#[ignore = "Too computationally intensive, perhaps underflows in `compute_full_solution_with_reduced_solution`."]
fn test_80BAU3B() {
    let result = solve("80BAU3B");
    let expected = 9.872241924e+05; // Gurobi
    assert!((result.objective_value - RB!(expected)).abs() < RB!(1e-5)); // Gurobi
}
