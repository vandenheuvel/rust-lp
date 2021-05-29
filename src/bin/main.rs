use std::convert::TryInto;
use std::path::Path;
use std::process::exit;

use clap::Clap;
use relp_num::RationalBig;

use relp::algorithm::{OptimizationResult, SolveRelaxation};
use relp::algorithm::two_phase::matrix_provider::MatrixProvider;
use relp::algorithm::two_phase::tableau::inverse_maintenance::carry::Carry;
use relp::algorithm::two_phase::tableau::inverse_maintenance::carry::lower_upper::LUDecomposition;
use relp::data::linear_program::elements::LinearProgramType;
use relp::data::linear_program::general_form::GeneralForm;
use relp::io::import;

/// A linear program solver written in rust.
#[derive(Clap)]
#[clap(version = "0.0.4", author = "Bram van den Heuvel <bram@vandenheuvel.online>")]
struct Opts {
    /// File containing the problem description
    problem_file: String,
    #[clap(long)]
    no_presolve: bool,
    #[clap(long)]
    no_scale: bool,
}

fn main() {
    let opts: Opts = Opts::parse();

    let path = Path::new(&opts.problem_file);
    println!("Reading problem file: \"{}\"...", path.to_string_lossy());

    let mps = import(path)
        .expect("Couldn't parse the file.");

    let mut general: GeneralForm<RationalBig> = mps.try_into()
        .expect("Problem is inconsistent");

    println!("Presolving...");
    let data = match general.derive_matrix_data(!opts.no_presolve, !opts.no_scale) {
        Ok(presolved_data) => presolved_data,
        Err(program_type) => {
            match program_type {
                LinearProgramType::FiniteOptimum(solution) => {
                    println!("Solution computed.\n{}", solution.to_string())
                },
                LinearProgramType::Infeasible => println!("Problem is not feasible."),
                LinearProgramType::Unbounded => println!("Problem is unbounded."),
            }
            exit(0);
        },
    };

    println!("Solving relaxation...");
    let result = data.solve_relaxation::<Carry<RationalBig, LUDecomposition<_>>>();

    println!("Solution computed:");
    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let reconstructed = data.reconstruct_solution(vector);
            let solution = general.compute_full_solution_with_reduced_solution(reconstructed);
            println!("{}", solution.to_string());
        },
        OptimizationResult::Infeasible => println!("Problem is not feasible."),
        OptimizationResult::Unbounded => println!("Problem is unbounded."),
    }
}
