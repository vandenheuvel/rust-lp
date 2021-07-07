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
use relp::data::linear_program::general_form::Scalable;
use relp::io::import;

/// An exact linear program solver written in rust.
#[derive(Clap)]
#[clap(version = "0.0.4", author = "Bram van den Heuvel <bram@vandenheuvel.online>")]
struct Opts {
    /// File containing the problem description
    problem_file: String,
    /// Disable presolving
    #[clap(long)]
    no_presolve: bool,
    /// Disable scaling
    #[clap(long)]
    no_scale: bool,
}

fn main() {
    let opts: Opts = Opts::parse();

    let path = Path::new(&opts.problem_file);
    println!("Reading problem file: \"{}\"...", path.to_string_lossy());

    let mps = import(path)
        .expect("Couldn't parse the file");

    let mut general: GeneralForm<RationalBig> = mps.try_into()
        .expect("Problem is inconsistent");

    if !opts.no_presolve {
        println!("Presolving...");
        if let Err(program_type) = general.presolve() {
            match program_type {
                LinearProgramType::FiniteOptimum(solution) => {
                    println!("Solution computed during presolve.\n{}", solution.to_string())
                },
                LinearProgramType::Infeasible => println!("Problem is not feasible."),
                LinearProgramType::Unbounded => println!("Problem is unbounded."),
            }
            exit(0);
        }
    }

    let constraint_type_counts = general.standardize();

    let scaling = if !opts.no_scale {
        println!("Scaling...");
        Some(general.scale())
    } else {
        None
    };
    println!("{:?}", scaling);

    let data = general.derive_matrix_data(constraint_type_counts);

    println!("Solving relaxation...");
    let result = data.solve_relaxation::<Carry<RationalBig, LUDecomposition<_>>>();

    println!("Solution computed:");
    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let mut reconstructed = data.reconstruct_solution(vector);
            if let Some(scaling) = scaling {
                scaling.scale_back(&mut reconstructed);
            }
            let solution = general.compute_full_solution_with_reduced_solution(reconstructed);
            println!("{}", solution.to_string());
        },
        OptimizationResult::Infeasible => println!("Problem is not feasible."),
        OptimizationResult::Unbounded => println!("Problem is unbounded."),
    }
}
