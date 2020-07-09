#![feature(or_patterns)]

extern crate num;

use std::convert::TryInto;
use std::env;
use std::fs::File;
use std::io::Read;
use std::path::Path;

use num_traits::BigInt;
use num_traits::rational::Ratio;

use lp::algorithm::{OptimizationResult, SolveRelaxation};
use lp::algorithm::two_phase::matrix_provider::MatrixProvider;
use lp::data::linear_program::general_form::GeneralForm;
use lp::io::import;
use rust_lp::algorithm::two_phase::matrix_provider::MatrixProvider;
use rust_lp::data::linear_program::general_form::GeneralForm;
use rust_lp::io::error::Import;
use rust_lp::io::mps;

/// Import a problem from a file.
///
/// Currently only supports the MPS filetype.
///
/// The `import` function takes a file path and returns, if successful, a struct which can be
/// converted to a linear program in general form.
///
/// # Errors
///
/// When a file extension is unknown, a file cannot be found or read, there is an inconsistency in
/// the problem file, etc. an error type is returned.
pub fn import(
    file_path: &Path
) -> Result<impl TryInto<GeneralForm<Rational64, Rational64>>, Import> {
    // Open and read the file
    let mut program = String::new();
    File::open(&file_path)
        .map_err(Import::IO)?
        .read_to_string(&mut program)
        .map_err(Import::IO)?;

    // Choose the right parser
    match file_path.extension() {
        Some(extension) => match extension.to_str() {
            Some("mps" | "SIF") => mps::parse(&program),
            Some(extension_string) => Err(Import::FileExtension(format!(
                "Could not recognise file extension \"{}\" of file: {:?}",
                extension_string, file_path
            ))),
            None => Err(Import::FileExtension(format!(
                "Could not convert OsStr to &str, probably invalid unicode: {:?}",
                extension
            ))),
        },
        None => Err(Import::FileExtension(format!(
            "Could not read extension from file path: {:?}",
            file_path
        ))),
    }
}

fn main() {
    type T = Ratio<BigInt>;

    let name = env::args().nth(1).unwrap();
    let pathname = format!("./test_files/{}.mps", name);
    let path = Path::new(&pathname);

    let result = import(path).unwrap();

    let mut general_form: GeneralForm<T, T> = result.try_into().ok().unwrap();
    let matrix_data = general_form.derive_matrix_data().ok().unwrap();

    let result = matrix_data.solve_relaxation();

    match result {
        OptimizationResult::FiniteOptimum(vector) => {
            let reconstructed = matrix_data.reconstruct_solution(vector);
            let solution = general_form.compute_full_solution_with_reduced_solution(reconstructed);
            println!("{:?}", solution.objective_value);
        },
        _ => panic!(),
    }
}
