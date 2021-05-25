//! # Reading and writing of linear programs
//!
//! This module provides read and write functionality for linear program formats.
use std::convert::TryInto;
use std::fs::File;
use std::io::Read;
use std::ops::{Add, Neg};
use std::path::Path;

use num_traits::{One, Zero};
use relp_num::NonZero;
use relp_num::Rational64;

use crate::data::linear_algebra::traits::Element;
use crate::data::linear_program::general_form::GeneralForm;
use crate::io::error::{Import, Inconsistency};
use crate::io::mps::MPS;

pub mod error;
pub mod mps;

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
pub fn import<F: From<Rational64> + Zero + One + Neg<Output=F> + Ord + Element + NonZero>(
    file_path: &Path
) -> Result<impl TryInto<GeneralForm<F>, Error=Inconsistency>, Import>
where
    for<'r> F: Add<&'r F, Output=F>,
    for<'r> &'r F: Neg<Output=F>,
{
    // Open and read the file
    let mut program = String::new();
    File::open(&file_path)
        .map_err(Import::IO)?
        .read_to_string(&mut program)
        .map_err(Import::IO)?;

    // Choose the right parser
    match file_path.extension() {
        Some(extension) => match extension.to_str() {
            Some("mps" | "SIF") => mps::parse(&program).map(DataTypes::MPS),
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

enum DataTypes {
    MPS(MPS<Rational64>),
}

impl<F> TryInto<GeneralForm<F>> for DataTypes
where
    F: From<Rational64> + NonZero + Zero + One + Neg<Output=F> + Ord + Element,
    for<'r> F: Add<&'r F, Output=F>,
    for<'r> &'r F: Neg<Output=F>,
{
    type Error = Inconsistency;

    fn try_into(self) -> Result<GeneralForm<F>, Self::Error> {
        match self {
            DataTypes::MPS(mps) => TryInto::try_into(mps),
        }
    }
}
