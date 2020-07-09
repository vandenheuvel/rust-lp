//! # Parsing MPS files
//!
//! First stage of importing linear programs. Checks if the MPS file is syntactically correct, but
//! doesn't do any consistency checks for e.g. undefined row names in the column section.
use std::convert::TryFrom;

use num::{FromPrimitive, One};

use crate::data::linear_program::elements::{ConstraintType, VariableType};
use crate::io::error::{FileLocation, Parse};
use crate::io::mps::BoundType;
use crate::io::mps::RowType;
use crate::io::mps::Section;
use crate::io::mps::token::{
    COLUMN_SECTION_MARKER, COMMENT_INDICATOR, END_OF_INTEGER, START_OF_INTEGER,
};

use self::Atom::{Number, Word};

/// Most fundamental element in an MPS text file.
///
/// Every part of the input string to end up in the final `MPS` struct is parsed as either a `Word`
/// or a `Number`.
#[derive(Clone, Debug, PartialEq)]
pub(crate) enum Atom<'a, F> {
    /// A token consisting of text which should not contain whitespace.
    Word(&'a str),
    /// A token representing a number.
    Number(F),
}

/// Convert an MPS program in string form to structured data.
///
/// # Arguments
///
/// * `program`: The input string
///
/// # Return value
///
/// All lines stored in a `Vec`. The first `(usize, &str)` tuple is used for creating of errors, it
/// contains the line number and line.
pub(crate) fn into_atom_lines<F: FromPrimitive>(
    program: &impl AsRef<str>,
) -> Vec<((u64, &str), Vec<Atom<F>>)> {
    program.as_ref()
        .lines()
        .enumerate()
        .map(|(number, line)| (number as u64 + 1_u64, line))
        .filter(|(_, line)| !line.trim_start().starts_with(COMMENT_INDICATOR))
        .filter(|(_, line)| !line.is_empty())
        .map(|(number, line)| ((number, line), into_atoms(line)))
        .collect()
}

/// Convert a line into `Atom`s by testing whether an atom is a number.
///
/// # Arguments
///
/// * `line`: The input string slice
///
/// # Return value
///
/// A `Vec` of words and numbers.
fn into_atoms<F: FromPrimitive>(line: &str) -> Vec<Atom<F>> {
    line.split_whitespace()
        .map(|atom| match atom.parse().map(F::from_f64) {
            Ok(Some(value)) => Number(value),
            _ => Word(atom),
        })
        .collect()
}

/// Used to gather all data of the `MPS`.
///
/// The struct holds all `MPS` data in an intermediate parse phase. It is mostly references to the
/// original problem string, or to the values in the parsed tokens.
///
/// # Note
///
/// The information contained in this struct is not necessarily consistent, e.g. some columns
/// might reference rows which are not declared.
#[derive(Debug, PartialEq)]
pub(crate) struct UnstructuredMPS<'a, F> {
    /// Name of the linear program.
    pub name: &'a str,
    /// Name of the cost row, or objective function.
    pub cost_row_name: &'a str,
    /// Names of all constraint rows and their constraint type, excluding the cost row.
    pub rows: Vec<UnstructuredRow<'a>>,
    /// Names and type (continuous, discrete) of all variables.
    pub columns: Vec<UnstructuredColumn<'a, F>>,
    /// Right-hand sides giving a numerical value to the constraints.
    pub rhss: Vec<UnstructuredRhs<'a, F>>,
    /// Activation of constraints has a range.
    pub ranges: Vec<UnstructuredRange<'a, F>>,
    /// Bounds of the problem (does not includes elements of the "RANGES" section.
    pub bounds: Vec<UnstructuredBound<'a, F>>,
}

/// Name and constraint type of a row.
#[derive(Debug, PartialEq, Eq)]
pub(crate) struct UnstructuredRow<'a> {
    pub name: &'a str,
    pub constraint_type: ConstraintType,
}
/// Name of a column, variable type of a column, name of a row and a value (a constraint matrix
/// entry in sparse description).
#[derive(Debug, PartialEq)]
pub(crate) struct UnstructuredColumn<'a, F> {
    pub name: &'a str,
    pub variable_type: VariableType,
    pub row_name: &'a str,
    pub value: F,
}
/// Name of the right-hand side, name of the variable and constraint value.
#[derive(Debug, PartialEq)]
pub(crate) struct UnstructuredRhs<'a, F> {
    pub name: &'a str,
    pub row_name: &'a str,
    pub value: F,
}
/// Name of the range, name of the variable and range value.
///
/// Overview of how the range of a constraint is defined, depending on the constraint type:
///
/// row type | sign of r |    h    |    u
/// ---------|-----------|---------|---------
/// G  (>=)  |  + or -   |    b    | b + |r|
/// L  (<=)  |  + or -   | b - |r| |   b
/// E  (==)  |     +     |    b    | b + |r|
/// E  (==)  |     -     | b - |r| |   b
#[derive(Debug, PartialEq)]
pub(crate) struct UnstructuredRange<'a, F> {
    pub name: &'a str,
    pub row_name: &'a str,
    pub value: F,
}
/// Name of the bound, bound type, name of the variable and constraint value.
#[derive(Debug, PartialEq)]
pub(crate) struct UnstructuredBound<'a, F> {
    pub name: &'a str,
    pub bound_type: BoundType<F>,
    pub column_name: &'a str,
}

#[cfg(test)]
mod test {
    use std::convert::TryFrom;

    use crate::data::linear_program::elements::{ConstraintType, VariableType};
    use crate::io::mps::BoundType;
    use crate::io::mps::parsing::{UnstructuredColumn, UnstructuredRow};
    use crate::io::mps::parsing::Atom::*;
    use crate::io::mps::parsing::Atom;
    use crate::io::mps::parsing::into_atom_lines;
    use crate::io::mps::parsing::into_atoms;
    use crate::io::mps::parsing::parse_column_line;
    use crate::io::mps::parsing::parse_row_line;
    use crate::io::mps::RowType;
    use crate::io::mps::Section;
    use crate::io::mps::token::{COLUMN_SECTION_MARKER, END_OF_INTEGER, START_OF_INTEGER};

    #[test]
    fn test_into_atom_lines() {
        let program = "SOMETEXT 1.2\n*comment\n\n   \t line before is empty".to_string();
        let result = into_atom_lines(&program);
        let expected = vec![
            ((1, "SOMETEXT 1.2"), vec![Word("SOMETEXT"), Number(1.2f64)]),
            ((4, "   \t line before is empty"), vec![Word("line"), Word("before"), Word("is"), Word("empty")]),
        ];
        assert_eq!(result, expected);
    }

    #[test]
    fn test_into_atoms() {
        macro_rules! test {
            ($line:expr, [$($words:expr), *]) => {
                let result = into_atoms::<f64>($line);
                let expected = vec![$($words), *];
                assert_eq!(result, expected);
            }
        }

        test!("NAME TESTPROB", [Word("NAME"), Word("TESTPROB")]);
        test!("ROWS", [Word("ROWS")]);
        test!("RHS     ", [Word("RHS")]);
        test!("NUMBER 134", [Word("NUMBER"), Number(134f64)]);
        test!("NUMBER 1.6734", [Word("NUMBER"), Number(1.6734f64)]);
        test!("    MARK0000  'MARKER'                 'INTORG'",
              [Word("MARK0000"), Word("'MARKER'"), Word("'INTORG'")]);
    }

    #[test]
    fn test_parse_row_line() {
        let initial_cost_row_name = Some("COST_ROW_NAME");

        macro_rules! test_positive {
            ([$($words:expr), *], [$($expected_row:expr), *], $expected_cost_row_name:expr) => {
                let line = vec![$($words), *];
                let mut cost_row_name = initial_cost_row_name.clone();
                let mut collector = Vec::new();

                let result = parse_row_line::<f64>(line, &mut cost_row_name, &mut collector);

                assert!(result.is_ok());
                assert_eq!(cost_row_name, $expected_cost_row_name);
                assert_eq!(collector, vec![$($expected_row), *]);
            }
        }

        test_positive!([Word("E"), Word("ROW_NAME")],
                       [UnstructuredRow { name: "ROW_NAME", constraint_type: ConstraintType::Equal, }],
                       Some("COST_ROW_NAME"));
        test_positive!([Word("N"), Word("NEW_COST_ROW_NAME")],
                       [], Some("NEW_COST_ROW_NAME"));

        macro_rules! test_negative {
            ([$($words:expr), *]) => {
                let line = vec![$($words), *];
                let mut cost_row_name = initial_cost_row_name.clone();
                let mut collector = Vec::new();

                let result = parse_row_line::<f64>(line, &mut cost_row_name, &mut collector);

                assert!(result.is_err());
            }
        }

        test_negative!([Word("UNKNOWN_ROW_TYPE"), Word("ROW_NAME")]);
        test_negative!([Word("JUST_ONE_WORD")]);
        test_negative!([Word("ONE"), Word("TWO"), Word("THREE")]);
        test_negative!([Word("ONE"), Word("TWO"), Word("THREE"), Word("FOUR")]);
    }

    #[test]
    fn test_parse_column_line() {
        macro_rules! test_positive {
            (
                [$($words:expr), *],
                [$($expected_data:expr, ) *],
                $initial_marker:expr,
                $expected_marker:expr
            ) => {
                let line = vec![$($words), *];
                let mut marker = $initial_marker;
                let mut collector = Vec::new();

                let result = parse_column_line::<i32>(line, &mut marker, &mut collector);

                assert!(result.is_ok());
                assert_eq!(collector, vec![$($expected_data), *]);
                assert_eq!(marker, $expected_marker);
            }
        }

        test_positive!([Word("CNAME"), Word("RNAME"), Number(5)],
            [
                UnstructuredColumn {
                    name: "CNAME",
                    variable_type: VariableType::Continuous,
                    row_name: "RNAME",
                     value: 5,
                 },
            ],
            VariableType::Continuous, VariableType::Continuous);
        test_positive!([Word("CNAME"), Word("RNAME"), Number(5)],
            [
                UnstructuredColumn {
                    name: "CNAME",
                    variable_type: VariableType::Integer,
                    row_name: "RNAME",
                     value: 5,
                },
            ],
            VariableType::Integer, VariableType::Integer);
        test_positive!([Word("CNAME1"), Word("RNAME1"), Number(1), Word("RNAME2"), Number(2)],
            [
                UnstructuredColumn {
                    name: "CNAME1",
                    variable_type: VariableType::Continuous,
                    row_name: "RNAME1",
                     value: 1,
                 },
                 UnstructuredColumn {
                     name: "CNAME1",
                     variable_type: VariableType::Continuous,
                     row_name: "RNAME2",
                      value: 2,
                  },
            ],
            VariableType::Continuous, VariableType::Continuous);
        test_positive!([Word("MARKER_NAME"), Word(COLUMN_SECTION_MARKER), Word(START_OF_INTEGER)],
            [], VariableType::Continuous, VariableType::Integer);
        test_positive!([Word("MARKER_NAME"), Word(COLUMN_SECTION_MARKER), Word(END_OF_INTEGER)],
            [], VariableType::Integer, VariableType::Continuous);
    }

    #[test]
    fn try_from_section() {
        macro_rules! test {
            ([$($words:expr), *], $expected:expr) => {
                let input: Vec<Atom<f64>> = vec![$($words), *];
                let result = Section::try_from(&input);
                assert_eq!(result, $expected);
            }
        }

        test!([Word("NAME"), Word("THENAME")], Ok(Section::Name("THENAME")));
        test!([Word("ROWS")], Ok(Section::Rows));
        test!([Word("COLUMNS")], Ok(Section::Columns(VariableType::Continuous)));
        test!([Word("RHS")], Ok(Section::Rhs));
        test!([Word("BOUNDS")], Ok(Section::Bounds));
        test!([Word("RANGES")], Ok(Section::Ranges));
        test!([Word("ENDATA")], Ok(Section::Endata));
        test!([Word("X")], Err(()));
        test!([Number(1.556f64)], Err(()));
        test!([], Err(()));
    }

    #[test]
    fn try_from_row_type() {
        macro_rules! test_positive {
            ($word:expr, $expected:expr) => {
                let result = RowType::try_from($word);
                assert!(result.is_ok());
                assert_eq!(result.unwrap(), $expected);
            };
        }
        test_positive!("N", RowType::Cost);
        test_positive!("L", RowType::Constraint(ConstraintType::Less));
        test_positive!("E", RowType::Constraint(ConstraintType::Equal));
        test_positive!("G", RowType::Constraint(ConstraintType::Greater));

        macro_rules! test_negative {
            ($word:expr) => {
                let result = RowType::try_from($word);
                assert!(result.is_err());
            };
        }
        test_negative!("X");
        test_negative!("");
        test_negative!("\t");
    }

    #[test]
    fn try_from_bound_type() {
        let result: Result<BoundType<f64>, _> = BoundType::try_from("FR");
        assert_eq!(result, Ok(BoundType::Free));

        macro_rules! test {
            ($word:expr, $expected:ident) => {
                let result: Result<BoundType<f64>, _> = BoundType::try_from($word);
                assert_eq!(result, Ok(BoundType::$expected));
            };
        }

        test!("MI", LowerMinusInfinity);
        test!("PL", UpperInfinity);
        test!("BV", Binary);
        assert!(BoundType::<f64>::try_from("X").is_err());
        assert!(BoundType::<f64>::try_from("").is_err());
        assert!(BoundType::<i32>::try_from("\t").is_err());

        macro_rules! test {
            ($word:expr, $expected:ident) => {
                let result = BoundType::try_from(($word, 4));
                assert_eq!(result, Ok(BoundType::$expected(4)));
            };
        }
        test!("LO", LowerContinuous);
        test!("UP", UpperContinuous);
        test!("FX", Fixed);
        test!("LI", LowerInteger);
        test!("UI", UpperInteger);
    }
}
