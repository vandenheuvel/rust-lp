use std::ops::Range;

use crate::io::error::Import as ImportError;
use crate::io::error::Parse as ParseError;
use crate::io::error::ParseResult;
use crate::io::mps::MPS;
use crate::io::mps::number::parse::Parse;
use crate::io::mps::parse::{ColumnLineContent, ColumnRetriever};
use crate::io::mps::parse::free::Free;
use crate::io::mps::token::COLUMN_SECTION_MARKER;

/// Read a linear program with the "fixed columns" parsing strategy.
pub fn parse<N: Parse>(program_string: &str) -> Result<MPS<N>, ImportError> {
    super::parse::<_, Fixed>(program_string)
}

pub(super) struct Fixed;
impl<'a> ColumnRetriever<'a> for Fixed {
    type RestOfLine = &'a str;

    fn two_or_three(line_after_name: &str) -> ParseResult<[&str; 1]> {
        // Defer to the flexible method, because some implementations print the program name in the
        // third column (rather than the second, as the lp solve authors specify). It's called only
        // once per MPS file anyway.
        Free::two_or_three(line_after_name)
    }

    fn one_and_two(line: &str) -> ParseResult<[&str; 2]> {
        if line.len() > FIELDS[2].start {
            let type_text = &line[FIELDS[1].clone()];

            let name_field_end = FIELDS[2].end.min(line.len());
            let name = line[FIELDS[2].start..name_field_end].trim_end();
            if name.is_empty() {
                Err(ParseError::new("Empty row name."))
            } else {
                Ok([type_text, name])
            }
        } else {
            Err(ParseError::new("Line is too short."))
        }
    }

    fn is_column_marker_line(line: &'a str) -> ParseResult<ColumnLineContent<Self::RestOfLine>> {
        if line.len() >= FIELDS[4].end {
            if &line[FIELDS[3].clone()] == COLUMN_SECTION_MARKER {
                if line.len() >= FIELDS[5].end {
                    Ok(ColumnLineContent::Marker([&line[FIELDS[5].clone()]]))
                } else {
                    Err(ParseError::new("Line is too short to be a marker line."))
                }
            } else {
                Ok(ColumnLineContent::Data([
                    &line[FIELDS[2].clone()].trim_end(),
                    &line[FIELDS[3].clone()].trim_end(),
                    &line[FIELDS[4].clone()].trim_start(),
                ], &line[FIELDS[4].end..]))
            }
        } else {
            Err(ParseError::new("Line is too short."))
        }
    }

    fn two_through_four(line: &'a str) -> ParseResult<([&str; 3], Self::RestOfLine)> {
        if line.len() >= FIELDS[4].end {
            Ok(([
                line[FIELDS[2].clone()].trim_end(),
                line[FIELDS[3].clone()].trim_end(),
                line[FIELDS[4].clone()].trim_start(),
            ], &line[FIELDS[4].end..]))
        } else {
            Err(ParseError::new("Line is too short."))
        }
    }

    fn five_and_six(line_after_first_four: Self::RestOfLine) -> ParseResult<Option<[&'a str; 2]>> {
        let relative_field_six_end = FIELDS[6].end - FIELDS[4].end;
        let result = if line_after_first_four.len() >= relative_field_six_end {
            let (five_start, five_end) = (FIELDS[5].start - FIELDS[4].end, FIELDS[5].end - FIELDS[4].end);
            let (six_start, six_end) = (FIELDS[6].start - FIELDS[4].end, FIELDS[6].end - FIELDS[4].end);

            let five = line_after_first_four[five_start..five_end].trim_end();
            let six = line_after_first_four[six_start..six_end].trim_start();
            if !five.is_empty() && !six.is_empty() {
                Some([
                    five,
                    six,
                ])
            } else { None }
        } else { None };

        Ok(result)
    }

    fn one_through_three(line: &'a str) -> ParseResult<([&str; 3], Self::RestOfLine)> {
        if line.len() >= FIELDS[3].start {
            let type_text = &line[FIELDS[1].clone()];
            let name = line[FIELDS[2].clone()].trim_end();
            let column = line[FIELDS[3].clone()].trim_end();

            let rest = &line[FIELDS[3].end..];

            Ok(([type_text, name, column], rest))
        } else {
            Err(ParseError::new("Line is too short."))
        }
    }

    fn four(rest_of_line: Self::RestOfLine) -> ParseResult<[&'a str; 1]> {
        let end_of_column = FIELDS[4].end - FIELDS[3].end;
        if rest_of_line.len() >= end_of_column {
            let start_of_column = FIELDS[4].start - FIELDS[3].end;
            let value_text = rest_of_line[start_of_column..end_of_column].trim_start();
            Ok([value_text])
        } else {
            Err(ParseError::new("Line doesn't contain another value, it's too short."))
        }
    }
}

/// Character ranges for the different fields of the column.
const FIELDS: [Range<usize>; 7] = [
    0..1,
    1..3,
    4..12,
    14..22,
    24..36,
    39..47,
    49..61,
];
