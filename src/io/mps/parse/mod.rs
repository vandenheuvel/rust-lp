//! # Parsing MPS files
//!
//! First stage of importing linear programs. Checks if the MPS file is syntactically correct, but
//! doesn't do any consistency checks for e.g. undefined row names in the column section.
use std::collections::{HashMap, HashSet};
use std::mem::take;
use std::str::FromStr;

use crate::data::linear_algebra::{SparseTuple, SparseTupleVec};
use crate::data::linear_program::elements::{ConstraintRelation, Objective, VariableType};
use crate::io::error::{FileLocation, Import as ImportError, Inconsistency, Parse as ParseError, ParseResult};
use crate::io::mps::{Bound, BoundType, Column, MPS, Range, Rhs, Row};
use crate::io::mps::number::parse::Parse;
use crate::io::mps::RowType;
use crate::io::mps::Section;
use crate::io::mps::token::{COMMENT_INDICATOR, END_OF_INTEGER, NAME, START_OF_INTEGER};

pub(crate) mod fixed;
pub(crate) mod free;

// TODO(ARCHITECTURE): Consider moving all methods in a trait to reduce the number of type
//  annotations. Might need impl trait in trait definitions?

/// Parse an MPS file in either the fixed or free format.
///
/// This method is written to be fast and is not "defensive" at all. E.g. when there is an extra
/// value on a line, this value will just be ignored.
///
/// # Arguments
///
/// * `program_string`: String holding the entire program.
///
/// # Errors
///
/// Parse Errors if the file format is found out to be wrong, Inconsistency errors if e.g. a row is
/// mentioned that wasn't declared in advance.
fn parse<'a, F: Parse, CR: ColumnRetriever<'a>>(
    program_string: &'a str,
) -> Result<MPS<F>, ImportError> {
    // Go through the file line by line
    let mut lines = into_lines(program_string);

    // Extract the program name
    let name = parse_program_name::<CR>(lines.next())
        .map_err(|e| e.wrap("Error while parsing reading the program name."))?;

    // Go through the file, line by line, section by section. We never look back.
    // The below call removes the line that indicates the start of the ROWS section.
    let objective = parse_objsense_section(&mut lines)?;

    let (cost_row_name, rows) = parse_row_section::<CR, _>(&mut lines)?;
    let (cost_row_name, rows) = check_row_section_consistency(cost_row_name, rows)?;
    let row_index = build_row_index(&rows, &cost_row_name)?;

    // All sections after this one are optional, so the next section is read by the previous method.

    let (columns, cost_values, next_section) =
        parse_column_section::<_, CR, _>(&mut lines, &cost_row_name, &row_index)?;
    let column_index = build_column_index(&columns);

    let (rhss, next_section) = if next_section == Section::Rhs {
        parse_optional_section::<_, CR, _, _, 2>(&mut lines, [Section::Ranges, Section::Bounds], &row_index)?
    } else { (Vec::with_capacity(0), next_section) };

    let (ranges, next_section) = if next_section == Section::Ranges {
        parse_optional_section::<_, CR, _, _, 1>(&mut lines, [Section::Bounds], &row_index)?
    } else { (Vec::with_capacity(0), next_section) };
    check_ranges_consistency(&ranges)?;

    let bounds = if next_section == Section::Bounds {
        parse_bounds_section::<_, CR, _>(&mut lines, &column_index)?
    } else { Vec::with_capacity(0) };

    if lines.next().is_some() {
        return Err(ParseError::new("File parsed successfully, but it has nonempty lines a\
        t the end").into());
    }

    Ok(MPS {
        name,
        objective,
        cost_row_name,
        cost_values,
        rows,
        columns,
        rhss,
        ranges,
        bounds,
    })
}

/// Split a linear program into numbered lines, skipping comments.
///
/// # Arguments
///
/// * `text`: The program string.
///
/// # Return value
///
/// An iterator over numbered lines.
fn into_lines<'a>(text: &'a str) -> impl Iterator<Item = FileLocation> + 'a {
    text.lines()
        .enumerate()
        .map(|(number, line)| (number + 1, line)) // Count from 1
        .filter(|(_, line)| !line.trim_start().starts_with(COMMENT_INDICATOR))
        .filter(|(_, line)| !line.is_empty())
}

trait ColumnRetriever<'a> {
    /// Value representing a line that is partially parsed.
    ///
    /// When some values might need to be read from the same line later, but we don't want to start
    /// searching from scratch, this type can be used to represent the point at the line where the
    /// searching can continue.
    type RestOfLine;

    /// Read the program name from the second or third column on the line.
    ///
    /// # Arguments
    ///
    /// * `line_after_first_column`: Rest of the line after the first column was removed. That is,
    /// this method should essentially return the first column relative to the input.
    fn two_or_three(line_after_name: &str) -> ParseResult<[&str; 1]>;

    /// Read the first and second column.
    ///
    /// These are the row type and name of the row.
    fn one_and_two(line: &str) -> ParseResult<[&str; 2]>;

    /// Check whether this is a column marker line.
    ///
    /// Depending on whether it is, either the marker field is returned (marker name is discarded)
    /// or the column name, row name and value text are returned (typically happens).
    /// TODO(OPTIMIZATION): Should the implementors of this method use branch prediction hints?
    fn is_column_marker_line(line: &'a str) -> ParseResult<ColumnLineContent<Self::RestOfLine>>;

    /// Read three columns from either the RHS or RANGES section.
    ///
    /// # Return value
    ///
    /// If successful, the group name, row name and value text, as well as the rest of the line.
    fn two_through_four(line: &'a str) -> ParseResult<([&str; 3], Self::RestOfLine)>;

    /// Try to read two more columns.
    ///
    /// TODO(CORRECTNESS): Let this fail when there is e.g. only one value?
    fn five_and_six(line_after_first_four: Self::RestOfLine) -> ParseResult<Option<[&'a str; 2]>>;

    fn one_through_three(line: &'a str) -> ParseResult<([&str; 3], Self::RestOfLine)>;

    /// While parsing bounds: if the bound type demands it, try to read one more field.
    fn four(rest_of_line: Self::RestOfLine) -> ParseResult<[&'a str; 1]>;
}
enum ColumnLineContent<'a, ROL> {
    Marker([&'a str; 1]),
    Data([&'a str; 3], ROL),
}

/// Read the program name.
///
/// # Arguments
///
/// `location`: First line of the program that should exist.
///
/// # Errors
///
/// If the first line of the program is not there, when the line doesn't start with the `NAME`
/// indicator, or when the name could not be read.
fn parse_program_name<'a, CR: ColumnRetriever<'a>>(
    location: Option<FileLocation<'a>>,
) -> ParseResult<String> {
    match location {
        None => Err(ParseError::new("No line to read, is the file empty?")),
        Some((number, line)) => {
            if line.len() < NAME.len() {
                return Err(ParseError::with_location("Line too short.", (number, line)));
            }

            let expected_name_string = &line[0..NAME.len()];
            if expected_name_string == NAME {
                let name = CR::two_or_three(&line[NAME.len()..])
                    .map_err(|e| e.wrap(format!("Could not read program name \
                    from line {}. Ths line looked like \"{}\"", number, line)))?[0];
                Ok(name.to_string())
            } else {
                Err(ParseError::with_location(
                    format!(
                        "Expected a \"{}\" indicator, found \"{}\" instead", NAME, expected_name_string,
                    ),
                    (number, line),
                ))
            }
        },
    }
}

/// Try to read the objective direction and consume the rows section indicator.
///
/// # Arguments
///
/// `lines`: Iterator yielding meaningful program lines.
///
/// # Errors
///
/// If the iterator is empty, or the next lines don't pretty much look like
/// ```compile_fail
/// OBJSENSE
///   MAXIMIZE
/// ```
/// where it is allowed to abbreviate the direction to the first three characters.
fn parse_objsense_section<'a, L: Iterator<Item = FileLocation<'a>>>(
    lines: &mut L,
) -> ParseResult<Objective> {
    match lines.next() {
        None => Err(ParseError::new("No line to read, is the program more than a name?")),
        Some((number, line)) => {
            match line.trim_end() {
                // The default: OBJSENSE is not mentioned and we minimize.
                "ROWS" => Ok(Objective::Minimize),
                "OBJSENSE" => {
                    match lines.next() {
                        None => Err(ParseError::new("Program can't end in the OBJSENSE section.")),
                        Some((number, line)) => {
                            // Throw away the next line, who's contents are "ROWS", check to be sure
                            if !matches!(lines.next(), Some((_, line)) if line.starts_with("ROWS")) {
                                return Err(ParseError::new("Expected the ROWS section next."));
                            }

                            match line.trim_end() {
                                "  MINIMIZE" | "  MIN" => Ok(Objective::Minimize),
                                "  MAXIMIZE" | "  MAX" => Ok(Objective::Maximize),
                                _ => Err(ParseError::with_location(
                                    format!("Can't read objective {}", line.trim_end()),
                                    (number, line),
                                )),
                            }
                        }
                    }
                },
                _ => Err(ParseError::with_location(
                    format!("Line contents \"{}\" were unexpected", line),
                    (number, line),
                )),
            }
        }
    }
}

fn parse_row_section<'a, CR: ColumnRetriever<'a>, L: Iterator<Item = FileLocation<'a>>>(
    lines: &mut L,
) -> ParseResult<(Option<String>, Vec<Row>)> {
    let mut collector = Vec::new();
    let mut cost_row_name = None;

    for (number, line) in lines {
        if is_part_of_same_section(line) {
            // Only indicators for new sections start at the first character, so this row must be
            // part of same section.
            let [row_type, name] = CR::one_and_two(line)?;
            let name = name.to_string();
            match RowType::from_str(row_type) {
                Ok(row_type) => {
                    // If this is the cost row, save it separately.
                    match row_type {
                        RowType::Cost => {
                            if cost_row_name.is_some() {
                                return Err(ParseError::with_location("Second cost row \
                                detected. This is not supported.", (number, line)));
                            }
                            cost_row_name = Some(name)
                        },
                        RowType::Constraint(constraint_type) => {
                            collector.push(Row { name, constraint_type, });
                        },
                    }
                },
                Err(error) => return Err(error.wrap(
                    format!("Couldn't parse row type on line {}: \"{}\"", number, line),
                )),
            }
        } else {
            // Section has ended.
            try_parse_next_section(line, [Section::Columns])?;
            return Ok((cost_row_name, collector));
        }
    }

    Err(ParseError::new("Section ended sooner than expected."))
}

fn check_row_section_consistency(
    cost_row_name: Option<String>,
    mut rows: Vec<Row>,
) -> Result<(String, Vec<Row>), Inconsistency> {
    if let Some(cost_name) = cost_row_name {
        rows.sort_unstable_by(|row1, row2| row1.name.cmp(&row2.name));
        if rows.binary_search_by_key(&cost_name.as_str(), |row| &row.name.as_str()).is_ok() {
            return Err(Inconsistency::new("Cost row name found in other rows."));
        }

        if let Some(name) = rows.windows(2)
            .find(|w| w[0].name == w[1].name).map(|w| &w[0].name) {
            return Err(Inconsistency::new(format!("Duplicate row name {} found.", name)));
        }

        Ok((cost_name, rows))
    } else {
        Err(Inconsistency::new("No cost name read."))
    }
}

/// Extract all row names, and assign to each row a fixed row index.
///
/// This index will be used throughout building the `MPS`, and ordered as in the original problem.
/// It is needed because constraint values for columns might appear in an unsorted order, so then a
/// lookup is practical.
///
/// # Arguments
///
/// * `unstructured_rows`: Collection of unstructured rows, cost row should not be in here.
///
/// # Return value
///
/// Names of all rows (this excludes the cost row) and a map from the name to the index at which
/// they are stored. This arbitrary order is the order in which rows will be stored in the final
/// MPS data structure.
///
/// # Errors
///
/// If there are rows with duplicate names, the MPS is considered inconsistent.
fn build_row_index<'a>(
    rows: &'a [Row],
    cost_row_name: &str,
) -> Result<HashMap<&'a str, usize>, Inconsistency> {
    let names = rows.iter().map(|row| row.name.as_str()).collect::<HashSet<_>>();
    debug_assert_eq!(names.len(), rows.len());
    debug_assert!(!names.contains(cost_row_name));

    let index = rows.iter().enumerate()
        .map(|(i, row)| (row.name.as_str(), i))
        .collect::<HashMap<_, _>>();

    if index.values().len() < rows.len() {
        Err(Inconsistency::new("Not all row names are unique"))
    } else {
        Ok(index)
    }
}

fn parse_column_section<'a, F: Parse, CR: ColumnRetriever<'a>, L: Iterator<Item = FileLocation<'a>>>(
    lines: &mut L,
    cost_row_name: &str,
    row_map: &HashMap<&str, usize>,
) -> Result<(Vec<Column<F>>, Vec<(usize, F)>, Section), ImportError> {
    // Collecting data that is returned from the function.
    let mut collector = Vec::new();
    // Collecting data per column. Data should by sorted by column, that makes this easier.
    let mut column_collector = Vec::new();
    // Collecting cost vector tuples.
    let mut cost_values_collector = Vec::new();
    // Column currently being read. If this value is a `Some`, the values in `column_collector`
    // still need to be saved in `collector`.
    let mut column_to_be_saved = None;
    // We switch back and forth with the variable type as the relevant marker is encountered.
    let mut active_variable_type = VariableType::Continuous;

    for (number, line) in lines {
        if is_part_of_same_section(line) {
            // Row is part of same section

            let content = CR::is_column_marker_line(line)
                .map_err(|e| e.wrap(format!(
                    "Could not determine whether line {} is a marker line: \"{}\"", number, line,
                )))?;
            match content {
                ColumnLineContent::Marker([marker_text]) => {
                    // A marker indicates a variable type change.

                    // When the variable type changes, a new column comes. Save before changing the
                    // variable type.
                    save_to_column_collector(
                        None,
                        &mut column_to_be_saved, &mut column_collector, &mut collector, active_variable_type,
                    )?;

                    // TODO(CORRECTNESS): Err when the type doesn't actually change?
                    let previous_variable_type = active_variable_type;
                    active_variable_type = match marker_text {
                        START_OF_INTEGER => VariableType::Integer,
                        END_OF_INTEGER => VariableType::Continuous,
                        _ => return Err(ParseError::with_location(
                            format!("Marker type \"{}\" unknown.", marker_text), (number, line),
                        ).into()),
                    };
                    debug_assert_ne!(previous_variable_type, active_variable_type);
                },
                ColumnLineContent::Data([column_name, row_name, value_text], rest_of_line) => {
                    // Read the name of the group and save the group if necessary
                    if let Some(active_group) = column_to_be_saved.as_ref() {
                        if active_group != column_name {
                            save_to_column_collector(
                                Some(column_name.to_string()),
                                &mut column_to_be_saved, &mut column_collector, &mut collector, active_variable_type,
                            )?;
                        }
                    } else {
                        column_to_be_saved = Some(column_name.to_string());
                    }

                    let mut save_pair = |row_name, value_text| {
                        let value = F::parse(value_text).map_err(|e| e.wrap(
                            "Couldn't parse (row name, value) pair",
                        ))?;

                        match row_map.get(row_name) {
                            None => if row_name == cost_row_name {
                                let column_index = collector.len();
                                cost_values_collector.push((column_index, value))
                            } else {
                                return Err(Inconsistency::new(format!("Row \"{}\" not known.", row_name)).into());
                            },
                            Some(row_index) => column_collector.push((*row_index, value)),
                        }

                        Ok(())
                    };

                    let _result: Result<_, ImportError> = save_pair(row_name, value_text);
                    _result?;
                    match CR::five_and_six(rest_of_line) {
                        Ok(Some([row_name, value_text])) => save_pair(row_name, value_text)?,
                        Ok(None) => {},
                        Err(e) => return Err(e.wrap(format!(
                            "Line {} contained an unexpected number of elements: \"{}\"", number, line,
                        )).into()),
                    }
                },
            }
        } else {
            let next_section = try_parse_next_section(
                line,
                [Section::Rhs, Section::Ranges, Section::Bounds],
            )?;
            // Save the values read of the last column before quiting.
            save_to_column_collector(
                None,
                &mut column_to_be_saved, &mut column_collector, &mut collector, active_variable_type,
            )?;

            debug_assert!(cost_values_collector.is_sorted_by_key(|&(j, _)| j));
            return Ok((collector, cost_values_collector, next_section));
        }
    }

    Err(ParseError::new("Section ended sooner than expected.").into())
}

/// Every time an INTEND or INTORG marker is encountered, or when a new column name is encountered,
/// the (row name, value) tuples collected in the `column_collector` are `saved` in the "over-all"
/// collector `collector` that is returned from the function.
fn save_to_column_collector<F>(
    new_column: Option<String>,
    column_to_be_saved: &mut Option<String>,
    column_collector: &mut SparseTupleVec<F>,
    collector: &mut Vec<Column<F>>,
    active_variable_type: VariableType,
) -> Result<(), Inconsistency> {
    if let Some(column_name) = column_to_be_saved.take() {
        let mut values = take(column_collector);
        values.sort_unstable_by_key(|&(i, _)| i);
        if values.windows(2).any(|w| w[0].0 == w[1].0) {
            return Err(Inconsistency::new(format!("Duplicate row for column \"{}\"", column_name)));
        }

        collector.push(Column {
            name: column_name,
            variable_type: active_variable_type,
            values,
        });
    }

    *column_to_be_saved = new_column;

    Ok(())
}

/// Assign to each column a fixed column index.
///
/// This index will be used to organize the bounds information.
///
/// # Arguments
///
/// * `columns`: Collection of the columns
///
/// # Return value
///
/// A map assigning to each column name a value.
///
/// # Note
///
/// This assignment is not a specific order.
fn build_column_index<F>(columns: &[Column<F>]) -> HashMap<&str, usize> {
    debug_assert_eq!(
        columns.iter().map(|column| column.name.as_str()).collect::<HashSet<_>>().len(),
        columns.len(),
    );

    columns.iter()
        .enumerate()
        .map(|(index, column)| (column.name.as_str(), index))
        .collect()
}

trait ListedInGroup {
    type ValueType;
    fn new(name: String, values: Vec<(usize, Self::ValueType)>) -> Self;
}
impl<F> ListedInGroup for Rhs<F> {
    type ValueType = F;

    fn new(name: String, values: Vec<(usize, Self::ValueType)>) -> Self {
        Self { name, values, }
    }
}
impl<F> ListedInGroup for Range<F> {
    type ValueType = F;

    fn new(name: String, values: Vec<(usize, Self::ValueType)>) -> Self {
        Self { name, values, }
    }
}
impl<F> ListedInGroup for Bound<F> {
    type ValueType = BoundType<F>;

    fn new(name: String, values: Vec<(usize, Self::ValueType)>) -> Self {
        Self { name, values, }
    }
}

fn parse_optional_section<
    'a,
    F: Parse,
    CR: ColumnRetriever<'a>,
    L: Iterator<Item = FileLocation<'a>>,
    T: ListedInGroup<ValueType = F>,
    const NRVNS: usize,
>(
    lines: &mut L,
    valid_next_sections: [Section; NRVNS],
    row_map: &HashMap<&str, usize>,
) -> Result<(Vec<T>, Section), ImportError> {
    let mut collector = Vec::new();
    let mut group_collector = Vec::new();
    let mut group_to_be_saved = None;

    for (_number, line) in lines {
        if is_part_of_same_section(line) {
            // Row is part of same section
            let ([group_name, row_name, value_text], rest_of_line) = CR::two_through_four(line)?;
            parse_value_section_line::<_, _, CR>(
                group_name,
                row_name,
                value_text,
                rest_of_line,
                &mut group_to_be_saved,
                &mut group_collector,
                &mut collector,
                row_map,
            )?;
        } else {
            let next_section = try_parse_next_section(line, valid_next_sections)?;
            save_to_group_collector::<_, false>(
                None,
                &mut group_to_be_saved, &mut group_collector, &mut collector,
            )?;
            return Ok((collector, next_section))
        }
    }

    Err(ParseError::new("Section \"COLUMNS\" ended sooner than expected.").into())
}

fn parse_value_section_line<'a, F: Parse, T: ListedInGroup<ValueType = F>, CR: ColumnRetriever<'a>>(
    group_name: &str,
    row_name: &str,
    value_text: &str,
    rest_of_line: CR::RestOfLine,
    group_to_be_saved: &mut Option<String>,
    group_collector: &mut Vec<(usize, F)>,
    collector: &mut Vec<T>,
    row_map: &HashMap<&str, usize>,
) -> Result<(), ImportError> {
    // Read the name of the group and save the group if necessary
    if let Some(active_group) = &*group_to_be_saved {
        if active_group != group_name {
            save_to_group_collector::<_, false>(
                Some(group_name.to_string()),
                group_to_be_saved, group_collector, collector,
            )?;
        }
    } else {
        *group_to_be_saved = Some(group_name.to_string());
    }

    let mut save_pair = |row_name, value_text| {
        let row_index: Result<_, Inconsistency> = row_map.get(row_name)
            .ok_or_else(|| Inconsistency::new(format!("Row \"{}\" not known.", row_name)));
        let value: Result<_, ParseError> = F::parse(value_text).map_err(|e| e.wrap(
            "Couldn't parse (row name, value) pair",
        ));

        group_collector.push((*row_index?, value?));
        Ok(())
    };

    save_pair(row_name, value_text)?;
    match CR::five_and_six(rest_of_line) {
        Ok(Some([row_name, value_text])) => save_pair(row_name, value_text),
        Ok(None) => Ok(()),
        Err(e) => Err(e.wrap("Line contained an unexpected number of elements").into()),
    }
}

fn save_to_group_collector<T: ListedInGroup, const CAN_HAVE_DUPLICATES: bool>(
    new_group: Option<String>,
    group_to_be_saved: &mut Option<String>,
    group_collector: &mut SparseTupleVec<T::ValueType>,
    collector: &mut Vec<T>,
) -> Result<(), Inconsistency> {
    if let Some(group_name) = group_to_be_saved.take() {
        let mut values = take(group_collector);
        values.sort_unstable_by_key(|&(i, _)| i);
        // TODO(CORRECTNESS): The level of defensiveness should be more consistent: why check for
        //  duplicates here?
        if !CAN_HAVE_DUPLICATES {
            if let Some(&[(row_id, _), _]) = values.windows(2)
                .find(|w| w[0].0 == w[1].0) {
                return Err(Inconsistency::new(format!(
                    "Duplicate row id \"{}\" for group \"{}\"", row_id, group_name)));
            }
        }

        collector.push(T::new(group_name, values));
    }

    *group_to_be_saved = new_group;

    Ok(())
}

fn check_ranges_consistency<F>(ranges: &[Range<F>]) -> Result<(), Inconsistency> {
    let mut unique = HashSet::with_capacity(ranges.iter().map(|r| r.values.len()).sum());
    let all_unique = ranges.iter().flat_map(|range| range.values.iter())
        .all(|&(row_index, _)| unique.insert(row_index));

    if all_unique { Ok(()) } else {
        Err(Inconsistency::new("Each row can have at most one range value"))
    }
}

fn parse_bounds_section<'a, F: Parse, CR: ColumnRetriever<'a>, L: Iterator<Item = FileLocation<'a>>>(
    lines: &mut L,
    column_index: &HashMap<&str, usize>,
) -> Result<Vec<Bound<F>>, ImportError> {
    let mut collector = Vec::new();
    let mut bound_collector = Vec::new();
    let mut bound_to_be_saved = None;

    for (_number, line) in lines {
        if is_part_of_same_section(line) {
            // Row is part of same section
            let ([bound_type_text, bound_name, column_name], rest_of_line) = CR::one_through_three(line)?;
            parse_bound_line::<_, CR>(
                bound_type_text,
                bound_name,
                column_name,
                rest_of_line,
                &mut bound_to_be_saved,
                &mut bound_collector,
                &mut collector,
                column_index,
            )?;
        } else {
            try_parse_next_section(line, [])?;
            save_to_group_collector::<_, true>(
                None,
                &mut bound_to_be_saved, &mut bound_collector, &mut collector,
            )?;
            return Ok(collector);
        }
    }

    Err(ParseError::new("Section \"COLUMNS\" ended sooner than expected.").into())
}

fn parse_bound_line<'a, F: Parse, CR: ColumnRetriever<'a>>(
    bound_type_text: &str,
    bound_name: &str,
    column_name: &str,
    rest_of_line: CR::RestOfLine,
    active_group: &mut Option<String>,
    group_collector: &mut Vec<SparseTuple<BoundType<F>>>,
    collector: &mut Vec<Bound<F>>,
    column_index: &HashMap<&str, usize>,
) -> Result<(), ImportError> {
    let column_index = *column_index.get(column_name).ok_or_else(
        || Inconsistency::new(format!("Column name \"{}\" unknown", column_name))
    )?;

    if let Some(active_column) = active_group {
        if active_column != bound_name {
            save_to_group_collector::<_, true>(
                Some(bound_name.to_string()),
                active_group, group_collector, collector,
            )?;
        }
    } else {
        *active_group = Some(bound_name.to_string());
    }

    let bound_type = match bound_type_text {
        "FR" => BoundType::Free,
        "MI" => BoundType::LowerMinusInfinity,
        "PL" => BoundType::UpperInfinity,
        "BV" => BoundType::Binary,
        "LO" | "UP" | "FX" | "LI" | "UI" => {
            let [value_text] = CR::four(rest_of_line)?;
            let value = F::parse(value_text).map_err(|e| e.wrap(
                "Couldn't parse (row name, value) pair",
            ))?;

            match bound_type_text {
                "LO" => BoundType::LowerContinuous(value),
                "UP" => BoundType::UpperContinuous(value),
                "FX" => BoundType::Fixed(value),
                "LI" => BoundType::LowerInteger(value),
                "UI" => BoundType::UpperInteger(value),
                _ => unreachable!(),
            }
        },
        "SC" => unimplemented!(),
        _ => return Err(ParseError::new(format!("Bound type \"{}\" unknown.", bound_type_text)).into()),
    };

    group_collector.push((column_index, bound_type));

    Ok(())
}

fn try_parse_next_section<const N: usize>(
    line: &str,
    acceptable_next_sections: [Section; N],
) -> ParseResult<Section> {
    let new_section = Section::from_str(line)
        .map_err(|e| e.wrap(format!("Could not parse new section header \"{}\".", line)))?;
    if new_section != Section::Endata && !acceptable_next_sections.contains(&new_section) {
        return Err(ParseError::new(format!(
            "Expected the {} section headers, found the {} section.", Section::Endata, new_section,
        )));
    }

    Ok(new_section)
}

fn is_part_of_same_section(line: &str) -> bool {
    debug_assert_ne!(line.len(), 0);

    line.starts_with(' ')
}

impl<'a> FromStr for Section {
    type Err = ParseError;

    /// Try to read a `Section` from a `Vec` slice of `Atom`s.
    ///
    /// # Arguments
    ///
    /// * `line`: The input line consisting of a sequence of `Atom`s.
    ///
    /// # Return value
    ///
    /// A `Section` variant describing the section this line announces, if one is recognized.
    ///
    /// # Errors
    ///
    /// A `()` error if no `Section` is recognized.
    fn from_str(text: &str) -> Result<Self, Self::Err> {
        match text {
            "ROWS"     => Ok(Section::Rows),
            "COLUMNS"  => Ok(Section::Columns),
            "RHS"      => Ok(Section::Rhs),
            "BOUNDS"   => Ok(Section::Bounds),
            "RANGES"   => Ok(Section::Ranges),
            "ENDATA"   => Ok(Section::Endata),
            _  => Err(ParseError::new(format!("Unknown section header \"{}\".", text))),
        }
    }
}

impl FromStr for RowType {
    type Err = ParseError;

    /// Try to read a `RowType` from a string slice.
    ///
    /// The type of a row is denoted by `N` if it's the cost row; this row is often unique. There is
    /// no defined behaviour for multiple cost rows. Constraint rows are indicated by `L`, `E` or
    /// `G`.
    ///
    /// # Arguments
    ///
    /// * `word`: The input `String` slice.
    ///
    /// # Return value
    ///
    /// A `RowType` variant if the `String` slice matches either `N`, `L`, `E` or `G`.
    ///
    /// # Errors
    ///
    /// Any `String` slices not equal to either `N`, `L`, `E` or `G` will fair to be parsed.
    fn from_str(word: &str) -> Result<RowType, Self::Err> {
        match &word[0..1] {
            "N" => Ok(RowType::Cost),
            "L" => Ok(RowType::Constraint(ConstraintRelation::Less)),
            "E" => Ok(RowType::Constraint(ConstraintRelation::Equal)),
            "G" => Ok(RowType::Constraint(ConstraintRelation::Greater)),
            _ => Err(ParseError::new(format!("Row type \"{}\" unknown.", word))),
        }
    }
}

#[cfg(test)]
mod test {
    use crate::data::linear_program::elements::ConstraintRelation;
    use crate::io::mps::parse::{check_row_section_consistency, ColumnRetriever, into_lines, parse_program_name, parse_row_section};
    use crate::io::mps::parse::fixed::Fixed;
    use crate::io::mps::parse::free::Free;
    use crate::io::mps::Row;

    const MPS_LITERAL_STRING: &str = "
* Start of the file
NAME          TESTPROB

ROWS
* This is the cost row
 N  COST
 L  LIM1
 G  LIM2
 E  MYEQN
COLUMNS
    XONE      COST                 1   LIM1                 1
    XONE      LIM2                 1
    MARKER0   'MARKER'      'INTORG'
    YTWO      COST                 4   LIM1                 1
    YTWO      MYEQN               -1
    MARKER0   'MARKER'      'INTEND'
    ZTHREE    COST                 9   LIM2                 1
    ZTHREE    MYEQN                1
RHS
    RHS1      LIM1                 5   LIM2                10
    RHS1      MYEQN                7
BOUNDS
 UP BND1      XONE                 4
 LO BND1      YTWO                -1
 UP BND1      YTWO                 1
ENDATA";

    #[test]
    fn test_into_lines() {
        let result = into_lines(MPS_LITERAL_STRING).collect::<Vec<_>>();
        let expected = vec![
            (3, "NAME          TESTPROB"),
            (5, "ROWS"),
            (7, " N  COST"),
        ];
        assert_eq!(result[..expected.len()], expected);
        assert_eq!(result.last(), Some(&(27, "ENDATA")));
    }

    fn test_parse_program_name<'a, CR: ColumnRetriever<'a>>() {
        let mut lines = into_lines(MPS_LITERAL_STRING);
        let result = parse_program_name::<CR>(lines.next());

        assert_eq!(result.ok(), Some("TESTPROB".to_string()));
    }
    #[test]
    fn test_parse_program_name_free() {
        test_parse_program_name::<Free>()
    }
    #[test]
    fn test_parse_program_name_fixed() {
        test_parse_program_name::<Fixed>()
    }

    fn test_parse_row_section<'a, CR: ColumnRetriever<'a>>() {
        let mut lines = into_lines(MPS_LITERAL_STRING);
        let _name = lines.next();
        let _section = lines.next();
        let result = parse_row_section::<CR, _>(&mut lines);

        assert_eq!(result.ok(), Some((Some("COST".to_string()), vec![
            Row {
                name: "LIM1".to_string(),
                constraint_type: ConstraintRelation::Less,
            },
            Row {
                name: "LIM2".to_string(),
                constraint_type: ConstraintRelation::Greater,
            },
            Row {
                name: "MYEQN".to_string(),
                constraint_type: ConstraintRelation::Equal,
            },
        ])));
        assert!(lines.next().unwrap().1.starts_with("    XONE"));
    }
    #[test]
    fn test_parse_row_section_free() {
        test_parse_row_section::<Free>()
    }
    #[test]
    fn test_parse_row_section_fixed() {
        test_parse_row_section::<Fixed>()
    }

    #[test]
    fn test_check_row_section_consistency() {
        assert!(check_row_section_consistency(None, vec![]).is_err());
        assert!(check_row_section_consistency(Some("a".to_string()), vec![
            Row {
                name: "a".to_string(),
                constraint_type: ConstraintRelation::Equal
            },
            Row {
                name: "a".to_string(),
                constraint_type: ConstraintRelation::Equal
            },
        ]).is_err());
        assert!(check_row_section_consistency(Some("a".to_string()), vec![
            Row {
                name: "a".to_string(),
                constraint_type: ConstraintRelation::Equal
            },
            Row {
                name: "b".to_string(),
                constraint_type: ConstraintRelation::Equal
            },
        ]).is_err());
        assert!(check_row_section_consistency(Some("a".to_string()), vec![
            Row {
                name: "b".to_string(),
                constraint_type: ConstraintRelation::Equal
            },
            Row {
                name: "c".to_string(),
                constraint_type: ConstraintRelation::Equal
            },
        ]).is_ok());
    }
}
