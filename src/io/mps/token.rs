//! # Tokens used in MPS files

/// Indicates the start of a comment.
pub const COMMENT_INDICATOR: &str = "*";

/// Should be on the start of the first line of the non comment lines.
pub const NAME: &str = "NAME";

/// Indicates a line denoting a change in variable type.
///
/// This change is either from continuous to integer, or vice versa.
///
/// # Note
///
/// This token is expected only in the COLUMN section.
pub const COLUMN_SECTION_MARKER: &str = "'MARKER'";

/// Marks the start of the integer variables.
///
/// # Note
///
/// Expected only on a line with a `COLUMN_SECTION_MARKER`.
pub const START_OF_INTEGER: &str = "'INTORG'";

/// Indicates the end of the integer variables.
///
/// More data for continuous variables may follow.
///
/// # Note
///
/// Expected only on a line with a `COLUMN_SECTION_MARKER`.
pub const END_OF_INTEGER: &str = "'INTEND'";
