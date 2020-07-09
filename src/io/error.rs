//! # Error reporting for reading of linear program files
//!
//! A collection of enums and structures describing any problems encountered during reading and
//! parsing.
use core::fmt::Display;
use std::convert::Into;
use std::error::Error;
use std::fmt;
use std::io;

/// An `ImportError` is created when an error was encountered during IO or parsing.
///
/// It is the highest error in the io error hierarchy.
#[derive(Debug)]
pub enum Import {
    /// The file extension of the provided file path is not known or supported.
    ///
    /// The contained `String` is a message for the end user.
    FileExtension(String),
    /// The file to read isn't found, or the reading of file couldn't start or was interrupted.
    IO(io::Error),
    /// Contents of the file could not be parsed into a linear program.
    ///
    /// # Note
    ///
    /// If the linear program is inconsistent, that will not be represented with this error. This
    /// variant should only be created for syntactically incorrect files.
    Parse(Parse),
    /// There is a logical inconsistency in the linear program described by a file.
    ///
    /// For example, a bound might be given for a variable which is not known.
    LinearProgram(Inconsistency),
}

impl From<Parse> for Import {
    fn from(error: Parse) -> Self {
        Self::Parse(error)
    }
}
impl From<Inconsistency> for Import {
    fn from(error: Inconsistency) -> Self {
        Self::LinearProgram(error)
    }
}

impl Display for Import {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Import::FileExtension(err) => err.fmt(f),
            Import::IO(error) => error.fmt(f),
            Import::Parse(error) => error.fmt(f),
            Import::LinearProgram(error) => error.fmt(f),
        }
    }
}

impl Error for Import {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            Import::FileExtension(_) => None,
            Import::IO(error) => error.source(),
            Import::Parse(error) => error.source(),
            Import::LinearProgram(_error) => None,
        }
    }
}

/// A `ParseError` represents all errors encountered during parsing.
///
/// It may recursively hold more ParseErrors to provide more detail. At the end of this chain, there
/// may be a file location containing a line number and line, at which the error was caused.
#[derive(Debug)]
pub struct Parse {
    description: String,
    source: Option<Source>,
}
pub type ParseResult<T> = Result<T, Parse>;

impl Display for Parse {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        writeln!(f, "ParseError: {}", self.description)
    }
}

impl Error for Parse {
    /// Find out what caused this error.
    ///
    /// # Return value
    ///
    /// Option<&Error> which may be a `ParseErrorCause`.
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        if let Some(Source::Other(ref error)) = self.source {
            Some(error.as_ref())
        } else { None }
    }

    /// Describe this error.
    ///
    /// # Return value
    ///
    /// The description as a `&str`.
    fn description(&self) -> &str {
        self.description.as_str()
    }
}

impl Parse {
    /// Create a new `ParseError` with only a description.
    ///
    /// # Arguments
    ///
    /// * `description`: What's wrong at the moment of creation.
    ///
    /// # Return value
    ///
    /// * A `ParseError` instance without a cause.
    pub fn new(description: impl Into<String>) -> Self {
        Self { description: description.into(), source: None, }
    }

    /// Create a new `ParseError` instance with a `FileLocation` as a cause.
    ///
    /// # Arguments
    ///
    /// * `description`: What's wrong at the moment of creation.
    /// * `file_location`: A reference to a line number and line that caused the error.
    ///
    /// # Return value
    ///
    /// A new `ParseError` instance with a `FileLocation` cause.
    pub fn with_location(
        description: impl Into<String>,
        file_location: FileLocation,
    ) -> Self {
        let (line_number, line) = file_location;
        Self {
            description: description.into(),
            source: Some(Source::FileLocation(line_number, line.to_string())),
        }
    }

    pub fn wrap(self, description: impl Into<String>) -> Self {
        Self {
            description: description.into(),
            source: Some(Source::Parse(Box::new(self))),
        }
    }

    /// Wrap a new `ParseError` around an existing one.
    ///
    /// # Arguments
    ///
    /// * `description`: What's wrong at the moment of creation.
    /// * `parse_error`: What caused this `ParseError`.
    ///
    /// # Return value
    ///
    /// A new `ParseError` instance with a `ParseError` cause.
    pub fn wrap_other(source: impl Error + 'static, description: impl Into<String>) -> Self {
        Self {
            description: description.into(),
            source: Some(Source::Other(Box::new(source))),
        }
    }

    /// Get all errors in the chain, leading up to this one.
    fn trace(&self) -> Vec<String> {
        let mut descriptions = if let Some(ref source) = self.source {
            match source {
                Source::FileLocation(line_number, line) => {
                    vec![format!("Caused at line {}: \t{}", line_number, line)]
                },
                Source::Parse(error) => error.trace(),
                Source::Other(error) => vec![error.to_string()],
            }
        } else {
            Vec::new()
        };

        descriptions.push(self.to_string());
        descriptions
    }
}

/// A `ParseErrorCause` can be used with a `ParseError` to describe its cause.
///
/// It can be either a file line number and line contents, or another `ParseError` with its own
/// description and optionally, a cause.
#[derive(Debug)]
enum Source {
    FileLocation(usize, String),
    Parse(Box<Parse>),
    Other(Box<dyn Error>),
}

/// A `FileLocation` references a line in the file by the line number of the file as originally
/// read from the disk. It contains a reference to the line itself.
pub(super) type FileLocation<'a> = (usize, &'a str);

/// A `LinearProgramError` is thrown when the linear program is inconsistently represented in the
/// file.
///
/// This `Error` is not returned when the linear program is unfeasible or unbounded. It is meant
/// only for descriptions of linear programs, and should not be used after the importing process.
#[derive(Debug, Eq, PartialEq)]
pub struct Inconsistency {
    description: String,
}

impl Inconsistency {
    /// Wrap a text in a `LinearProgramError`.
    ///
    /// # Arguments
    ///
    /// * `description`: A human-readable text mean for the end user.
    ///
    /// # Return value
    ///
    /// A `LinearProgramError`.
    pub fn new(description: impl Into<String>) -> Inconsistency {
        Inconsistency { description: description.into(), }
    }
}

impl Display for Inconsistency {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ProgramError: {}", self.to_string())
    }
}

impl Error for Inconsistency {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        Some(self)
    }
}
