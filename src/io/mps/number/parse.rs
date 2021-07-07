//! # Number parsing
//!
//! Reading numbers from strings.
use std::convert::TryFrom;

use relp_num::Rational64;
use relp_num::Sign;

use crate::io::error::Parse as ParseError;
use crate::io::error::ParseResult;

/// Parsing a number read from an MPS file.
pub trait Parse: Sized {
    /// Read a string representation from a decimal (abc.xyz) number.
    ///
    /// # Errors
    ///
    /// When the number could not be parsed, an error.
    fn parse(text: &str) -> ParseResult<Self>;
}

impl Parse for f64 {
    fn parse(text: &str) -> ParseResult<Self> {
        let value: Self = text.parse()
            .map_err(|error| ParseError::wrap_other(
                error,
                format!("Failed to parse value text \"{}\" into f64", text),
            ))?;

        let minimum_magnitude = 1e-10;
        let maximum_magnitude = 1e10;
        if value.abs() < minimum_magnitude {
            Err(ParseError::new(format!(
                "Parsed value was {}, but can't be smaller than {}.", value, minimum_magnitude,
            )))
        } else if value.abs() > maximum_magnitude {
            Err(ParseError::new(format!(
                "Parsed value was {}, but can't be larger than {}.", value, maximum_magnitude,
            )))
        } else {
            Ok(value)
        }
    }
}

impl Parse for Rational64 {
    fn parse(text: &str) -> ParseResult<Self> {
        let raw = Raw::try_from(text)?;
        let value = raw.into();

        Ok(value)
    }
}

impl From<Raw> for Rational64 {
    fn from(value: Raw) -> Self {
        let Raw { sign, integer, decimal_steps_from_right } = value;
        debug_assert!(sign != Sign::Zero || integer == 0);

        let numerator = integer;
        let denominator = 10_u64.pow(decimal_steps_from_right);

        Self::new_signed(sign, numerator, denominator)
    }
}

/// Intermediate form of read number.
///
/// TODO(PERFORMANCE): Reduce the size of this struct to 64 bits.
#[derive(Eq, PartialEq, Clone, Debug)]
pub(crate) struct Raw {
    sign: Sign, // Need 2 bits
    integer: u64, // We need about 40 bits to represent
    decimal_steps_from_right: u32, // Need about 5 bits
}

impl TryFrom<&str> for Raw {
    type Error = ParseError;

    fn try_from(text: &str) -> Result<Self, Self::Error> {
        let (has_minus, text) = match &text[0..1] {
            "-" => (true, &text[1..]),
            _ => (false, text),
        };

        let parse = |text: &str, number_part| {
            if text.is_empty() {
                Ok(0)
            } else {
                text.parse()
                    .map_err(|error| ParseError::wrap_other(
                        error,
                        format!("Failed to parse {} \"{}\" as u64.", number_part, text),
                    ))
            }
        };

        let (integer, decimal_steps_from_right): (u64, _) = match text.find('.') {
            None => (parse(text, "entire value")?, 0),
            Some(index) => {
                let from_right = u32::try_from(text.len() - index - 1)
                    .expect("It wouldn't make sense for the number text to be so large.");

                let integer_part = parse(&text[0..index], "integer part")?;
                let mantissa_part = parse(&text[(index + 1)..], "mantissa part")?;

                (integer_part * 10_u64.pow(from_right) + mantissa_part, from_right)
            },
        };
        
        let sign = match (has_minus, integer) {
            (_, 0) => Sign::Zero,
            (false, _) => Sign::Positive,
            (true, _) => Sign::Negative,
        };

        Ok(Self { sign, integer, decimal_steps_from_right, })
    }
}

#[cfg(test)]
mod test {
    use std::convert::TryFrom;

    use relp_num::{R64, Rational64};
    use relp_num::Sign;

    use crate::io::mps::number::parse::Raw;

    #[test]
    fn parse() {
        // dot location, single digit
        assert_eq!(Raw::try_from("1").unwrap(), Raw {
            sign: Sign::Positive,
            integer: 1,
            decimal_steps_from_right: 0,
        });
        assert_eq!(Raw::try_from("2.").unwrap(), Raw {
            sign: Sign::Positive,
            integer: 2,
            decimal_steps_from_right: 0,
        });
        assert_eq!(Raw::try_from(".3").unwrap(), Raw {
            sign: Sign::Positive,
            integer: 3,
            decimal_steps_from_right: 1,
        });

        // sign
        assert_eq!(Raw::try_from("-1").unwrap(), Raw {
            sign: Sign::Negative,
            integer: 1,
            decimal_steps_from_right: 0,
        });
        assert_eq!(Raw::try_from("-2.").unwrap(), Raw {
            sign: Sign::Negative,
            integer: 2,
            decimal_steps_from_right: 0,
        });
        assert_eq!(Raw::try_from("-.3").unwrap(), Raw {
            sign: Sign::Negative,
            integer: 3,
            decimal_steps_from_right: 1,
        });

        // dot location, multiple digits
        assert_eq!(Raw::try_from("16456").unwrap(), Raw {
            sign: Sign::Positive,
            integer: 16456,
            decimal_steps_from_right: 0,
        });
        assert_eq!(Raw::try_from("64896848.").unwrap(), Raw {
            sign: Sign::Positive,
            integer: 64_896_848,
            decimal_steps_from_right: 0,
        });
        assert_eq!(Raw::try_from(".984654684").unwrap(), Raw {
            sign: Sign::Positive,
            integer: 984_654_684,
            decimal_steps_from_right: 9,
        });
        assert_eq!(Raw::try_from("-.95").unwrap(), Raw {
            sign: Sign::Negative,
            integer: 95,
            decimal_steps_from_right: 2,
        });

        // multiple digits, both sides
        assert_eq!(Raw::try_from("15465.2").unwrap(), Raw {
            sign: Sign::Positive,
            integer: 154652,
            decimal_steps_from_right: 1,
        });
        assert_eq!(Raw::try_from("1234.56789").unwrap(), Raw {
            sign: Sign::Positive,
            integer: 123_456_789,
            decimal_steps_from_right: 5,
        });
        assert_eq!(Raw::try_from("1.24654").unwrap(), Raw {
            sign: Sign::Positive,
            integer: 124_654,
            decimal_steps_from_right: 5,
        });

        // Zero
        assert_eq!(Raw::try_from("0").unwrap(), Raw {
            sign: Sign::Zero,
            integer: 0,
            decimal_steps_from_right: 0,
        });
        assert_eq!(Raw::try_from("0.").unwrap(), Raw {
            sign: Sign::Zero,
            integer: 0,
            decimal_steps_from_right: 0,
        });
        assert_eq!(Raw::try_from(".0").unwrap(), Raw {
            sign: Sign::Zero,
            integer: 0,
            decimal_steps_from_right: 1,
        });
    }

    #[test]
    fn convert() {
        // dot location, single digit
        assert_eq!(Rational64::from(Raw {
            sign: Sign::Positive,
            integer: 1,
            decimal_steps_from_right: 0,
        }), R64!(1));
        assert_eq!(Rational64::from(Raw {
            sign: Sign::Positive,
            integer: 2,
            decimal_steps_from_right: 0,
        }), R64!(2));
        assert_eq!(Rational64::from(Raw {
            sign: Sign::Positive,
            integer: 3,
            decimal_steps_from_right: 1,
        }), R64!(3, 10));

        // sign
        assert_eq!(Rational64::from(Raw {
            sign: Sign::Negative,
            integer: 1,
            decimal_steps_from_right: 0,
        }), R64!(-1, 1));
        assert_eq!(Rational64::from(Raw {
            sign: Sign::Negative,
            integer: 2,
            decimal_steps_from_right: 0,
        }), R64!(-2, 1));
        assert_eq!(Rational64::from(Raw {
            sign: Sign::Negative,
            integer: 3,
            decimal_steps_from_right: 1,
        }), R64!(-3, 10));

        // dot location, multiple digits
        assert_eq!(Rational64::from(Raw {
            sign: Sign::Positive,
            integer: 16456,
            decimal_steps_from_right: 0,
        }), R64!(16456, 1));
        assert_eq!(Rational64::from(Raw {
            sign: Sign::Positive,
            integer: 64_896_848,
            decimal_steps_from_right: 0,
        }), R64!(64_896_848, 1));
        assert_eq!(Rational64::from(Raw {
            sign: Sign::Positive,
            integer: 984_654_684,
            decimal_steps_from_right: 9,
        }), R64!(246_163_671, 250_000_000));

        // multiple digits, both sides
        assert_eq!(Rational64::from(Raw {
            sign: Sign::Positive,
            integer: 12,
            decimal_steps_from_right: 1,
        }), R64!(12, 10));
        assert_eq!(Rational64::from(Raw {
            sign: Sign::Positive,
            integer: 123_456_789,
            decimal_steps_from_right: 5,
        }), R64!(123_456_789, 100_000));
        assert_eq!(Rational64::from(Raw {
            sign: Sign::Positive,
            integer: 124_654,
            decimal_steps_from_right: 5,
        }), R64!(62_327, 50_000));

        // Zero
        assert_eq!(Rational64::from(Raw {
            sign: Sign::Zero,
            integer: 0,
            decimal_steps_from_right: 0,
        }), R64!(0, 1));
        assert_eq!(Rational64::from(Raw {
            sign: Sign::Zero,
            integer: 0,
            decimal_steps_from_right: 1,
        }), R64!(0, 1));
    }
}
