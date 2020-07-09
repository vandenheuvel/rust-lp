use std::convert::TryFrom;

use crate::io::error::Parse as ParseError;
use crate::io::error::ParseResult;
use crate::data::number_types::rational::small::Rational64;

pub trait Parse: Sized {
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

        let unsigned_numerator = integer;
        let denominator = 10_u64.pow(decimal_steps_from_right as u32);

        let signed_numerator = match sign {
            Sign::Positive => unsigned_numerator,
            Sign::Negative => -unsigned_numerator,
        };

        Self::new(signed_numerator, denominator)
    }
}

pub(crate) struct Raw {
    sign: Sign,
    integer: i64, // We need 40 bits to represent TODO: Nonzero, size efficiency, non negative
    decimal_steps_from_right: u8,
}

pub(crate) enum Sign {
    Positive,
    Negative,
}

impl TryFrom<&str> for Raw {
    type Error = ParseError;

    fn try_from(text: &str) -> Result<Self, Self::Error> {
        let (sign, text) = match &text[0..1] {
            "-" => (Sign::Negative, &text[1..]),
            _ => (Sign::Positive, text),
        };

        let parse = |text: &str, number_part| text.parse()
            .map_err(|error| ParseError::wrap_other(
                error,
                format!("Failed to parse {} \"{}\" as u64.", number_part, text),
            ));

        let (integer, decimal_steps_from_right) = match text.find('.') {
            None => (parse(text, "entire value")?, 0),
            Some(index) => {
                let from_right = (text.len() - index) as u8;

                let integer_part = parse(&text[0..index], "integer part")?;
                let mantissa_part = parse(&text[(index + 1)..], "mantissa part")?;

                (integer_part * 10_i64.pow(from_right as u32) + mantissa_part, from_right)
            },
        };

        Ok(Self { sign, integer, decimal_steps_from_right, })
    }
}
