use crate::data::number_types::rational::small::Rational64;

/// Shorthand for creating a rational number in tests.
#[macro_export]
macro_rules! R64 {
    ($value:expr) => {
        Rational64::from_f64($value as f64).unwrap()
    };
    ($numer:expr, $denom:expr) => {
        Rational64::new($numer, $denom)
    };
}

/// Shorthand for creating a rational number in tests.
#[macro_export]
macro_rules! BR {
    ($value:expr) => {
        Rational64::from_f64($value as f64).unwrap()
    };
    ($numer:expr, $denom:expr) => {
        Rational64::new($numer, $denom)
    };
}
