//! # Identity column
//!
//! The trait is useful for artificial tableaus: if the first phase is to be used, one needs to
//! represent the artificial columns, which are identity columns, somehow.
//!
//! Mostly useful for debugging, this module also contains a trivial identity column implementation
//! using a custom type for the constant `1`. As computing an identity column w.r.t. the current
//! basis is equivalent to just inverting the basis explicitly.
use std::{fmt, iter};
use std::fmt::Display;
use std::ops::{Add, AddAssign, Div};
use std::ops::Mul;

use crate::algorithm::two_phase::matrix_provider::column::{Column, OrderedColumn};

/// Identity columns are needed for artificial matrices.
///
/// When a matrix provider is to be used in the first phase, it should be possible to represent
/// identity columns.
pub trait IdentityColumn: Column {
    /// Create an identity column, placing a "1" at a certain index and "0"'s otherwise.
    ///
    /// # Arguments
    ///
    /// * `i`: Index at which the "1" should be placed.
    /// * `len`: Length of the column. Might not be used in an actual implementation.
    fn identity(i: usize, len: usize) -> Self;
}

/// A simple identity column useful for debugging, mostly.
#[derive(Clone)]
pub struct IdentityColumnStruct(pub (usize, One));

impl Column for IdentityColumnStruct {
    type F = One;
    type Iter<'a> = std::iter::Once<&'a (usize, Self::F)>;

    fn iter(&self) -> Self::Iter<'_> {
        iter::once(&self.0)
    }

    fn index_to_string(&self, i: usize) -> String {
        if i == self.0.0 {
            "1"
        } else {
            "0"
        }.to_string()
    }
}

impl OrderedColumn for IdentityColumnStruct {
}

macro_rules! define_ops {
    ($primitive:ident) => {
        impl From<One> for $primitive {
            fn from(_: One) -> Self {
                1
            }
        }

        impl From<&One> for $primitive {
            fn from(_: &One) -> Self {
                1
            }
        }

        impl Add<One> for $primitive {
            type Output = Self;

            fn add(self, _: One) -> Self::Output {
                self + 1
            }
        }

        impl Add<&One> for $primitive {
            type Output = Self;

            fn add(self, _: &One) -> Self::Output {
                self + 1
            }
        }

        impl AddAssign<&One> for $primitive {
            fn add_assign(&mut self, _: &One) {
                *self += 1;
            }
        }

        impl Mul<&One> for $primitive {
            type Output = Self;

            fn mul(self, _: &One) -> Self::Output {
                self
            }
        }

        impl Mul<&One> for &$primitive {
            type Output = $primitive;

            fn mul(self, _: &One) -> Self::Output {
                self.clone()
            }
        }

        impl Div<&One> for $primitive {
            type Output = Self;

            fn div(self, _: &One) -> Self::Output {
                self
            }
        }
    }
}

define_ops!(i32);
define_ops!(i64);
define_ops!(i128);
define_ops!(u32);
define_ops!(u64);
define_ops!(u128);
