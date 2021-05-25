//! # Identity column
//!
//! The trait is useful for artificial tableaus: if the first phase is to be used, one needs to
//! represent the artificial columns, which are identity columns, somehow.
//!
//! Mostly useful for debugging, this module also contains a trivial identity column implementation
//! using a custom type for the constant `1`. As computing an identity column w.r.t. the current
//! basis is equivalent to just inverting the basis explicitly.
use std::iter;

use relp_num::One;

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
