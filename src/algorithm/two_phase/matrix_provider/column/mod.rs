//! # Columns
//!
//! Items yielded by a matrix provider. Used to iterate over when computing what they look like
//! w.r.t. the latest basis.
use std::fmt::{Debug, Display};

use crate::data::linear_algebra::SparseTuple;

pub mod identity;

/// Columns represent part of a (virtual) data matrix.
///
/// This column is sparse.
///
/// A column might be expensive to compute. It can store computed values and once GATs work
/// (better), also references to items stored in the matrix provider that yields instances of this
/// trait.
///
/// It can't necessarily be iterated over directly. That needs to happen many times, for example
/// when computing many inner products with a data matrix. The trait has an associated type to be
/// used for iteration, that should be cheaply cloneable and probably not store any values itself.
/// Rather, it should describe how this column should be iterated over.
///
/// TODO(ARCHITECTURE): Many basis inverse maintenance algorithms require reallocation of the
///  column. Is this more complex set-up worth it?
// TODO(ARCHITECTURE): Once GATs work, consider giving this trait a lifetime parameter.
pub trait Column: Clone {
    /// Input data type.
    ///
    /// Items of this type get read and used in additions and multiplications often.
    // TODO(ENHANCEMENT): Don't work with a field type directly, but an `Into<F>` type to separate.
    type F: ColumnNumber;

    /// Type of struct to iterate over this column.
    ///
    /// It should be somewhat cheaply cloneable and as such not be too large.
    type Iter<'a, G: 'a>: Iterator<Item = &'a SparseTuple<G>> + Clone;

    /// Derive the iterator object.
    ///
    /// Because this column might need to be iterated over many times, it doesn't consume the
    /// column but instead produces a struct that might keep references to this column.
    fn iter(&self) -> Self::Iter<'_, Self::F>;

    /// Format an index of the column.
    ///
    /// Note that this index might not be explicitly stored due to the column being sparse.
    fn index_to_string(&self, i: usize) -> String;
}

/// Basic operations that should be possible with the type of the values of the column.
///
/// This trait was introduced to avoid some verbose trait bounds throughout the codebase.
pub trait ColumnNumber =
    Eq +
    PartialEq +

    Display +
    Clone +
    Debug +
;

/// Column that can be iterated over in-order.
///
/// This trait is simply a marker trait to be used in specialization.
///
/// TODO(ENHANCEMENT): At the time of writing, it is not possible to specialize the generic
///  arguments of trait methods. That is why this trait and the standard `Column` trait are
///  currently both needed.
///
// TODO(ARCHITECTURE): Once GATs work, consider giving this trait a lifetime parameter.
pub trait OrderedColumn: Column {
}

/// Wrapping a sparse vector of tuples.
#[derive(Clone)]
#[allow(missing_docs)]
pub struct SparseColumn<F> {
    pub inner: Vec<SparseTuple<F>>
}
impl<F: ColumnNumber> Column for SparseColumn<F> {
    type F = F;
    type Iter<'a, G: 'a> = impl Iterator<Item = &'a SparseTuple<G>> + Clone;

    fn iter(&self) -> Self::Iter<'_, Self::F> {
        self.inner.iter()
    }

    fn index_to_string(&self, i: usize) -> String {
        match self.inner.iter().find(|&&(ii, _)| ii == i) {
            None => "0".to_string(),
            Some((_, v)) => v.to_string(),
        }
    }
}
impl<F: 'static + ColumnNumber> OrderedColumn for SparseColumn<F> {
}
