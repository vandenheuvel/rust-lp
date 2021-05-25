//! # Traits for elements in sparse datastructures
//! 
//! Once a sparse data structure contains references to values, it is not obvious what value should
//! be returned for a zero value, that is not stored. It is also not clear, how one should compare
//! the elements contained in the sparse structure with the zero element (mostly for debug 
//! purposes).
//! 
//! One idea, implemented here, is to have three types related to a sparse data structure:
//! 
//! * The first is the type that is stored many times in the data structure
//! * The second is a type that can be zero (and is not a reference), ideally is small (Copy) and 
//! not stored behind a reference (like `RationalBig`).
//! * The third is the type that both can be dereferenced to. This is used to create a row-major
//! copy of the constraint matrix using references, rather than the actual values.
use std::borrow::Borrow;
use std::fmt::{Debug, Display};

use relp_num::NonZero;

/// Element of a `Vector` of `Matrix` type.
///
/// This is an alias for a traits that are needed to derive a few practical traits for the
/// aforementioned types.
pub trait Element =
    PartialEq +
    Clone +
    Display +
    Debug +
;

/// Element of a sparse data structure.
///
/// Needs to be borrowable as a type that can be used for comparison with the zero type, as well as
/// reference types.
pub trait SparseElement<Comparator> =
    Borrow<Comparator> +
    Element +
;

/// Element to do comparisons between vectors of different "reference levels".
///
/// We might have an inner product between a `SparseVector<F>` and `SparseVector<&F>`. Then this
/// comparator type would be `F`, such that the values can be compared as `&F`'s.
pub trait SparseComparator =
    PartialEq +
    Element +
    NonZero +
;
