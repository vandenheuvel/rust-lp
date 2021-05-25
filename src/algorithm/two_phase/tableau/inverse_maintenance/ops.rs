//! # Number operations
//!
//! During the interaction with each of a `MatrixProvider`'s associated types it should be possible
//! to do certain operations.
use std::fmt::{Debug, Display};
use std::iter::Sum;
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

use num_traits::{One, Zero};
use relp_num;
use relp_num::NonZero;

use crate::data::linear_algebra::traits::SparseElement;

/// Operations done by the number type within the inverse maintenance algorithm.
pub trait Field =
    Zero +
    NonZero +
    One +

    Neg<Output=Self> +

    for<'r> Add<&'r Self, Output=Self> +
    AddAssign +
    for<'r> AddAssign<&'r Self> +

    Sub<Self, Output=Self> +
    for<'r> Sub<&'r Self, Output=Self> +
    SubAssign +
    for<'r> SubAssign<&'r Self> +

    for<'r> Mul<&'r Self, Output=Self> +
    MulAssign +
    for<'r> MulAssign<&'r Self> +

    Div<Self, Output=Self> +
    for<'r> Div<&'r Self, Output=Self> +
    DivAssign<Self> +
    for<'r> DivAssign<&'r Self> +

    Sum +

    Column<relp_num::One> +

    Eq +
    PartialEq +
    Ord +
    PartialOrd +

    SparseElement<Self> +

    Clone +
    Debug +
    Display +
;

// TODO(ARCHITECTURE): Once HRTB are propagated like normal associated type trait bounds, remove
//  this trait by integrating the requirements into `InverseMaintenance::F`'s trait bounds.
#[allow(clippy::type_repetition_in_bounds)]
pub trait FieldHR =
where
    for<'r> &'r Self: Neg<Output=Self>,
    for<'r> &'r Self: Mul<&'r Self, Output=Self>,
    for<'r> &'r Self: Div<&'r Self, Output=Self>,
;

/// Operations with the values in the columns.
pub trait Column<Rhs> =
    for<'r> AddAssign<&'r Rhs> +
    for<'r> Add<&'r Rhs, Output=Self> +

    From<Rhs> +
    for<'r> From<&'r Rhs> +
where
    for<'r> &'r Self: Mul<&'r Rhs, Output=Self>,
;

/// Operations with the cost type.
pub trait Cost<Rhs> =
    Add<Rhs, Output=Self> +
    Mul<Rhs, Output=Self> +
where
    for<'r> &'r Self: Mul<Rhs, Output=Self>,
;

/// Operations with the right-hand side.
pub trait Rhs<Rhs> =
    for<'r> AddAssign<&'r Rhs> +
    for<'r> Add<&'r Rhs, Output=Self> +

    From<Rhs> +
;
