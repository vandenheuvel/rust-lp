use std::fmt::{Debug, Display};

use num_traits::{Zero, One};
use std::ops::{Add, Div, Neg, Mul, AddAssign, DivAssign};
use crate::data::number_types::traits::Signed;

/// Can be zero
///
/// Used for -z, -pi
pub trait Field =
    Zero +

    Neg<Output = Self> +
    NegAssign +

    Eq +
    PartialEq +

    Clone +
    Display +
    Debug +
;

pub trait FieldNZ<NZ> =
    Add<NZ, Output = Self> +
    AddAssign<NZ> +

    Sized +

// where
//     for<'r> &'r Self: FieldRefNZ<NZ, Self>,
;

// pub trait FieldRefNZ<NZ, Output> =
//     Mul<NZ, Output = Output>
// ;

/// Sparse basis inverse (never zero because of sparse storage)
///
/// Used for B^{-1}
pub trait NonZero =
    One +

    PartialOrd +
    Eq +
    PartialEq +

    Clone +
    Display +
    Debug +

where
    for<'r> &'r Self: NonZeroRef,
;

pub trait NonZeroRef =
    Signed +
;

pub trait NonZeroNZ<NZ> =
    Mul<NZ, Output = Self> +

    Sized +

where
    for<'r, 's> &'r Self: NonZeroRefNZ<&'s NZ, Self>,
;

pub trait NonZeroRefNZ<NZ, OutputType> =
    Mul<NZ, Output = OutputType> +
;


/// Constraint vector (always non negative because two-phase) (can be zero)
///
/// Used for b
pub trait NonNegative =
    Zero +

    Eq +
    PartialEq +

    Clone +
    Display +
    Debug +
;
pub trait NonNegativeNZ<NZ> =
    DivAssign<NZ> + // TODO: Problematic as sign might not be positive (normalization of pivot row)
    // Mul<>
;
pub trait NonNegativeRef<NZ> =
    for<'r> Div<&'r NZ, Output = NZ> +
    for<'r> Mul<&'r NZ, Output = NZ> +
;
pub trait NonNegativeRefWithOutput<NZ, Output> =
    for<'r> Div<&'r NZ, Output = Output> +
    for<'r> Mul<&'r NZ, Output = Output> +
;

pub trait NegAssign {
    fn neg_assign(&mut self);
}
