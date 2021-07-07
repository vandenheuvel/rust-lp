//! # A linear program solver
//!
//! Linear programs are solved using the Simplex Method as described in the book Combinatorial
//! Optimization by Christos H. Papadimitriou and Kenneth Steiglitz.
#![warn(missing_docs)]

#![allow(incomplete_features)]
#![feature(is_sorted)]
#![feature(drain_filter)]
#![feature(result_cloned)]
#![feature(specialization)]
#![feature(min_type_alias_impl_trait)]
#![feature(associated_type_bounds)]
#![feature(trait_alias)]
#![feature(generic_associated_types)]
#![feature(const_generics)]
#![feature(map_first_last)]
#![feature(unchecked_math)]
#![feature(result_into_ok_or_err)]
#![feature(label_break_value)]

pub mod algorithm;
pub mod data;
pub mod io;

#[cfg(test)]
mod tests;
