//! # A linear program solver
//!
//! Linear programs are solved using the Simplex Method as described in the book Combinatorial
//! Optimization by Christos H. Papadimitriou and Kenneth Steiglitz.
#![warn(missing_docs)]

#![feature(is_sorted)]
#![feature(or_patterns)]
#![feature(drain_filter)]
#![feature(result_cloned)]
#![feature(specialization)]
#![feature(type_alias_impl_trait)]
#![feature(associated_type_bounds)]
#![feature(trait_alias)]
#![feature(generic_associated_types)]
#![feature(const_generics)]

extern crate daggy;
extern crate num_traits;
#[macro_use]
extern crate smallvec;

pub mod algorithm;
pub mod data;
pub mod io;

#[cfg(test)]
mod tests;
