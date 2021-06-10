//! # Sparse vector
//!
//! Wrapping a `Vec<(usize, _)>`, fixed size.
use std::{fmt, mem};
use std::borrow::Borrow;
use std::cmp::Ordering;
use std::collections::HashSet;
use std::fmt::{Debug, Display};
use std::iter::Sum;
use std::marker::PhantomData;
use std::ops::{Add, AddAssign, DivAssign, Mul, MulAssign, Neg};
use std::slice::Iter;

use num_traits::{One, Zero};
use relp_num::NonZero;

use crate::algorithm::utilities::remove_sparse_indices;
use crate::data::linear_algebra::SparseTuple;
use crate::data::linear_algebra::traits::{SparseComparator, SparseElement};
use crate::data::linear_algebra::vector::{DenseVector, Vector};

/// A sparse vector using a `Vec` with (row, value) combinations as back-end. Indices start at
/// `0`.
///
/// TODO(ENHANCEMENT): Consider making this backed by a `HashMap`.
#[derive(Eq, PartialEq, Clone, Debug)]
pub struct Sparse<F, C> {
    data: Vec<SparseTuple<F>>,
    len: usize,

    /// The level that comparison is done at: a single reference to the underlying data.
    phantom_comparison_type: PhantomData<C>,
}

impl<F, C> Sparse<F, C> {
    fn get_data_index(&self, i: usize) -> Result<usize, usize> {
        self.data.binary_search_by_key(&i, |&(index, _)| index)
    }

    fn set_zero(&mut self, i: usize) {
        if let Ok(index) = self.get_data_index(i) {
            self.data.remove(index);
        }
    }

    /// Increase the length of the vector by passing with zeros.
    pub fn extend(&mut self, extra_len: usize) {
        self.len += extra_len;
    }

    /// Convert the inner data structure into an iterator, consuming the struct.
    pub fn into_iter(self) -> impl Iterator<Item=SparseTuple<F>> {
        self.data.into_iter()
    }
}

impl<F, C> Vector<F> for Sparse<F, C>
where
    F: SparseElement<C>,
    C: SparseComparator,
{
    type Inner = SparseTuple<F>;

    /// Create a vector of length `len` from `data`.
    ///
    /// Requires that values close to zero are already filtered.
    fn new(data: Vec<Self::Inner>, len: usize) -> Self {
        debug_assert!(data.iter().all(|&(i, _)| i < len));
        debug_assert!(data.is_sorted_by_key(|&(i, _)| i));
        debug_assert!(data.iter().all(|(_, v)| v.borrow().is_not_zero()));
        debug_assert_ne!(len, 0);
        debug_assert!(data.len() <= len);

        Self {
            data,
            len,

            phantom_comparison_type: PhantomData,
        }
    }

    fn sparse_inner_product<'a, H, G: 'a, I: Iterator<Item=&'a SparseTuple<G>>>(&self, column: I) -> H
    where
        H: Zero + AddAssign<F> + Display + Debug,
        G: Display + Debug,
        for<'r> &'r F: Mul<&'r G, Output=F>,
    {
        let mut total = H::zero();

        let mut i = 0;
        for (index, value) in column {
            while i < self.data.len() && self.data[i].0 < *index {
                i += 1;
            }

            if i < self.data.len() && self.data[i].0 == *index {
                total += &self.data[i].1 * value;
                i += 1;
            }

            if i == self.len {
                break;
            }
        }

        total
    }

    /// Append a non-zero value.
    fn push_value(&mut self, value: F) {
        debug_assert!(value.borrow().is_not_zero());
        self.data.push((self.len, value));
        self.len += 1;
    }

    /// Set the value at index `i` to `value`.
    ///
    /// # Arguments
    ///
    /// * `i`: Index of the value. New tuple will be inserted, potentially causing many values to
    /// be shifted.
    /// * `value`: Value to be taken at index `i`. Should not be very close to zero to avoid
    /// memory usage and numerical error build-up.
    fn set(&mut self, i: usize, value: F) {
        debug_assert!(i < self.len);
        debug_assert!(value.borrow().is_not_zero());

        match self.get_data_index(i) {
            Ok(index) => self.data[index].1 = value,
            Err(index) => self.data.insert(index, (i, value)),
        }
    }

    fn get(&self, index: usize) -> Option<&F> {
        debug_assert!(index < self.len);

        self.get_data_index(index).ok().map(|i| &self.data[i].1)
    }

    /// Remove elements.
    ///
    /// # Arguments
    ///
    /// * `indices` is assumed sorted.
    fn remove_indices(&mut self, indices: &[usize]) {
        debug_assert!(indices.is_sorted());
        // All values are unique
        debug_assert!(indices.iter().collect::<HashSet<_>>().len() == indices.len());
        debug_assert!(indices.iter().all(|&i| i < self.len));
        debug_assert!(indices.len() < self.len);

        remove_sparse_indices(&mut self.data, indices);
        self.len -= indices.len();
    }

    fn iter_values(&self) -> Iter<Self::Inner> {
        self.data.iter()
    }

    /// The length of this vector.
    fn len(&self) -> usize {
        self.len
    }

    /// Whether this vector has zero size.
    fn is_empty(&self) -> bool {
        self.len == 0
    }

    /// The size of this vector in memory.
    fn size(&self) -> usize {
        self.data.len()
    }
}

impl<F, C> Sparse<F, C>
where
    F: SparseElement<C>,
    C: SparseComparator,
{
    /// Create a `SparseVector` representation of standard basis unit vector e_i.
    ///
    /// # Arguments
    ///
    /// * `i`: Only index where there should be a 1. Note that indexing starts at zero, and runs
    /// until (not through) `len`.
    /// * `len`: Size of the `SparseVector`.
    #[must_use]
    pub fn standard_basis_vector(i: usize, len: usize) -> Self
    where
        F: One + Clone,
    {
        debug_assert!(i < len);

        Self::new(vec![(i, F::one())], len)
    }

    /// Add the multiple of another row to this row.
    ///
    /// # Arguments
    ///
    /// * `multiple`: Non-zero constant that all elements of the `other` vector are multiplied with.
    /// * `other`: Vector to add a multiple of to this vector.
    ///
    /// # Return value
    ///
    /// A new `SparseVector`.
    ///
    /// # Note
    ///
    /// The implementation of this method doesn't look pretty, but it seems to be reasonably fast.
    /// If this method is too slow, it might be wise to consider the switching of the `SparseVector`
    /// storage backend from a `Vec` to a `HashMap`.
    pub fn add_multiple_of_row<H>(&mut self, multiple: &F, other: &Sparse<F, C>)
    where
        H: Zero + Add<F, Output=H>,
        F: Add<F, Output=H> + From<H>,
        for<'r> &'r F: Mul<&'r F, Output=F>,
    {
        debug_assert_eq!(other.len(), self.len());
        debug_assert!(multiple.borrow().is_not_zero());

        let mut new_tuples = Vec::new();

        let mut j = 0;  // data index
        let old_data = mem::replace(&mut self.data, Vec::with_capacity(0));
        for (i, value) in old_data {
            while j < other.data.len() && other.data[j].0 < i {
                let new_value = multiple * &other.data[j].1;
                new_tuples.push((other.data[j].0, new_value.into()));
                j += 1;
            }

            if j < other.data.len() && i == other.data[j].0 {
                let new_value = value + multiple * &other.data[j].1;
                if !new_value.is_zero() {
                    new_tuples.push((i, new_value.into()));
                }
                j += 1;
            } else {
                new_tuples.push((i, value));
            }
        }
        for (j, value) in &other.data[j..] {
            new_tuples.push((*j, multiple * value));
        }

        self.data = new_tuples;
    }

    /// Set the value at index `i` to `value`.
    ///
    /// # Arguments
    ///
    /// * `i`: Index of the value. New tuple will be inserted, potentially causing many values to
    /// be shifted.
    /// * `value`: Value to be taken at index `i`.
    pub fn shift_value<G>(&mut self, i: usize, value: G)
    where
        F: PartialEq<G> + AddAssign<G> + From<G>,
        for<'r> &'r G: Neg<Output=G>,
        G: NonZero,
    {
        debug_assert!(i < self.len);

        if value.is_not_zero() {
            match self.get_data_index(i) {
                Ok(index) => {
                    if self.data[index].1 == -&value {
                        self.set_zero(i);
                    } else {
                        self.data[index].1 += value;
                    }
                },
                Err(index) => self.data.insert(index, (i, From::from(value))),
            }
        }
    }

    /// Multiply each element of the vector by a value.
    pub fn element_wise_multiply(&mut self, value: &F)
    where
        for<'r> F: NonZero + MulAssign<&'r F>,
    {
        debug_assert!(value.is_not_zero());

        for (_, v) in &mut self.data {
            *v *= value;
        }
    }

    /// Divide each element of the vector by a value.
    pub fn element_wise_divide(&mut self, value: &F)
    where
       for<'r> F: NonZero + DivAssign<&'r F>,
    {
        debug_assert!(value.is_not_zero());

        for (_, v) in &mut self.data {
            *v /= value;
        }
    }
}

impl<F, C> Sparse<F, C>
where
    F: SparseElement<C>,
    C: SparseComparator,
{
    /// Calculate the inner product between two vectors.
    ///
    /// # Arguments
    ///
    /// * `other`: Vector to calculate inner product with.
    ///
    /// # Return value
    ///
    /// The inner product.
    #[must_use]
    pub fn inner_product_with_dense<'a, F2, O: 'a>(&'a self, other: &'a DenseVector<F2>) -> O
    where
        F: Borrow<O>,
        F2: Borrow<O> + PartialEq + Display + Debug,
        O: Sum,
        &'a O: Mul<&'a O, Output=O>,
    {
        debug_assert_eq!(other.len(), self.len());

        self.data.iter().map(|(i, value)| other[*i].borrow() * value.borrow()).sum()
    }

    /// Calculate the inner product between two vectors.
    ///
    /// # Arguments
    ///
    /// * `other`: Vector to calculate inner product with.
    ///
    /// # Return value
    ///
    /// The inner product.
    #[must_use]
    pub fn inner_product<'a, O, F2>(&'a self, other: &'a Sparse<F2, C>) -> O
    where
        O: Zero + AddAssign<C>,
        F2: SparseElement<C>,
        // We choose to have multiplication output at the C level, because it would also be nonzero
        // if both F and F2 values are not zero.
        &'a C: Mul<&'a C, Output=C>,
    {
        debug_assert_eq!(other.len(), self.len());
        debug_assert!(self.data.iter().all(|(_, v)| v.borrow().is_not_zero()));
        debug_assert!(other.data.iter().all(|(_, v)| v.borrow().is_not_zero()));

        let mut self_lowest = 0;
        let mut other_lowest = 0;

        let mut total = O::zero();
        while self_lowest < self.data.len() && other_lowest < other.data.len() {
            let self_sought = self.data[self_lowest].0;
            let other_sought = other.data[other_lowest].0;
            match self_sought.cmp(&other_sought) {
                Ordering::Less => {
                    match self.data[self_lowest..].binary_search_by_key(&other_sought, |&(i, _)| i) {
                        Err(diff) => {
                            self_lowest += diff;
                            other_lowest += 1;
                        },
                        Ok(diff) => {
                            total += self.data[self_lowest + diff].1.borrow() * other.data[other_lowest].1.borrow();
                            self_lowest += diff + 1;
                            other_lowest += 1;
                        },
                    }
                },
                Ordering::Greater => {
                    match other.data[other_lowest..].binary_search_by_key(&self_sought, |&(i, _)| i) {
                        Err(diff) => {
                            self_lowest += 1;
                            other_lowest += diff;
                        },
                        Ok(diff) => {
                            total += self.data[self_lowest].1.borrow() * other.data[other_lowest + diff].1.borrow();
                            self_lowest += 1;
                            other_lowest += diff + 1;
                        },
                    }
                },
                Ordering::Equal => {
                    total += self.data[self_lowest].1.borrow() * other.data[other_lowest].1.borrow();
                    self_lowest += 1;
                    other_lowest += 1;
                },
            }
        }

        total
    }
}

impl<F: SparseElement<C>, C: SparseComparator> Display for Sparse<F, C> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for (index, value) in &self.data {
            writeln!(f, "({} {}), ", index, value)?;
        }
        writeln!(f)
    }
}
