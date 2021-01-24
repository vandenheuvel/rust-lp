//! # Vector types for linear programs
//!
//! Sparse and dense vectors. These were written by hand, because a certain specific set of
//! operations needs to be done quickly with these types.
use std::fmt::{Debug, Display};
use std::ops::{AddAssign, Mul};
use std::slice::Iter;

use num::Zero;

pub use dense::Dense as DenseVector;
pub use sparse::Sparse as SparseVector;

use crate::data::linear_algebra::SparseTuple;
use crate::data::number_types::nonzero::Nonzero;

mod dense;
mod sparse;

/// Defines basic ways to create or change a vector, regardless of back-end.
pub trait Vector<F>: PartialEq + Display + Debug {
    /// Items stored internally.
    type Inner;

    /// Create a new instance.
    ///
    /// # Arguments
    ///
    /// * `data`: Internal data values. Will not be changed and directly used for creation.
    /// * `len`: Length of the vector represented (and not necessarily of the internal data struct
    /// ure).
    ///
    /// # Return value
    ///
    /// Input data wrapped inside a vector.
    fn new(data: Vec<Self::Inner>, len: usize) -> Self;
    /// Compute the inner product with a column vector from a matrix.
    fn sparse_inner_product<'a, H, G: 'a, V: Iterator<Item=&'a SparseTuple<G>>>(&self, column: V) -> H
    where
        H: Zero + AddAssign<F>,
        for<'r> &'r F: Mul<&'r G, Output=F>,
    ;
    /// Make a vector longer by one, by adding an extra value at the end of this vector.
    fn push_value(&mut self, value: F) where F: Nonzero;
    /// Set the value at an index.
    ///
    /// Depending on internal representation, this can be an expensive operation (for `SparseVector`
    /// 's, the cost depends on the (lack of) sparsity.
    fn set(&mut self, index: usize, value: F) where F: Nonzero;
    /// Retrieve the value at an index.
    ///
    /// # Returns
    ///
    /// `None` if the representation is `Sparse` and the value at the index is zero.
    fn get(&self, index: usize) -> Option<&F>;
    /// Remove the items at the specified indices.
    fn remove_indices(&mut self, indices: &[usize]);
    /// Iterate over the internal values.
    fn iter_values(&self) -> Iter<Self::Inner>;
    /// Number of items represented by the vector.
    fn len(&self) -> usize;
    /// Whether the vector is empty.
    fn is_empty(&self) -> bool;
    /// Get the size of the internal data structure (and not of the represented vector).
    fn size(&self) -> usize;
}

#[cfg(test)]
#[allow(dead_code)]
pub mod test {
    //! Contains also some helper methods to be used in other test modules.

    use std::f64::EPSILON;

    use num::{FromPrimitive, NumCast, ToPrimitive, Zero};

    use crate::data::linear_algebra::traits::{SparseComparator, SparseElement};
    use crate::data::linear_algebra::vector::{DenseVector, SparseVector, Vector};
    use crate::data::number_types::traits::{Field, FieldRef};
    use crate::F;

    pub trait TestVector<F>: Vector<F> {
        fn from_test_data<T: ToPrimitive + Zero>(data: Vec<T>) -> Self;
    }
    impl<F: Field + FromPrimitive> TestVector<F> for DenseVector<F> {
        /// Create a `DenseVector` from the provided data.
        fn from_test_data<T: ToPrimitive + Zero>(data: Vec<T>) -> Self {
            Self {
                data: data.into_iter()
                    .map(|v| F::from_f64(v.to_f64().unwrap()).unwrap())
                    .collect(),
            }
        }
    }
    impl<F: SparseElement<C>, C: SparseComparator> TestVector<F> for SparseVector<F, C>
    where
        F: Field + FromPrimitive,
    {
        /// Create a `SparseVector` from the provided data.
        fn from_test_data<T: ToPrimitive + Zero>(data: Vec<T>) -> Self {
            debug_assert_ne!(data.len(), 0);

            let size = data.len();

            Self::new(
                data.into_iter()
                    .enumerate()
                    .filter(|(_, v)| !v.is_zero())
                    .map(|(i, v)| (i, F::from_f64(v.to_f64().unwrap()).unwrap()))
                    .collect(),
                size,
            )
        }
    }

    impl<F: SparseElement<C>, C: SparseComparator> SparseVector<F, C>
    where
        F: FromPrimitive,
    {
        /// Create a `SparseVector` from the provided data.
        pub fn from_test_tuples<T: NumCast + Copy>(data: Vec<(usize, T)>, len: usize) -> Self {
            debug_assert!(data.iter().all(|&(i, _)| i < len));
            debug_assert!(data.is_sorted_by_key(|&(i, _)| i));
            debug_assert_ne!(len, 0);
            debug_assert!(data.len() <= len);
            debug_assert!(data.iter().all(|&(_, v)| v.to_f64().unwrap().abs() > EPSILON));

            Self::new(
                data.into_iter().map(|(i, v)| {
                    (i, F::from_f64(v.to_f64().unwrap()).unwrap())
                }).collect(),
                len,
            )
        }
    }

    /// A test vector used in tests.
    fn get_test_vector<F: Field + FromPrimitive, V: TestVector<F>>() -> V {
        V::from_test_data(vec![0, 5, 6])
    }

    /// Test
    fn push_value<F: Field + FromPrimitive, V: TestVector<F>>() where for<'r> &'r F: FieldRef<F> {
        let mut v = get_test_vector::<F, V>();
        let len = v.len();
        let new_v = F!(1);
        v.push_value(new_v.clone());
        assert_eq!(v.len(), len + 1);
        let x: Option<&F> = v.get(v.len() - 1);
        assert_eq!(x, Some(&new_v));
    }

    /// Test
    fn get_set<F: Field + FromPrimitive, V: TestVector<F>>() {
        let mut v = get_test_vector::<F, V>();

        // Getting a nonzero value
        assert_eq!(v.get(1), Some(&F!(5)));

        // Setting to the same value doesn't change
        let value = v.get(2).unwrap().clone();
        v.set(2, value.clone());
        assert_eq!(v.get(2), Some(&value));

        // Changing a value
        let value = F!(3);
        v.set(1, value.clone());
        assert_eq!(v.get(1), Some(&value));

        // Changing a value
        let value = F!(3);
        v.set(0, value.clone());
        assert_eq!(v.get(1), Some(&value));
    }

    /// Test
    fn out_of_bounds_get<F: Field + FromPrimitive, V: TestVector<F>>() {
        let v = get_test_vector::<F, V>();

        &v.get(400);
    }

    /// Test
    fn out_of_bounds_set<F: Field + FromPrimitive, V: TestVector<F>>() {
        let mut v = get_test_vector::<F, V>();

        v.set(400, F!(45));
    }

    /// Test
    fn len<F: Field + FromPrimitive, V: TestVector<F>>() {
        let v = get_test_vector::<F, V>();

        assert_eq!(v.len(), 3);
    }

    #[cfg(test)]
    mod dense_vector {
        use num::traits::FromPrimitive;
        use num::Zero;

        use crate::{F, R32};
        use crate::data::linear_algebra::vector::{DenseVector, Vector};
        use crate::data::linear_algebra::vector::test::{get_set, get_test_vector, len, out_of_bounds_get, out_of_bounds_set, push_value, TestVector};
        use crate::data::number_types::rational::{Rational32, Rational64};
        use crate::data::number_types::traits::Field;

        type T = Rational32;

        fn new<F: Field + FromPrimitive>() {
            let d = vec![0, 5, 6].into_iter().map(|v| F!(v)).collect::<Vec<_>>();
            let len = d.len();
            let v = DenseVector::<F>::new(d, len);

            assert_eq!(v[0], F::zero());
        }

        #[test]
        fn test_new() {
            new::<Rational32>();
            new::<Rational64>();
        }

        #[test]
        fn test_push_value() {
            push_value::<T, DenseVector<T>>();
        }

        #[test]
        fn test_get_set() {
            let v = get_test_vector::<T, DenseVector<_>>();
            // Getting a zero value
            assert_eq!(v.get(0), Some(&T::zero()));

            get_set::<T, DenseVector<T>>();
        }

        #[test]
        #[should_panic]
        fn test_out_of_bounds_get() {
            out_of_bounds_get::<T, DenseVector<T>>();
        }

        #[test]
        #[should_panic]
        fn test_out_of_bounds_set() {
            out_of_bounds_set::<T, DenseVector<T>>();
        }

        #[test]
        fn remove_indices() {
            let mut v = DenseVector::<T>::from_test_data(vec![0, 1, 2]);
            v.remove_indices(&vec![1]);
            assert_eq!(v, DenseVector::<T>::from_test_data(vec![0, 2]));

            let mut v = DenseVector::<T>::from_test_data(vec![3, 0, 0]);
            v.remove_indices(&vec![0]);
            assert_eq!(v, DenseVector::<T>::from_test_data(vec![0, 0]));

            let vs = vec![0, 0, 2, 3, 0, 5, 0, 0, 0, 9];
            let mut v = DenseVector::<T>::from_test_data(vs);
            v.remove_indices(&vec![3, 4, 6]);
            let vs = vec![0, 0, 2, 5, 0, 0, 9];
            assert_eq!(v, DenseVector::<T>::from_test_data(vs));
        }

        #[test]
        fn test_len() {
            len::<T, DenseVector<T>>()
        }

        #[test]
        fn extend_with_values() {
            let mut v = DenseVector::<T>::from_test_data(vec![0, 1, 2]);
            v.extend_with_values(vec![]);
            assert_eq!(v, DenseVector::<T>::from_test_data(vec![0, 1, 2]));

            let mut v = DenseVector::<T>::from_test_data(vec![0, 1, 2]);
            v.extend_with_values(vec![R32!(3)]);
            assert_eq!(v, DenseVector::<T>::from_test_data(vec![0, 1, 2, 3]));

            let mut v = DenseVector::<T>::from_test_data(vec![0, 1, 2]);
            v.extend_with_values(vec![R32!(3), R32!(4)]);
            assert_eq!(v, DenseVector::<T>::from_test_data(vec![0, 1, 2, 3, 4]));
        }
    }

    #[cfg(test)]
    mod sparse_vector {
        use num::FromPrimitive;

        use crate::data::linear_algebra::vector::{SparseVector, Vector};
        use crate::data::linear_algebra::vector::test::{get_set, get_test_vector, len, out_of_bounds_get, out_of_bounds_set, push_value, TestVector};
        use crate::data::number_types::rational::Rational32;
        use crate::R32;

        type T = Rational32;

        #[test]
        fn new() {
            let d = vec![(1, T::from_i32(5).unwrap()), (2, T::from_i32(6).unwrap())];
            let len = 3;
            let v = SparseVector::<T, T>::new(d, len);

            assert_eq!(v.get(0), None);
        }

        #[test]
        fn test_push_value() {
            push_value::<T, SparseVector<T, T>>();
        }

        #[test]
        fn test_get_set() {
            let v = get_test_vector::<T, SparseVector<_, _>>();
            // Getting a zero value
            assert_eq!(v.get(0), None);

            get_set::<T, SparseVector<T, T>>();
        }

        #[test]
        #[should_panic]
        fn test_out_of_bounds_get() {
            out_of_bounds_get::<T, SparseVector<T, T>>();
        }

        #[test]
        #[should_panic]
        fn test_out_of_bounds_set() {
            out_of_bounds_set::<T, SparseVector<T, T>>();
        }

        #[test]
        fn inner_product() {
            let v = get_test_vector::<T, SparseVector<T, T>>();
            let u = get_test_vector::<T, SparseVector<T, T>>();
            assert_eq!(v.inner_product::<Rational32, _>(&u), R32!(5 * 5 + 6 * 6));

            let v = SparseVector::<T, T>::from_test_data(vec![3]);
            let w = SparseVector::<T, T>::from_test_data(vec![5]);
            assert_eq!(v.inner_product::<Rational32, _>(&w), R32!(15));

            let v = SparseVector::<T, T>::from_test_data(vec![0]);
            let w = SparseVector::<T, T>::from_test_data(vec![0]);
            assert_eq!(v.inner_product::<Rational32, _>(&w), R32!(0));

            let v = SparseVector::<T, T>::from_test_data(vec![0, 2]);
            let w = SparseVector::<T, T>::from_test_data(vec![0, 3]);
            assert_eq!(v.inner_product::<Rational32, _>(&w), R32!(6));

            let v = SparseVector::<T, T>::from_test_data(vec![2, 0]);
            let w = SparseVector::<T, T>::from_test_data(vec![0, 3]);
            assert_eq!(v.inner_product::<Rational32, _>(&w), R32!(0));

            let v = SparseVector::<T, T>::from_test_data(vec![0, 2]);
            let w = SparseVector::<T, T>::from_test_data(vec![3, 0]);
            assert_eq!(v.inner_product::<Rational32, _>(&w), R32!(0));

            let v = SparseVector::<T, T>::from_test_data(vec![0, 2]);
            let w = SparseVector::<T, T>::from_test_data(vec![0, 3]);
            assert_eq!(v.inner_product::<Rational32, _>(&w), R32!(6));

            let v = SparseVector::<T, T>::from_test_data(vec![2, 3]);
            let w = SparseVector::<T, T>::from_test_data(vec![5, 7]);
            assert_eq!(v.inner_product::<Rational32, _>(&w), R32!(31));

            let v = SparseVector::<T, T>::from_test_data(vec![0, 0, 0]);
            let w = SparseVector::<T, T>::from_test_data(vec![0, 3, 7]);
            assert_eq!(v.inner_product::<Rational32, _>(&w), R32!(0));

            let v = SparseVector::<T, T>::from_test_data(vec![0, 2, 0]);
            let w = SparseVector::<T, T>::from_test_data(vec![0, 3, 0]);
            assert_eq!(v.inner_product::<Rational32, _>(&w), R32!(6));

            let v = SparseVector::<T, T>::from_test_data(vec![0, 0, 0]);
            let w = SparseVector::<T, T>::from_test_data(vec![5, 7, 0]);
            assert_eq!(v.inner_product::<Rational32, _>(&w), R32!(0));

            let v = SparseVector::<T, T>::from_test_data(vec![0, 0, 2]);
            let w = SparseVector::<T, T>::from_test_data(vec![5, 7, 0]);
            assert_eq!(v.inner_product::<Rational32, _>(&w), R32!(0));

            let v = SparseVector::<T, T>::from_test_data(vec![2, 3, 0]);
            let w = SparseVector::<T, T>::from_test_data(vec![5, 7, 0]);
            assert_eq!(v.inner_product::<Rational32, _>(&w), R32!(31));

            let v = SparseVector::<T, T>::from_test_data(vec![2, 0, 0]);
            let w = SparseVector::<T, T>::from_test_data(vec![5, 7, 0]);
            assert_eq!(v.inner_product::<Rational32, _>(&w), R32!(10));

            let v = SparseVector::<T, T>::from_test_data(vec![0, 2, 0]);
            let w = SparseVector::<T, T>::from_test_data(vec![5, 7, 0]);
            assert_eq!(v.inner_product::<Rational32, _>(&w), R32!(14));

            let v = SparseVector::<T, T>::from_test_data(vec![0, 2, 0]);
            let w = SparseVector::<T, T>::from_test_data(vec![5, 0, 7]);
            assert_eq!(v.inner_product::<Rational32, _>(&w), R32!(0));

            let v = SparseVector::<T, T>::from_test_data(vec![1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0]);
            let w = SparseVector::<T, T>::from_test_data(vec![-1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0]);
            assert_eq!(v.inner_product::<Rational32, _>(&w), R32!(0));
        }

        #[test]
        fn sparse_inner_product() {
            let v = get_test_vector::<T, SparseVector<T, T>>();
            let w = vec![(1, R32!(5)), (2, R32!(6)), ];
            assert_eq!(v.sparse_inner_product::<Rational32, _, _>(w.iter()), R32!(5 * 5 + 6 * 6));

            let v = SparseVector::<T, T>::from_test_data(vec![3]);
            let w = vec![(0, R32!(5))];
            assert_eq!(v.sparse_inner_product::<Rational32, _, _>(w.iter()), R32!(15));

            let v = SparseVector::<T, T>::from_test_data(vec![0]);
            let w = vec![(0, R32!(0))];
            assert_eq!(v.sparse_inner_product::<Rational32, _, _>(w.iter()), R32!(0));

            let v = SparseVector::<T, T>::from_test_data(vec![0, 2]);
            let w = vec![(1, R32!(3))];
            assert_eq!(v.sparse_inner_product::<Rational32, _, _>(w.iter()), R32!(6));

            let v = SparseVector::<T, T>::from_test_data(vec![2, 0]);
            let w = vec![(1, R32!(3))];
            assert_eq!(v.sparse_inner_product::<Rational32, _, _>(w.iter()), R32!(0));

            let v = SparseVector::<T, T>::from_test_data(vec![0, 2]);
            let w = vec![(0, R32!(3))];
            assert_eq!(v.sparse_inner_product::<Rational32, _, _>(w.iter()), R32!(0));

            let v = SparseVector::<T, T>::from_test_data(vec![0, 2]);
            let w = vec![(1, R32!(3))];
            assert_eq!(v.sparse_inner_product::<Rational32, _, _>(w.iter()), R32!(6));

            let v = SparseVector::<T, T>::from_test_data(vec![0, 0, 0]);
            let w = vec![(1, R32!(3)), (2, R32!(7))];
            assert_eq!(v.sparse_inner_product::<Rational32, _, _>(w.iter()), R32!(0));

            let v = SparseVector::<T, T>::from_test_data(vec![0, 2, 0]);
            let w = vec![(0, R32!(5)), (1, R32!(7))];
            assert_eq!(v.sparse_inner_product::<Rational32, _, _>(w.iter()), R32!(14));

            let v = SparseVector::<T, T>::from_test_data(vec![0, 2, 0]);
            let w = vec![(0, R32!(5)), (2, R32!(7))];
            assert_eq!(v.sparse_inner_product::<Rational32, _, _>(w.iter()), R32!(0));

            let v = SparseVector::<T, T>::from_test_data(vec![1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0]);
            let w = vec![(0, R32!(-1)), (1, R32!(1)), (2, R32!(1)), (12, R32!(1))];
            assert_eq!(v.sparse_inner_product::<Rational32, _, _>(w.iter()), R32!(0));
        }

        #[test]
        fn add_multiple_of_row() {
            let mut v = SparseVector::<T, T>::from_test_data(vec![3]);
            let w = SparseVector::<T, T>::from_test_data(vec![5]);
            v.add_multiple_of_row(&R32!(4), &w);
            assert_eq!(v, SparseVector::<T, T>::from_test_data(vec![23]));

            let mut v = SparseVector::<T, T>::from_test_data(vec![0]);
            let w = SparseVector::<T, T>::from_test_data(vec![0]);
            v.add_multiple_of_row(&R32!(4), &w);
            assert_eq!(v, SparseVector::<T, T>::from_test_data(vec![0]));

            let mut v = SparseVector::<T, T>::from_test_data(vec![0, 3]);
            let w = SparseVector::<T, T>::from_test_data(vec![5, 0]);
            v.add_multiple_of_row(&R32!(4), &w);
            assert_eq!(v, SparseVector::<T, T>::from_test_data(vec![20, 3]));

            let mut v = SparseVector::<T, T>::from_test_data(vec![12]);
            let w = SparseVector::<T, T>::from_test_data(vec![3]);
            v.add_multiple_of_row(&R32!(-4), &w);
            assert_eq!(v, SparseVector::<T, T>::from_test_data(vec![0]));

            let mut v = SparseVector::<T, T>::from_test_data(vec![12, 0, 0]);
            let w = SparseVector::<T, T>::from_test_data(vec![0, 0, 3]);
            v.add_multiple_of_row(&R32!(-4), &w);
            assert_eq!(v, SparseVector::<T, T>::from_test_data(vec![12, 0, -12]));
        }

        #[test]
        fn remove_indices() {
            let mut v = SparseVector::<T, T>::from_test_tuples(vec![(0, 3)], 3);
            v.remove_indices(&vec![1]);
            assert_eq!(v, SparseVector::<T, T>::from_test_tuples(vec![(0, 3)], 2));

            let mut v = SparseVector::<T, T>::from_test_tuples(vec![(0, 3)], 3);
            v.remove_indices(&vec![0]);
            assert_eq!(v, SparseVector::<T, T>::from_test_tuples::<f64>(vec![], 2));

            let vs = vec![(2, 2), (3, 3), (5, 5), (10, 10)];
            let mut v = SparseVector::<T, T>::from_test_tuples(vs, 11);
            v.remove_indices(&vec![3, 4, 6]);
            assert_eq!(v, SparseVector::<T, T>::from_test_tuples(vec![(2, 2), (3, 5), (7, 10)], 8));
        }

        #[test]
        fn test_len() {
            len::<T, SparseVector<T, T>>()
        }
    }
}
