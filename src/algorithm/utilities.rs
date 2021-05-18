//! # Utilities
//!
//! Helper functions for algorithms.
use std::cmp::Ordering;
use std::collections::HashSet;

use crate::data::number_types::nonzero::Nonzero;

/// Reduce the size of the vector by removing values.
///
/// # Arguments
///
/// * `vector`: `Vec` to remove indices from.
/// * `indices`: A set of indices to remove from the vector, assumed sorted.
pub(crate) fn remove_indices<T>(vector: &mut Vec<T>, indices: &[usize]) {
    debug_assert!(indices.len() <= vector.len());
    debug_assert!(indices.is_sorted());
    // All values are unique
    debug_assert!(indices.iter().collect::<HashSet<_>>().len() == indices.len());
    debug_assert!(indices.iter().all(|&i| i < vector.len()));

    let mut i = 0;
    let mut j = 0;
    vector.retain(|_| {
        if i < indices.len() && j < indices[i] {
            j += 1;
            true
        } else if i < indices.len() { // must have j == to_remove[i]
            j += 1;
            i += 1;
            false
        } else { // i == to_remove.len()
            j += 1;
            true
        }
    });
}

/// Reduce the size of the vector by removing values.
///
/// There is another version of this algorithm implemented on `DenseVector` for `Vec<T>`, where `T`
/// is not necessarily `Copy`.
///
/// The method operates in place.
///
/// # Arguments
///
/// * `vector`: `Vec` to remove indices from.
/// * `indices`: A set of indices to remove from the vector, assumed sorted.
pub(crate) fn remove_sparse_indices<T>(vector: &mut Vec<(usize, T)>, indices: &[usize]) {
    debug_assert!(indices.is_sorted());
    // All values are unique
    debug_assert!(indices.iter().collect::<HashSet<_>>().len() == indices.len());

    if indices.is_empty() || vector.is_empty() {
        return;
    }

    let mut nr_skipped_before = 0;
    vector.drain_filter(|(i, _)| {
        while nr_skipped_before < indices.len() && indices[nr_skipped_before] < *i {
            nr_skipped_before += 1;
        }

        if nr_skipped_before < indices.len() && indices[nr_skipped_before] == *i {
            true
        } else {
            *i -= nr_skipped_before;
            false
        }
    });
}

#[cfg(test)]
mod test {
    use std::convert::identity;
    use std::ops::Add;

    use crate::algorithm::utilities::{merge_sparse_indices, remove_indices, remove_sparse_indices};

    #[test]
    fn test_remove_indices() {
        let mut v = vec![0f64, 1f64, 2f64];
        remove_indices(&mut v, &vec![1]);
        assert_eq!(v, vec![0f64, 2f64]);

        let mut v = vec![3f64, 0f64, 0f64];
        remove_indices(&mut v, &vec![0]);
        assert_eq!(v, vec![0f64, 0f64]);

        let mut v = vec![0f64, 0f64, 2f64, 3f64, 0f64, 5f64, 0f64, 0f64, 0f64, 9f64];
        remove_indices(&mut v,&vec![3, 4, 6]);
        assert_eq!(v, vec![0f64, 0f64, 2f64, 5f64, 0f64, 0f64, 9f64]);

        let mut v = vec![0f64];
        remove_indices(&mut v, &vec![0]);
        assert_eq!(v, vec![]);

        let mut v: Vec<i32> = vec![];
        remove_indices(&mut v, &vec![]);
        assert_eq!(v, vec![]);
    }

    #[test]
    fn test_remove_sparse_indices() {
        // Removing value not present
        let mut tuples = vec![(0, 3f64)];
        let indices = vec![1];
        remove_sparse_indices(&mut tuples, &indices);
        assert_eq!(tuples, vec![(0, 3f64)]);

        // Removing present value, index should be adjusted
        let mut tuples = vec![(0, 0f64), (2, 2f64)];
        let indices = vec![0];
        remove_sparse_indices(&mut tuples, &indices);
        assert_eq!(tuples, vec![(1, 2f64)]);

        // Empty vec
        let mut tuples: Vec<(usize, i32)> = vec![];
        let indices = vec![0, 1, 1_000];
        remove_sparse_indices(&mut tuples, &indices);
        assert_eq!(tuples, vec![]);

        // Empty vec, removing nothing
        let mut tuples: Vec<(usize, i32)> = vec![];
        let indices = vec![];
        remove_sparse_indices(&mut tuples, &indices);
        assert_eq!(tuples, vec![]);

        // Non-empty vec, removing nothing
        let mut tuples = vec![(1_000, 1f64)];
        let indices = vec![];
        remove_sparse_indices(&mut tuples, &indices);
        assert_eq!(tuples, vec![(1_000, 1f64)]);
    }

    #[test]
    fn test_merge_sparse_indices() {
        // Empty
        let left: Vec<(i8, i16)> = vec![];
        let right = vec![];

        let result = merge_sparse_indices(left.into_iter(), right.into_iter(), Add::add, identity, identity);
        let expected = vec![];
        assert_eq!(result, expected);

        // One empty
        let left: Vec<(i8, i16)> = vec![(2, 1)];
        let right = vec![];

        let result = merge_sparse_indices(left.into_iter(), right.into_iter(), Add::add, identity, identity);
        let expected = vec![(2, 1)];
        assert_eq!(result, expected);

        // Not related
        let left = vec![(1, 6)].into_iter();
        let right = vec![(4, 9)].into_iter();

        let result = merge_sparse_indices(left, right, Add::add, identity, identity);
        let expected = vec![(1, 6), (4, 9)];
        assert_eq!(result, expected);

        // Same index
        let left = vec![(1, 6)].into_iter();
        let right = vec![(1, 9)].into_iter();

        let result = merge_sparse_indices(left, right, Add::add, identity, identity);
        let expected = vec![(1, 15)];
        assert_eq!(result, expected);

        // Alternating
        let left = vec![(1, 6), (3, 4)].into_iter();
        let right = vec![(2, 9)].into_iter();

        let result = merge_sparse_indices(left, right, Add::add, identity, identity);
        let expected = vec![(1, 6), (2, 9), (3, 4)];
        assert_eq!(result, expected);
    }
}
