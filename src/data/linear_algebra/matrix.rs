//! # Matrix implementations
//!
//! The `Matrix` trait defines a set of operations available for all matrix types defined in this
//! module. These were written by hand, because a certain specific set of operations needs to be
//! done quickly with these types.
use std::borrow::Borrow;
use std::collections::HashSet;
use std::iter::Iterator;
use std::marker::PhantomData;

use relp_num::Field;
use relp_num::NonZero;

use crate::algorithm::utilities::{remove_indices, remove_sparse_indices};
use crate::data::linear_algebra::{SparseTuple, SparseTupleVec};
use crate::data::linear_algebra::traits::{SparseComparator, SparseElement};

/// Indices start at `0`.
/// TODO(OPTIMIZATION): What data structure is best suited to back this struct? How are allocations
///  avoided (e.g. flattening) avoided?
#[allow(missing_docs)]
#[derive(Eq, PartialEq, Debug, Clone)]
pub struct Sparse<F, C, O: Order> {
    pub data: Vec<SparseTupleVec<F>>,
    major_dimension_size: usize,
    minor_dimension_size: usize,

    phantom_comparison: PhantomData<C>,
    phantom_ordering: PhantomData<O>,
}

/// The backend of a `SparseMatrix` can either be column or row major.
pub trait Order: Sized {
    /// Create a new matrix with this ordering, using specified data.
    ///
    /// Note that the dimensions specified should match the "dimensions" of the data given, and that
    /// there can be no duplicate indices in the same collection of values in the given data.
    fn new<F, C>(
        data: Vec<SparseTupleVec<F>>,
        nr_rows: usize,
        nr_columns: usize,
    ) -> Sparse<F, C, Self>
    where
        F: SparseElement<C>,
        C: SparseComparator,
    ;

    /// Utility method that creates a sparse matrix of this ordering for tests.
    ///
    /// Note that the numerics might not be exact due to intermediate casting to floats, for
    /// convenience in other places of the code base.
    fn from_test_data<F: From<IT>, C, IT: NonZero + Clone>(
        rows: &[Vec<IT>],
        nr_columns: usize,
    ) -> Sparse<F, C, Self>
    where
        F: SparseElement<C>,
        C: SparseComparator,
    ;

    /// Identity matrix of specified field type with this ordering.
    #[must_use]
    fn identity<F, C>(n: usize) -> Sparse<F, C, Self>
    where
        F: Field + Borrow<C> + SparseElement<C>,
        C: SparseComparator,
    {
        Sparse::identity(n)
    }
}

/// Row major sparse matrix ordering.
#[derive(Eq, PartialEq, Copy, Clone, Debug)]
pub struct RowMajor;
/// Column major sparse matrix ordering.
#[derive(Eq, PartialEq, Copy, Clone, Debug)]
pub struct ColumnMajor;
impl Order for RowMajor {
    fn new<F, C>(
        rows: Vec<SparseTupleVec<F>>,
        nr_rows: usize,
        nr_columns: usize,
    ) -> Sparse<F, C, Self>
    where
        F: SparseElement<C>,
        C: SparseComparator,
    {
        Sparse::from_major_ordered_tuples(rows, nr_rows, nr_columns)
    }

    fn from_test_data<F: From<IT>, C, IT: NonZero + Clone>(
        rows: &[Vec<IT>],
        nr_columns: usize,
    ) -> Sparse<F, C, Self>
    where
        F: Borrow<C>,
        C: SparseComparator,
    {
        debug_assert!(rows.iter().all(|v| v.len() == nr_columns));

        let nr_rows = rows.len();

        let mut data: Vec<_> = rows.iter().map(Vec::len).map(Vec::with_capacity).collect();
        for (row_index, row) in rows.iter().enumerate() {
            for (column_index, value) in row.iter().enumerate() {
                if value.is_not_zero() {
                    let new_value = F::from(value.clone());
                    data[row_index].push((column_index, new_value));
                }
            }
        }

        Sparse {
            data,
            major_dimension_size: nr_rows,
            minor_dimension_size: nr_columns,

            phantom_comparison: PhantomData,
            phantom_ordering: PhantomData,
        }
    }
}

impl Order for ColumnMajor {
    fn new<F, C>(
        columns: Vec<SparseTupleVec<F>>,
        nr_rows: usize,
        nr_columns: usize,
    ) -> Sparse<F, C, Self>
    where
        F: SparseElement<C>,
        C: SparseComparator,
    {
        Sparse::from_major_ordered_tuples(columns, nr_columns, nr_rows)
    }

    fn from_test_data<F: From<IT>, C, IT: NonZero + Clone>(
        rows: &[Vec<IT>],
        nr_columns: usize,
    ) -> Sparse<F, C, Self>
    where
        F: SparseElement<C>,
        C: SparseComparator,
    {
        debug_assert!(rows.iter().all(|v| v.len() == nr_columns));

        let nr_rows = rows.len();

        let mut data = vec![vec![]; nr_columns];
        for (row_index, row) in rows.iter().enumerate() {
            for (column_index, value) in row.iter().enumerate() {
                if value.is_not_zero() {
                    let new_value = F::from(value.clone());
                    data[column_index].push((row_index, new_value));
                }
            }
        }

        Sparse {
            data,
            major_dimension_size: nr_columns,
            minor_dimension_size: nr_rows,

            phantom_comparison: PhantomData,
            phantom_ordering: PhantomData,
        }
    }
}

impl<F> Sparse<F, F, RowMajor>
where
    F: SparseElement<F> + SparseComparator,
{
    /// A copy in a different ordering by reference.
    #[must_use]
    pub fn from_column_major_ordered_matrix_although_this_is_expensive(
        data: &Sparse<F, F, ColumnMajor>,
    ) -> Sparse<&F, F, RowMajor> {
        Sparse::from_minor_ordered_tuples(&data.data, data.nr_rows())
    }
}

impl<F, C> Sparse<F, C, RowMajor>
where
    F: SparseElement<C>,
    C: SparseComparator,
{
    /// Concatenate two row oriented matrices.
    ///
    /// The new matrix will have as many rows as the two previous matrices had combined. This is a
    /// cheap operation.
    #[must_use]
    pub fn concatenate_vertically(self, other: Self) -> Self {
        debug_assert_eq!(self.nr_columns(), other.nr_columns());

        self.concatenate_major_indices(other)
    }

    /// Remove rows from this row oriented matrix.
    ///
    /// In contrast to removing columns, this is a very cheap operation.
    pub fn remove_rows(&mut self, indices: &[usize]) {
        self.remove_major_indices(indices);
    }

    /// Remove columns from this row oriented matrix.
    ///
    /// This is an expensive computation that should be avoided almost always, hence the awkward
    /// name.
    pub fn remove_columns_although_this_matrix_is_row_ordered(&mut self, indices: &[usize]) {
        self.remove_minor_indices(indices)
    }

    /// Iterate over all rows by reference.
    pub fn iter_rows(&self) -> impl Iterator<Item = &SparseTupleVec<F>> {
        self.data.iter()
    }

    /// Iterate over a row by reference.
    pub fn iter_row(&self, i: usize) -> impl Iterator<Item = &SparseTuple<F>> {
        debug_assert!(i < self.major_dimension_size);

        self.iter_major_index(i)
    }

    /// Set the value at coordinate (`i`, `j`) to `value`.
    ///
    /// # Arguments
    ///
    /// * `i`: Row index
    /// * `j`: Column index
    /// * `value`: Float that should not be too close to zero to avoid memory usage and numerical
    /// imprecision.
    pub fn set_value(&mut self, row: usize, column: usize, value: F) {
        debug_assert!(row < self.nr_rows());
        debug_assert!(column < self.nr_columns());

        self.inner_set_value(row, column, value)
    }

    /// The number of rows (major dimension size).
    #[must_use]
    pub fn nr_rows(&self) -> usize {
        self.major_dimension_size
    }

    /// The number of columns (minor dimension size).
    #[must_use]
    pub fn nr_columns(&self) -> usize {
        self.minor_dimension_size
    }
}

impl<F: Field, C> Sparse<F, C, ColumnMajor>
where
    F: SparseElement<C>,
    C: SparseComparator,
{
    /// Flip the index of all items in a row.
    ///
    /// # Arguments
    ///
    /// * `rows_to_change`: Indices of rows to change the sign of.
    pub fn change_row_signs(&mut self, rows_to_change: &HashSet<usize>) {
        // TODO(ENHANCEMENT): Consider using a row-major constraint representation to skip a lot of
        //  columns.
        for column in &mut self.data {
            for (i, value) in column {
                if rows_to_change.contains(i) {
                    *value *= -F::one();
                } 
            }
        }
    }
}

impl<F: SparseElement<C>, C: SparseComparator> Sparse<F, C, ColumnMajor> {
    /// Create a row major ordered version of this `SparseMatrix`.
    ///
    /// # Arguments
    ///
    /// * `rows`: List of rows.
    /// * `current_nr_columns`: Amount of columns in the original matrix. Necessary, because the maximum
    /// value in `rows` is not necessarily the last column (could have zero columns at the end).
    ///
    /// # Return value
    ///
    /// A column-major copy.
    #[must_use]
    pub fn from_row_ordered_tuples_although_this_is_expensive(
        rows: &[SparseTupleVec<F>],
        current_nr_columns: usize,
    ) -> Sparse<&F, F, ColumnMajor>
    where
        F: SparseComparator,
    {
        Sparse::from_minor_ordered_tuples(rows, current_nr_columns)
    }

    /// Remove columns from the matrix.
    ///
    /// # Arguments
    ///
    /// * `indices`: Columns to be removed, is assumed sorted.
    pub fn remove_columns(&mut self, indices: &[usize]) {
        debug_assert!(indices.len() <= self.data.len());
        debug_assert!(indices.is_sorted());
        // All values are unique
        debug_assert!(indices.iter().collect::<HashSet<_>>().len() == indices.len());
        debug_assert!(indices.iter().all(|&i| i < self.data.len()));

        self.remove_major_indices(indices);
    }

    /// Remove rows from the matrix.
    ///
    /// # Arguments
    ///
    /// * `indices`: Rows to be removed, is assumed sorted.
    pub fn remove_rows_although_this_matrix_is_column_major(&mut self, indices: &[usize]) {
        debug_assert!(indices.is_sorted());

        self.remove_minor_indices(indices)
    }

    /// Concatenate two matrices.
    ///
    /// # Arguments
    ///
    /// * `other`: A column major ordered sparse matrix with the same number of rows as this
    /// matrix.
    #[must_use]
    pub fn concatenate_horizontally(self, other: Self) -> Self {
        debug_assert_eq!(self.nr_rows(), other.nr_rows());

        self.concatenate_major_indices(other)
    }

    /// Iterating over it's inner `SparseTuple`'s.
    pub fn iter_columns(&self) -> impl Iterator<Item = &SparseTupleVec<F>> {
        self.data.iter()
    }

    /// Get all (`row`, `value`) tuples of column `j`.
    pub fn iter_column(&self, j: usize) -> impl Iterator<Item = &SparseTuple<F>> {
        debug_assert!(j < self.nr_columns());

        self.iter_major_index(j)
    }

    /// Set the value at coordinate (`i`, `j`) to `value`.
    ///
    /// # Arguments
    ///
    /// * `i`: Row index
    /// * `j`: Column index
    /// * `value`: Float that should not be too close to zero to avoid memory usage and numerical
    /// imprecision.
    pub fn set_value(&mut self, row: usize, column: usize, value: F) {
        debug_assert!(column < self.major_dimension_size);
        debug_assert!(row < self.minor_dimension_size);

        self.inner_set_value(column, row, value)
    }

    /// The number of rows (minor dimension size).
    #[must_use]
    pub fn nr_rows(&self) -> usize {
        self.minor_dimension_size
    }

    /// The number of columns (major dimension size).
    #[must_use]
    pub fn nr_columns(&self) -> usize {
        self.major_dimension_size
    }
}

impl<F: SparseElement<C>, C: SparseComparator, MO: Order> Sparse<F, C, MO> {
    /// Create a new instance.
    ///
    /// # Arguments
    ///
    /// * `data`: Collection of `(minor_index, value)` tuples that should already be filtered for
    /// non-zero values.
    /// * `nr_rows`: Number of rows this matrix is large. Couldn't be derived from `columns`,
    /// because the final row(s) might be zero, so no column would have a value in that row.
    /// * `nr_columns`: Number of columns this matrix is large. Could be derived from `columns`,
    /// but not done for consistency. Having to fill in this number often also helps clarify what is
    /// happening in algorithms.
    fn from_major_ordered_tuples(
        data: Vec<SparseTupleVec<F>>,
        major_dimension_size: usize,
        minor_dimension_size: usize,
    ) -> Self {
        debug_assert_eq!(data.len(), major_dimension_size);
        debug_assert!(data.iter().all(|v| v.is_sorted_by_key(|&(i, _)| i)));
        debug_assert!(data.iter()
            .map(|c| c.iter().map(|&(i, _)| i).max())
            .all(|m| m.map_or(true, |max_minor_index| max_minor_index < minor_dimension_size)));
        debug_assert!(data.iter().all(|minor| minor.iter().all(|(_, v)| v.borrow().is_not_zero())));

        Sparse {
            data,
            major_dimension_size,
            minor_dimension_size,

            phantom_comparison: PhantomData,
            phantom_ordering: PhantomData,
        }
    }
}

impl<'a, F, MO> Sparse<&'a F, F, MO>
where
    F: SparseElement<F> + 'a,
    F: SparseComparator, // Implied
    MO: Order,
{
    /// Transpose a sparse matrix.
    pub fn from_minor_ordered_tuples(
        data: &'a [SparseTupleVec<F>],
        current_minor_dimension_size: usize,
    ) -> Self {
        debug_assert!(data.iter().all(|major| major.is_sorted_by_key(|&(i, _)| i)));
        let new_major_dimension_size = current_minor_dimension_size;
        debug_assert!(data.iter().all(|major| major.iter().all(|&(i, _)| i < current_minor_dimension_size)));
        let new_minor_dimension_size = data.len();

        let mut major_ordered = vec![Vec::new(); new_major_dimension_size];
        for (i, minor) in data.iter().enumerate() {
            for (j, value) in minor.iter() {
                major_ordered[*j].push((i, value));
            }
        }

        Sparse::from_major_ordered_tuples(
            major_ordered,
            new_major_dimension_size,
            new_minor_dimension_size,
        )
    }
}

impl<F, C, MO> Sparse<F, C, MO>
where
    F: SparseElement<C>,
    C: SparseComparator,
    MO: Order,
{
    /// Concatenate `SparseMatrix` instances along the major order direction.
    ///
    /// The horizontal direction is the "direction of the rows", that is, the columns will be
    /// stacked together.
    ///
    /// # Arguments
    ///
    /// * `other`: `SparseMatrix` that should have the same number of rows as this matrix.
    ///
    /// # Return value
    ///
    /// `SparseMatrix` with the same number of rows as either input matrix, and the number of
    /// columns equal to the sum of the number of columns of the two matrices.
    fn concatenate_major_indices(self, other: Self) -> Self {
        debug_assert_eq!(other.minor_dimension_size, self.minor_dimension_size);

        Sparse::from_major_ordered_tuples(
            self.data.into_iter().chain(other.data.into_iter()).collect(),
            self.major_dimension_size + other.major_dimension_size,
            self.minor_dimension_size,
        )
    }

    /// Remove columns from the matrix.
    ///
    /// # Arguments
    ///
    /// * `indices`: Columns to be removed, is assumed sorted.
    fn remove_major_indices(&mut self, indices: &[usize]) {
        debug_assert!(indices.len() <= self.major_dimension_size);
        debug_assert!(indices.is_sorted());
        // All values are unique
        debug_assert!(indices.iter().collect::<HashSet<_>>().len() == indices.len());
        debug_assert!(indices.iter().all(|&i| i < self.major_dimension_size));

        remove_indices(&mut self.data, indices);
        self.major_dimension_size -= indices.len();
    }

    /// Remove rows from the matrix.
    ///
    /// # Arguments
    ///
    /// * `indices`: Rows to be removed, is assumed sorted.
    fn remove_minor_indices(&mut self, indices: &[usize]) {
        debug_assert!(indices.len() <= self.minor_dimension_size);
        debug_assert!(indices.is_sorted());
        // All values are unique
        debug_assert!(indices.iter().collect::<HashSet<_>>().len() == indices.len());
        debug_assert!(indices.iter().all(|&i| i < self.minor_dimension_size));

        for j in 0..self.major_dimension_size {
            remove_sparse_indices(&mut self.data[j], indices);
        }
        self.minor_dimension_size -= indices.len();
    }

    fn iter_major_index(&self, major_index: usize) -> impl Iterator<Item = &SparseTuple<F>> {
        debug_assert!(major_index < self.major_dimension_size);

        self.data[major_index].iter()
    }

    /// Set the value at coordinate (`i`, `j`) to `value`.
    ///
    /// # Arguments
    ///
    /// * `i`: Row index
    /// * `j`: Column index
    /// * `value`: Float that should not be too close to zero to avoid memory usage and numerical
    /// imprecision.
    fn inner_set_value(&mut self, major_index: usize, minor_index: usize, value: F) {
        debug_assert!(major_index < self.major_dimension_size);
        debug_assert!(minor_index < self.minor_dimension_size);

        match self.data[major_index].binary_search_by_key(&minor_index, |&(index, _)| index) {
            Ok(index) => self.data[major_index][index].1 = value,
            Err(index) => self.data[major_index].insert(index, (minor_index, value)),
        }
    }

    /// Get the number of non-zero values in this matrix.
    #[must_use]
    pub fn size(&self) -> usize {
        self.data.iter().map(Vec::len).sum()
    }
}

impl<F: Field, C, MO> Sparse<F, C, MO>
where
    F: SparseElement<C>,
    C: SparseComparator,
    MO: Order,
{
    /// Create a dense square identity matrix of size `len`.
    fn identity(len: usize) -> Self {
        debug_assert_ne!(len, 0);

        Sparse::from_major_ordered_tuples(
            (0..len)
                .map(|i| vec![(i, F::one())])
                .collect(),
            len,
            len,
        )
    }

}

#[cfg(test)]
pub mod test {
    use num_traits::FromPrimitive;
    use relp_num::Field;
    use relp_num::R32;
    use relp_num::Rational32;

    use crate::data::linear_algebra::matrix::{ColumnMajor, Order, Sparse};
    use crate::data::linear_algebra::traits::{SparseComparator, SparseElement};

    type T = Rational32;

    fn get_test_matrix<F: From<u8>, C>() -> Sparse<F, C, ColumnMajor>
    where
        F: Field + SparseElement<C>,
        C: SparseComparator,
    {
        ColumnMajor::from_test_data(&vec![
            vec![1, 2, 0],
            vec![0, 5, 6],
        ], 3)
    }

    #[test]
    #[should_panic]
    fn out_of_bounds_set() {
        let mut m = get_test_matrix::<T, T>();

        m.set_value(2, 0, R32!(4));
    }

    #[test]
    fn column_tuples() {
        let m = get_test_matrix::<T, T>();

        assert_eq!(
            m.iter_column(2).nth(0).unwrap(),
            &(1, R32!(6)),
        );
        assert_eq!(
            m.iter_column(1).map(|&(_, value)| value).sum::<T>(),
            T::from_i32(2 + 5).unwrap(),
        );
    }

    #[test]
    fn remove_row() {
        for remove_row in 0..3 {
            let mut data = vec![
                vec![1, 2, 3, 4],
                vec![5, 6, 7, 8],
                vec![9, 10, 11, 12],
            ];
            let mut m = ColumnMajor::from_test_data::<T, T, _>(&data, 4);
            data.remove(remove_row);
            let expected = ColumnMajor::from_test_data::<T, T, _>(&data, 4);

            let mut to_remove = Vec::new();
            to_remove.push(remove_row);
            m.remove_minor_indices(&to_remove);
            let result = m;

            assert_eq!(result, expected);
        }
    }

    #[test]
    fn concatenate_horizontally() {
        let m1 = ColumnMajor::from_test_data::<T, T, _>(&vec![
            vec![1, 2, 3, 4],
            vec![5, 6, 7, 8],
            vec![9, 10, 11, 12],
        ], 4);
        let m2 = ColumnMajor::from_test_data::<T, T, _>(&vec![
            vec![1, 3, 4],
            vec![5, 7, 8],
            vec![9, 11, 12],
        ], 3);
        let result = m1.concatenate_horizontally(m2);
        assert_eq!(result.nr_columns(), 4 + 3);
        assert_eq!(result.nr_rows(), 3);

        let m1 = ColumnMajor::from_test_data::<T, T, _>(&vec![
            vec![1, 2],
            vec![9, 10],
        ], 2);
        let m2 = ColumnMajor::from_test_data::<T, T, _>(&vec![
            vec![5, 7],
            vec![9, 11],
        ], 2);
        let result = m1.concatenate_horizontally(m2);
        assert_eq!(result.nr_columns(), 2 + 2);
        assert_eq!(result.nr_rows(), 2);
    }

    #[test]
    fn remove_columns() {
        // Remove a middle column
        let data = vec![
            vec![1, 2, 3, 4],
            vec![5, 6, 7, 8],
            vec![9, 10, 11, 12],
        ];
        let mut m = ColumnMajor::from_test_data::<T, T, _>(&data, 4);
        let data = vec![
            vec![1, 3, 4],
            vec![5, 7, 8],
            vec![9, 11, 12],
        ];
        let expected = ColumnMajor::from_test_data::<T, T, _>(&data, 3);
        let to_remove = vec![1];
        m.remove_columns(&to_remove);
        let result = m;
        assert_eq!(result, expected);

        // Remove the first column
        let data = vec![
            vec![1, 2, 3, 4],
            vec![5, 6, 7, 8],
            vec![9, 10, 11, 12],
        ];
        let mut m = ColumnMajor::from_test_data::<T, T, _>(&data, 4);
        let data = vec![
            vec![2, 3, 4],
            vec![6, 7, 8],
            vec![10, 11, 12],
        ];
        let expected = ColumnMajor::from_test_data::<T, T, _>(&data, 3);
        let to_remove = vec![0];
        m.remove_columns(&to_remove);
        let result = m;
        assert_eq!(result, expected);

        // Remove the last column
        let data = vec![
            vec![1, 2, 3, 4],
            vec![5, 6, 7, 8],
            vec![9, 10, 11, 12],
        ];
        let mut m = ColumnMajor::from_test_data::<T, T, _>(&data, 4);
        let data = vec![
            vec![1, 2, 3],
            vec![5, 6, 7],
            vec![9, 10, 11],
        ];
        let expected = ColumnMajor::from_test_data::<T, T, _>(&data, 3);
        let to_remove = vec![3];
        m.remove_major_indices(&to_remove);
        let result = m;
        assert_eq!(result, expected);

        // Remove two
        let data = vec![
            vec![1, 2, 3, 4],
            vec![5, 6, 7, 8],
            vec![9, 10, 11, 12],
        ];
        let mut m = ColumnMajor::from_test_data::<T, T, _>(&data, 4);
        let data = vec![
            vec![1, 4],
            vec![5, 8],
            vec![9, 12],
        ];
        let expected = ColumnMajor::from_test_data::<T, T, _>(&data, 2);
        let to_remove = vec![1, 2];
        m.remove_major_indices(&to_remove);
        let result = m;
        assert_eq!(result, expected);
    }

    #[test]
    fn remove_rows() {
        // Remove a middle row
        let data = vec![
            vec![1, 2, 3, 4],
            vec![5, 6, 7, 8],
            vec![9, 10, 11, 12],
        ];
        let mut m = ColumnMajor::from_test_data::<T, T, _>(&data, 4);
        let data = vec![
            vec![1, 2, 3, 4],
            vec![9, 10, 11, 12],
        ];
        let expected = ColumnMajor::from_test_data::<T, T, _>(&data, 4);
        let to_remove = vec![1];
        m.remove_minor_indices(&to_remove);
        let result = m;
        assert_eq!(result, expected);

        // Remove the last row
        let data = vec![
            vec![1, 2, 3, 4],
            vec![5, 6, 7, 8],
            vec![9, 10, 11, 12],
        ];
        let mut m = ColumnMajor::from_test_data::<T, T, _>(&data, 4);
        let data = vec![
            vec![1, 2, 3, 4],
            vec![5, 6, 7, 8],
        ];
        let expected = ColumnMajor::from_test_data::<T, T, _>(&data, 4);
        let to_remove = vec![2];
        m.remove_minor_indices(&to_remove);
        let result = m;
        assert_eq!(result, expected);

        // Remove the first row
        let data = vec![
            vec![1, 2, 3, 4],
            vec![5, 6, 7, 8],
            vec![9, 10, 11, 12],
        ];
        let mut m = ColumnMajor::from_test_data::<T, T, _>(&data, 4);
        let data = vec![
            vec![5, 6, 7, 8],
            vec![9, 10, 11, 12],
        ];
        let expected = ColumnMajor::from_test_data::<T, T, _>(&data, 4);
        let to_remove = vec![0];
        m.remove_minor_indices(&to_remove);
        let result = m;
        assert_eq!(result, expected);

        // Remove two
        let data = vec![
            vec![1, 2, 3, 4],
            vec![5, 6, 7, 8],
            vec![9, 10, 11, 12],
        ];
        let mut m = ColumnMajor::from_test_data::<T, T, _>(&data, 4);
        let data = vec![
            vec![9, 10, 11, 12],
        ];
        let expected = ColumnMajor::from_test_data::<T, T, _>(&data, 4);
        let to_remove = vec![0, 1];
        m.remove_minor_indices(&to_remove);
        let result = m;
        assert_eq!(result, expected);
    }

    #[test]
    fn change_row_signs() {
        let data = vec![
            vec![1, 2, 3, 4],
            vec![5, 6, 7, 8],
            vec![9, 10, 11, 12],
        ];
        let mut m = ColumnMajor::from_test_data(&data, 4);
        let data = vec![
            vec![1, 2, 3, 4],
            vec![-5, -6, -7, -8],
            vec![9, 10, 11, 12],
        ];
        let expected = ColumnMajor::from_test_data::<T, T, _>(&data, 4);
        m.change_row_signs(&vec![1].into_iter().collect());
        assert_eq!(m, expected);
    }

    #[test]
    fn identity() {
        let m = ColumnMajor::identity::<T, T>(1);
        let expected = ColumnMajor::from_test_data::<T, T, _>(&vec![vec![1]], 1);
        assert_eq!(m, expected);

        let m = ColumnMajor::identity::<T, T>(2);
        let expected = ColumnMajor::from_test_data::<T, T, _>(&vec![vec![1, 0], vec![0, 1]], 2);
        assert_eq!(m, expected);
    }
}
