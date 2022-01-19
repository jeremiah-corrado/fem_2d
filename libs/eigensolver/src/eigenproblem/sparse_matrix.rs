use super::AIJMatrix;
use std::collections::BTreeMap;

/// Wrapper around a BTreeMap to store square-symmetric matrices in a sparse data structure
#[derive(Clone)]
pub struct SparseMatrix {
    /// Size of the square matrix
    pub dimension: usize,
    /// Matrix Entries
    entries: BTreeMap<[u32; 2], f64>,
}

impl SparseMatrix {
    pub fn new(dimension: usize) -> Self {
        assert!(
            dimension <= (std::u32::MAX as usize),
            "Matrix Dimension cannot exceed the size of a u32!"
        );

        Self {
            dimension,
            entries: BTreeMap::new(),
        }
    }

    pub fn num_entries(&self) -> usize {
        let num_diag = self.entries.keys().filter(|[i, j]| i == j).count();
        2 * self.entries.len() - num_diag
    }

    /// Insert a value into the matrix. Assumes symmetry: row/col order does not matter.
    pub fn insert(&mut self, [row_idx, col_idx]: [usize; 2], value: f64) {
        debug_assert!(
            row_idx < self.dimension,
            "row_idx exceeded matrix dimension; cannot insert value!"
        );
        debug_assert!(
            col_idx < self.dimension,
            "col_idx exceeded matrix dimension; cannot insert value!"
        );

        let coordinates = if row_idx <= col_idx {
            [
                row_idx.try_into().expect("Row Idx was too large!"),
                col_idx.try_into().expect("Col Idx was too large!"),
            ]
        } else {
            [
                col_idx.try_into().expect("Col Idx was too large!"),
                row_idx.try_into().expect("Row Idx was too large!"),
            ]
        };

        if let Some(current_value) = self.entries.get_mut(&coordinates) {
            *current_value += value;
        } else {
            self.entries.insert(coordinates, value);
        }
    }

    // Remove the entries from the matrix, replacing them with an empty BTreeMap.
    fn take_entries(&mut self) -> BTreeMap<[u32; 2], f64> {
        std::mem::replace(&mut self.entries, BTreeMap::new())
    }

    /// Consume the entries from another sparse matrix leaving it empty.
    pub fn consume_matrix(&mut self, other: &mut Self) {
        assert!(
            self.dimension == other.dimension,
            "Sparse Matrices have different dimensions; cannot consume matrix!"
        );
        let new_entries = other.take_entries();

        for (coordinates, value) in new_entries.iter() {
            if let Some(current_value) = self.entries.get_mut(coordinates) {
                *current_value += *value;
            } else {
                self.entries.insert(*coordinates, *value);
            }
        }
    }

    /// Iterate over the upper triangle of the matrix.
    pub fn iter_upper_tri(&self) -> impl Iterator<Item = ([usize; 2], f64)> + '_ {
        self.entries
            .iter()
            .map(|(coords, value)| ([coords[0] as usize, coords[1] as usize], *value))
    }
}

impl Into<AIJMatrix> for SparseMatrix {
    fn into(mut self) -> AIJMatrix {
        // number of entries in each row (indices offset by 1)
        let mut row_counts = vec![0; self.dimension + 1];

        for [r, c] in self.entries.keys() {
            if r == c {
                row_counts[*r as usize + 1] += 1;
            } else {
                row_counts[*r as usize + 1] += 1;
                row_counts[*c as usize + 1] += 1;
            }
        }

        // prefix sum on row_counts
        let mut i = vec![0; self.dimension + 1];
        for (r, r_count) in row_counts.drain(0..).enumerate().skip(1) {
            i[r] = r_count + i[r - 1];
        }

        // upper and lower triangles of matrix; sorted by row then column
        let mut full_matrix: BTreeMap<[u32; 2], f64> = self
            .entries
            .iter()
            .map(|([r, c], v)| ([*c, *r], *v))
            .collect();
        full_matrix.append(&mut self.entries);

        // matrix entries and their associated columns
        let (j, a) = full_matrix
            .iter()
            .map(|([_, c], v)| (*c as i32, *v))
            .unzip();

        AIJMatrix {
            a,
            i,
            j,
            dim: self.dimension,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn value_insertion() {
        let mut sm = SparseMatrix::new(10);

        sm.insert([0, 0], 1.0);
        sm.insert([0, 0], 1.0);
        sm.insert([9, 9], 10.0);
        sm.insert([4, 3], 0.25);
        sm.insert([0, 8], 0.125);
        sm.insert([8, 0], 0.125);

        let raw_entries = sm.take_entries();

        assert!((raw_entries.get(&[0, 0]).unwrap() - 2.0).abs() < 1e-15);
        assert!((raw_entries.get(&[9, 9]).unwrap() - 10.0).abs() < 1e-15);
        assert!((raw_entries.get(&[3, 4]).unwrap() - 0.25).abs() < 1e-15);
        assert!((raw_entries.get(&[0, 8]).unwrap() - 0.25).abs() < 1e-15);

        assert!(raw_entries.get(&[4, 3]).is_none());
        assert!(raw_entries.get(&[8, 0]).is_none());
    }

    #[test]
    fn consume_another_matrix() {
        let mut sm_a = SparseMatrix::new(5);
        let mut sm_b = SparseMatrix::new(5);

        sm_a.insert([0, 0], 1.0);
        sm_a.insert([1, 1], 2.0);
        sm_a.insert([2, 2], 3.0);
        sm_a.insert([3, 3], 4.0);
        sm_a.insert([4, 4], 5.0);
        sm_a.insert([0, 4], 0.5);
        sm_a.insert([3, 1], 0.5);

        sm_b.insert([0, 0], 5.0);
        sm_b.insert([1, 1], 4.0);
        sm_b.insert([2, 2], 3.0);
        sm_b.insert([3, 3], 2.0);
        sm_b.insert([4, 4], 1.0);
        sm_b.insert([4, 0], -0.5);
        sm_b.insert([2, 3], -0.5);

        sm_a.consume_matrix(&mut sm_b);

        assert_eq!(sm_b.num_entries(), 0);

        let sm_a_entries = sm_a.take_entries();

        for i in 0..5 {
            assert!((sm_a_entries.get(&[i, i]).unwrap() - 6.0).abs() < 1e-15);
        }

        assert!((sm_a_entries.get(&[0, 4]).unwrap()).abs() < 1e-15);
        assert!((sm_a_entries.get(&[1, 3]).unwrap() - 0.5).abs() < 1e-15);
        assert!((sm_a_entries.get(&[2, 3]).unwrap() + 0.5).abs() < 1e-15);

        assert!(sm_a_entries.get(&[4, 0]).is_none());
        assert!(sm_a_entries.get(&[3, 1]).is_none());
    }

    #[test]
    #[should_panic]
    fn consume_matrix_of_different_dim() {
        let mut sm_a = SparseMatrix::new(5);
        let mut sm_b = SparseMatrix::new(6);

        sm_a.consume_matrix(&mut sm_b);
    }

    #[test]
    #[should_panic]
    fn oversize_matrix_construction() {
        let _ = SparseMatrix::new((std::u32::MAX as usize) + 1);
    }

    #[test]
    #[should_panic]
    fn out_of_bounds_insertion() {
        let mut sm = SparseMatrix::new(10);
        sm.insert([10, 2], 1.0);
    }
}
