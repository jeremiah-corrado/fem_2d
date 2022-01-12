use std::collections::BTreeMap;

/// Wrapper around a BTreeMap to store square-symmetric matrices in a sparse data structure
pub struct SparseMatrix {
    /// Size of the square matrix
    pub dimension: usize,
    // Matrix Entries
    entries: BTreeMap<[u32; 2], f64>,
}

impl SparseMatrix {
    pub fn new(dimension: usize) -> Self {
        assert!(
            dimension <= std::mem::size_of::<u32>(),
            "Matrix Dimension cannot exceed the size of a u32!"
        );

        Self {
            dimension,
            entries: BTreeMap::new(),
        }
    }

    /// Insert a value into the matrix. Assumes symmetry: row/col order does not matter.
    pub fn insert(&mut self, [row_idx, col_idx]: [usize; 2], value: f64) {
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

    /// Remove the entries from the matrix, replacing them with an empty BTreeMap.
    /// The map only contains keys located in the upper triangle of the matrix
    pub fn take_entries(&mut self) -> BTreeMap<[u32; 2], f64> {
        std::mem::replace(&mut self.entries, BTreeMap::new())
    }

    /// Insert a series of entries from a BTreeMap.
    /// This function assumes that all key-arrays refer to an entry in the upper triangle (row_idx <= col_idx)
    pub fn insert_entries(&mut self, new_entries: BTreeMap<[u32; 2], f64>) {
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

    /// construct an [AIJMatrix] representation from the values in this matrix
    pub fn to_aij_format(&self) -> AIJMatrix {
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
        full_matrix.append(&mut self.entries.clone());

        // matrix entries and their associated columns
        let (j, a) = full_matrix.iter().map(|([_, c], v)| {
            (*c as i32, *v)
        }).unzip();

        AIJMatrix {
            a,
            i,
            j,
            dim: self.dimension,
        }
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
        let (j, a) = full_matrix.iter().map(|([_, c], v)| {
            (*c as i32, *v)
        }).unzip();

        AIJMatrix {
            a,
            i,
            j,
            dim: self.dimension,
        }
    }
}

pub struct AIJMatrix {
    pub a: Vec<f64>,
    pub i: Vec<i32>,
    pub j: Vec<i32>,
    pub dim: usize,
}
