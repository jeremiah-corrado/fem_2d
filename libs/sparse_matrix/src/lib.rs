mod aij_matrix;
mod sparse_matrix_btree;

pub use aij_matrix::AIJMatrix;
pub use sparse_matrix_btree::SparseMatrix;

#[cfg(test)]
mod tests {
    use super::*;

    const VALS_BY_COL: [f64; 6] = [1.0, 2.0, 2.0, 4.0, 3.0, 4.0];

    #[test]
    fn into_aij_representation() {
        let mut sm = SparseMatrix::new(4);
        sm.insert([0, 0], 1.0);
        sm.insert([1, 0], 2.0);
        sm.insert([2, 2], 3.0);
        sm.insert([3, 1], 4.0);

        let aij_mat: AIJMatrix = sm.into();

        assert_eq!(vec![0, 2, 4, 5, 6], aij_mat.i);
        assert_eq!(vec![0, 1, 0, 3, 2, 1], aij_mat.j);

        for (v, v_test) in aij_mat.a.iter().zip(VALS_BY_COL.iter()) {
            assert!((*v_test - *v).abs() < 1e-15);
        }
    }

    #[test]
    fn construct_aij_representation() {
        let mut sm = SparseMatrix::new(4);
        sm.insert([0, 0], 1.0);
        sm.insert([1, 0], 2.0);
        sm.insert([2, 2], 3.0);
        sm.insert([3, 1], 4.0);

        let aij_mat = sm.construct_aij_matrix();

        assert_eq!(sm.num_entries(), VALS_BY_COL.len());

        assert_eq!(vec![0, 2, 4, 5, 6], aij_mat.i);
        assert_eq!(vec![0, 1, 0, 3, 2, 1], aij_mat.j);

        for (v, v_test) in aij_mat.a.iter().zip(VALS_BY_COL.iter()) {
            assert!((*v_test - *v).abs() < 1e-15);
        }
    }
}
