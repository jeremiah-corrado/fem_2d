mod aij_matrix;
mod sparse_matrix;

pub use aij_matrix::AIJMatrix;
pub use sparse_matrix::SparseMatrix;

use rayon::prelude::*;
use std::sync::mpsc::channel;

/// Generalized Eigenproblem 
/// 
/// Au = λBu
pub struct GEP {
    /// A Matrix
    pub a: SparseMatrix,
    /// B Matrix
    pub b: SparseMatrix,
}

impl GEP {
    pub fn new(num_dofs: usize) -> Self {
        Self {
            a: SparseMatrix::new(num_dofs),
            b: SparseMatrix::new(num_dofs),
        }
    }
}


impl ParallelExtend<[SparseMatrix; 2]> for GEP {
    fn par_extend<I>(&mut self, elem_matrices_iter: I)
    where
        I: IntoParallelIterator<Item = [SparseMatrix; 2]>,
    {
        let (sender, receiver) = channel();

        elem_matrices_iter
            .into_par_iter()
            .for_each_with(sender, |s, mut elem_matrices| {
                s.send([
                    std::mem::replace(&mut elem_matrices[0], SparseMatrix::new(0)),
                    std::mem::replace(&mut elem_matrices[0], SparseMatrix::new(0)),
                ])
                .expect("Failed to send sub-matrices over MSPC channel; cannot construct Matrices!")
            });

        receiver
            .iter()
            .for_each(|[mut elem_a_mat, mut elem_b_mat]| {
                self.a.consume_matrix(&mut elem_a_mat);
                self.b.consume_matrix(&mut elem_b_mat);
            });
    }
}