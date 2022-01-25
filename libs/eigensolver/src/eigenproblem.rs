mod sparse_matrix;

pub use sparse_matrix::{SparseMatrix, AIJMatrixBinary};

use crate::slepc_wrapper::slepc_bridge::AIJMatrix;
use rayon::prelude::*;
use std::sync::mpsc::channel;

/// Generalized Eigenproblem
///
/// Au = λBu
#[derive(Clone)]
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

    pub(crate) fn to_aij_mats(self) -> [AIJMatrix; 2] {
        [self.a.into(), self.b.into()]
    }

    pub fn to_petsc_binary_files(self, prefix: impl AsRef<str>) -> std::io::Result<()> {
        let [a, b] : [AIJMatrixBinary; 2] = [self.a.into(), self.b.into()];
        a.to_petsc_binary_format(format!("{}_a.dat", prefix.as_ref()))?;
        b.to_petsc_binary_format(format!("{}_b.dat", prefix.as_ref()))
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
            .for_each_with(sender, |s, elem_matrices| {
                s.send(elem_matrices).expect(
                    "Failed to send sub-matrices over MSPC channel; cannot construct Matrices!",
                )
            });

        receiver
            .iter()
            .for_each(|[mut elem_a_mat, mut elem_b_mat]| {
                self.a.consume_matrix(&mut elem_a_mat);
                self.b.consume_matrix(&mut elem_b_mat);
            });
    }
}
