/// Use Nalgebra's Eigen decomposition to solve a GEP (not recommended)
pub mod nalgebra_solve;
/// Use External SLEPC solver to solve a GEP
pub mod slepc_solve;
/// Sparsely Packed Matrix
pub mod sparse_matrix;

use nalgebra::DMatrix;
use rayon::prelude::*;
use sparse_matrix::{AIJMatrixBinary, SparseMatrix};
use std::sync::mpsc::channel;

/// Generalized Eigenvalue Problem
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

    pub fn print_to_petsc_binary_files(
        self,
        dir: impl AsRef<str>,
        prefix: impl AsRef<str>,
    ) -> std::io::Result<()> {
        let [a, b]: [AIJMatrixBinary; 2] = [self.a.into(), self.b.into()];
        a.print_to_petsc_binary_file(format!("{}/tmp/{}_a.dat", dir.as_ref(), prefix.as_ref()))?;
        b.print_to_petsc_binary_file(format!("{}/tmp/{}_b.dat", dir.as_ref(), prefix.as_ref()))
    }

    pub fn to_nalgebra_dense_mats(self) -> [DMatrix<f64>; 2] {
        [self.a.into(), self.b.into()]
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

/// Solution to an Eigenvalue Problem
pub struct EigenPair {
    /// Eigenvalue
    pub value: f64,
    /// Eigenvector
    pub vector: Vec<f64>,
}

impl EigenPair {
    /// L2 normalized vector
    pub fn normalized_eigenvector(&self) -> Vec<f64> {
        let norm = self.vector.iter().map(|x| x.powi(2)).sum::<f64>().sqrt();
        self.vector.iter().map(|x| x / norm).collect()
    }
}
