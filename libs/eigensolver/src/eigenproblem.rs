mod sparse_matrix;

use bytes::Bytes;
pub use sparse_matrix::{SparseMatrix, AIJMatrixBinary};

use crate::EigenPair;
use crate::slepc_wrapper::slepc_bridge::AIJMatrix;
use rayon::prelude::*;
use std::sync::mpsc::channel;

use bytes::{BytesMut, Buf};
use std::fs::File;
use std::io::Read;
use std::io::BufReader;


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

pub fn retrieve_solution(prefix: impl AsRef<str>) -> std::io::Result<EigenPair> {
    let evec = retrieve_eigenvector(format!("{}_evec.dat", prefix.as_ref()))?;
    let eval = retrieve_eigenvalue(format!("{}_eval.dat", prefix.as_ref()))?;

    Ok(EigenPair {
        value: eval,
        vector: evec,
    })
}

fn retrieve_eigenvector(path: String) -> std::io::Result<Vec<f64>> {
    let mut vec_file = File::open(path)?;

    let mut header_bytes = BytesMut::new();
    header_bytes.resize(8, 0);
    vec_file.read_exact(&mut header_bytes)?;

    assert_eq!(header_bytes.get_i32(), 1211214_i32); // ensure this is a PETSC vector file
    let m = header_bytes.get_i32() as usize; 

    let mut value_bytes = BytesMut::new();
    value_bytes.resize(m * 8, 0);
    vec_file.read_exact(&mut value_bytes)?;

    let mut values = Vec::with_capacity(m);
    for _ in 0..m {
        values.push(value_bytes.get_f64());
    }

    Ok(values)
}

fn retrieve_eigenvalue(path: String) -> std::io::Result<f64> {
    let mut eval_file = File::open(path)?;
    
    let mut file_bytes = BytesMut::new();
    file_bytes.resize(8, 0);
    eval_file.read_exact(&mut file_bytes)?;

    let value = file_bytes.get_f64();

    Ok(value)
}