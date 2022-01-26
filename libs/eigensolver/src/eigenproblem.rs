mod sparse_matrix;

use bytes::Bytes;
pub use sparse_matrix::{SparseMatrix, AIJMatrixBinary};

use crate::EigenPair;
use crate::slepc_wrapper::slepc_bridge::AIJMatrix;
use rayon::prelude::*;
use std::hash::Hash;
use std::sync::mpsc::channel;

use bytes::{BytesMut, Buf};
use std::fs::{File, canonicalize};
use std::io::Read;
use std::env::var_os;
use std::process::{Command, ExitStatus};
use std::time::SystemTime;
use std::hash::Hasher;
use std::collections::hash_map::DefaultHasher;
use std::fmt;

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

    pub fn to_petsc_binary_files(self, dir: impl AsRef<str>, prefix: impl AsRef<str>) -> std::io::Result<()> {
        let [a, b] : [AIJMatrixBinary; 2] = [self.a.into(), self.b.into()];
        a.to_petsc_binary_format(format!("{}/tmp/{}_a.dat", dir.as_ref(), prefix.as_ref()))?;
        b.to_petsc_binary_format(format!("{}/tmp/{}_b.dat", dir.as_ref(), prefix.as_ref()))
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

pub fn solve_eigenproblem(gep: GEP, target_eigenvalue: f64) -> Result<EigenPair, Box<dyn std::error::Error>> {
    if let Some(esolve_dir) = var_os("GEP_SOLVE_DIR") {
        let dir = esolve_dir.to_str().unwrap();
        let prefix = unique_prefix();

        // Write the matrices to files
        gep.to_petsc_binary_files(&dir, &prefix)?;

        // Run the solver
        let esolve_exit_status = Command::new("mpiexec")
            .arg("-np").arg("1")
            .arg("-q")
            .arg("./solve_gep")
            .arg("-te").arg(&format!("{:.10}", target_eigenvalue))
            .arg("-fp").arg(&prefix)
            .current_dir(dir)
            .status();

        match esolve_exit_status {
            Ok(status) => {
                if status.success() {
                    let solution = retrieve_solution(&dir, &prefix)?;
                    clean_directory(&dir, &prefix)?;
                    Ok(solution)
                } else {
                    clean_directory(&dir, &prefix)?;

                    match status.code() {
                        Some(1) => Err(Box::new(EigenSolverError::FailedToInitializeSlepc)),
                        Some(2) => Err(Box::new(EigenSolverError::BadArguments)),
                        Some(3) => Err(Box::new(EigenSolverError::FailedToInitializeMatrices)),
                        Some(4 | 5 | 6) => Err(Box::new(EigenSolverError::FailedToInitializeEPS)),
                        Some(7) => Err(Box::new(EigenSolverError::FailedToConverge)),
                        Some(8) => Err(Box::new(EigenSolverError::FailedToRetreiveSolution)),
                        _ => Err(Box::new(EigenSolverError::UnknownError)),
                    }
                }
            } Err(_) => {
                clean_directory(&dir, &prefix)?;

                Err(Box::new(EigenSolverError::FailedToExecute))
            }
        }
    } else {
        Err(Box::new(EigenSolverError::SolverNotFound))
    }
}

#[derive(Debug, Clone)]
pub enum EigenSolverError {
    SolverNotFound,
    FailedToExecute,
    FailedToInitializeSlepc,
    BadArguments,
    FailedToInitializeMatrices,
    FailedToInitializeEPS,
    FailedToConverge,
    FailedToRetreiveSolution,
    UnknownError
}

impl std::fmt::Display for EigenSolverError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            EigenSolverError::SolverNotFound => write!(f, "Solver not found; please set the GEP_SOLVE_DIR environment variable to the directory containing the solver executable"),
            EigenSolverError::FailedToExecute => write!(f, "Failed to execute solve_gep with MPIEXEC!"),
            EigenSolverError::FailedToInitializeSlepc => write!(f, "Failed to initialize Slepc!"),
            EigenSolverError::BadArguments => write!(f, "Bad arguments passed to solve_gep!"),
            EigenSolverError::FailedToInitializeMatrices => write!(f, "Failed to initialize matrices!"),
            EigenSolverError::FailedToInitializeEPS => write!(f, "Failed to initialize EPS object!"),
            EigenSolverError::FailedToConverge => write!(f, "Failed to converge on the Target Eigenvalue!"),
            EigenSolverError::FailedToRetreiveSolution => write!(f, "Failed to retrieve the solution from solve_gep!"),
            EigenSolverError::UnknownError => write!(f, "Unknown solve_gep error!"),
        }
    }
}

impl std::error::Error for EigenSolverError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        None
    }
}

pub fn retrieve_solution(dir: impl AsRef<str>, prefix: impl AsRef<str>) -> std::io::Result<EigenPair> {
    let evec = retrieve_eigenvector(format!("{}/tmp/{}_evec.dat", dir.as_ref(), prefix.as_ref()))?;
    let eval = retrieve_eigenvalue(format!("{}/tmp/{}_eval.dat", dir.as_ref(), prefix.as_ref()))?;

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

fn unique_prefix() -> String {
    let t_now = SystemTime::now().duration_since(SystemTime::UNIX_EPOCH).unwrap();
    let mut hasher = DefaultHasher::new();
    t_now.hash(&mut hasher);
    format!("p_{}", hasher.finish().to_string().split_at(8).0)
}

fn clean_directory(dir: impl AsRef<str>, prefix: impl AsRef<str>) -> std::io::Result<()> {
    for file in std::fs::read_dir(&format!("{}/tmp/", dir.as_ref()))? {
        let file = file?;
        let file_name_os = file.file_name();
        let file_name = String::from(file_name_os.to_str().unwrap());

        if file_name.contains(prefix.as_ref()) {
            std::fs::remove_file(file.path())?;
        }
    }

    Ok(())
}