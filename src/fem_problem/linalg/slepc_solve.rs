use super::{EigenPair, GEP};
use std::fmt;

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::time::SystemTime;

use bytes::{Buf, BytesMut};
use std::env::var_os;
use std::fs::File;
use std::io::Read;
use std::process::Command;

pub fn slepc_solve_gep(
    gep: GEP,
    target_eigenvalue: f64,
) -> Result<EigenPair, Box<dyn std::error::Error>> {
    if let Some(esolve_dir) = var_os("GEP_SOLVE_DIR") {
        let dir = esolve_dir.to_str().unwrap();
        let prefix = unique_prefix();

        // Write the matrices to files
        gep.print_to_petsc_binary_files(&dir, &prefix)?;

        // Run the solver
        let esolve_exit_status = Command::new("mpiexec")
            .arg("-np")
            .arg("1")
            .arg("-q")
            .arg("./solve_gep")
            .arg("-te")
            .arg(&format!("{:.10}", target_eigenvalue))
            .arg("-fp")
            .arg(&prefix)
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
                        Some(1) => Err(Box::new(SlepcGEPError::FailedToInitializeSlepc)),
                        Some(2) => Err(Box::new(SlepcGEPError::BadArguments)),
                        Some(3) => Err(Box::new(SlepcGEPError::FailedToInitializeMatrices)),
                        Some(4 | 5 | 6) => Err(Box::new(SlepcGEPError::FailedToInitializeEPS)),
                        Some(7) => Err(Box::new(SlepcGEPError::FailedToConverge)),
                        Some(8) => Err(Box::new(SlepcGEPError::FailedToReturnSolution)),
                        _ => Err(Box::new(SlepcGEPError::UnknownError)),
                    }
                }
            }
            Err(_) => {
                clean_directory(&dir, &prefix)?;

                Err(Box::new(SlepcGEPError::FailedToExecute))
            }
        }
    } else {
        println!("Solver not found; please set the GEP_SOLVE_DIR environment variable to the directory containing the solver executable!");
        Err(Box::new(SlepcGEPError::SolverNotFound))
    }
}

#[derive(Debug, Clone)]
/// Error type for the SlepcGEP solver
pub enum SlepcGEPError {
    SolverNotFound,
    FailedToExecute,
    FailedToInitializeSlepc,
    BadArguments,
    FailedToInitializeMatrices,
    FailedToInitializeEPS,
    FailedToConverge,
    FailedToReturnSolution,
    UnknownError,
}

impl std::fmt::Display for SlepcGEPError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::SolverNotFound => write!(f, "Solver not found; please set the GEP_SOLVE_DIR environment variable to the directory containing the solver executable"),
            Self::FailedToExecute => write!(f, "Failed to execute solve_gep with MPIEXEC!"),
            Self::FailedToInitializeSlepc => write!(f, "Slepc failed to initialize!"),
            Self::BadArguments => write!(f, "Bad arguments passed to solve_gep!"),
            Self::FailedToInitializeMatrices => write!(f, "Slepc Failed to initialize matrices!"),
            Self::FailedToInitializeEPS => write!(f, "Slepc Failed to initialize Eigenproblem object!"),
            Self::FailedToConverge => write!(f, "Slepc Failed to converge on the Target Eigenvalue!"),
            Self::FailedToReturnSolution => write!(f, "Slepc Failed to return solution files!"),
            Self::UnknownError => write!(f, "Unknown solve_gep error!"),
        }
    }
}

impl std::error::Error for SlepcGEPError {}

fn retrieve_solution(dir: impl AsRef<str>, prefix: impl AsRef<str>) -> std::io::Result<EigenPair> {
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
    let t_now = SystemTime::now()
        .duration_since(SystemTime::UNIX_EPOCH)
        .unwrap();
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
