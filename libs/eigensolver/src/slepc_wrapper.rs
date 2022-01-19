use super::GEP;

/// Solution to an Eigenvalue Problem
pub struct EigenPair {
    /// Eigenvalue
    pub value: f64,
    /// Eigenvector
    pub vector: Vec<f64>,
}

impl EigenPair {
    /// L2 normalized vector
    pub fn eigenvector_l2(&self) -> Vec<f64> {
        let norm = self.vector.iter().map(|x| x.powi(2)).sum::<f64>().sqrt();
        self.vector.iter().map(|x| x / norm).collect()
    }
}

/// Solve a Generalized Eigenvalue Problem. 
/// 
/// A `target_eigenvalue` is used for spectral decomposition and solution selection
pub fn solve_gep(gep: GEP, target_eigenvalue: f64) -> Result<EigenPair, String> {
    let [a_aij, b_aij] = gep.to_aij_mats();

    let sol = slepc_bridge::slepc_eigenproblem(
        target_eigenvalue,
        a_aij,
        b_aij,
    );

    match sol.status {
        0 => Ok(EigenPair {
            value: sol.eigenvalue,
            vector: Vec::from(sol.eigenvector.as_slice()),
        }),
        i => Err(format!("SLEPC EigenSolver failed with error code: {}", i))
    }
}

#[cxx::bridge(namespace = slepc_wrapper)]
pub mod slepc_bridge {

    struct AIJMatrix {
        pub a: Vec<f64>,
        pub i: Vec<i32>,
        pub j: Vec<i32>,
        pub dim: usize,
    }
    
    struct EigenSolutionInternal {
        status: i32,
        eigenvalue: f64,
        eigenvector: UniquePtr<CxxVector<f64>>,
    }
    
    
    unsafe extern "C++" {
        include!("eigensolver/cpp_src/slepc_wrapper.h");

        fn slepc_eigenproblem(
            target_eigenvalue: f64,
            a_mat: AIJMatrix,
            b_mat: AIJMatrix,
        ) -> EigenSolutionInternal;
    }
}