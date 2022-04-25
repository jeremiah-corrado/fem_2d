use super::{EigenPair, GEP};
use nalgebra::SymmetricEigen;
use std::fmt;

// TODO: use Nalgebra's Sparse crate
const MAX_DENSE_SIZE: usize = 1000;

/// This function is only recommended in scenarios where the problem size is small and the B-matrix is known to be very well conditioned
///
/// This function directly inverts the B-matrix using Nalgebra's Cholesky Decomposition which does not work well when B is ill-conditioned.
/// It also casts the sparse-matrices as dense matrix objects which uses a very large amount of memory when the matrices are large.
///
/// For larger or more difficult problems the SLEPC Solver is recommended.
pub fn nalgebra_solve_gep(gep: GEP, target_eigenvalue: f64) -> Result<EigenPair, NalgebraGEPError> {
    if gep.a.dimension > MAX_DENSE_SIZE {
        return Err(NalgebraGEPError::ProblemTooLarge);
    }
    let [a_mat, b_mat] = gep.to_nalgebra_dense_mats();
    if let Some(cholesky_decomp) = b_mat.cholesky() {
        let b_inverse = cholesky_decomp.inverse();
        let ba_product = b_inverse * a_mat;
        let ba_se_decomp = SymmetricEigen::new(ba_product);

        if ba_se_decomp.eigenvalues.iter().all(|e| e.abs() < 1e-12) {
            return Err(NalgebraGEPError::SpuriouslyConverged);
        }

        let mut delta_min = f64::MAX;
        let mut delta_min_idx = 0;
        for (eval_idx, eval) in ba_se_decomp.eigenvalues.iter().enumerate() {
            let delta = (eval - target_eigenvalue).abs();
            if delta < delta_min {
                delta_min_idx = eval_idx;
                delta_min = delta;
            }
        }

        Ok(EigenPair {
            value: *ba_se_decomp.eigenvalues.get(delta_min_idx).unwrap(),
            vector: ba_se_decomp
                .eigenvectors
                .column(delta_min_idx)
                .iter()
                .cloned()
                .collect(),
        })
    } else {
        Err(NalgebraGEPError::FailedToInvertB)
    }
}

#[derive(Debug, Clone)]
/// Error type for the SlepcGEP solver
pub enum NalgebraGEPError {
    FailedToInvertB,
    SpuriouslyConverged,
    ProblemTooLarge,
}

impl std::fmt::Display for NalgebraGEPError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::FailedToInvertB => write!(
                f,
                "Failed to invert B-matrix (via cholesky); likely ill-conditioned!"
            ),
            Self::SpuriouslyConverged => write!(f, "Only spurious modes were found!"),
            Self::ProblemTooLarge => write!(
                f,
                "Matrices Exceeded Maximum Size ({}x{}); Cannot Solve!",
                MAX_DENSE_SIZE, MAX_DENSE_SIZE
            ),
        }
    }
}
