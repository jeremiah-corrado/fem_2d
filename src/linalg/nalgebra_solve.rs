use super::{EigenPair, GEP};
use nalgebra::SymmetricEigen;

/// This function is only recommended in scenarios where the problem size is small and the B-matrix is known to be very well conditioned
///
/// This function directly inverts the B-matrix using Nalgebra's Cholesky Decomposition which does not work well when B is ill-conditioned.
/// It also casts the sparse-matrices as dense matrix objects which uses a very large amount of memory for large matrices.
///
/// For larger or more difficult problems the SLEPC Solver is recommended.
pub fn nalgebra_solve_gep(gep: GEP, target_eigenvalue: f64) -> Result<EigenPair, String> {
    let [a_mat, b_mat] = gep.to_nalgebra_dense_mats();
    if let Some(cholesky_decomp) = b_mat.cholesky() {
        let b_inverse = cholesky_decomp.inverse();
        let ba_product = b_inverse * a_mat;
        let ba_se_decomp = SymmetricEigen::new(ba_product);

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
        Err(String::from("Unable to compute inverse of B Matrix!"))
    }
}
