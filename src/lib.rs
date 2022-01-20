extern crate basis;
extern crate fem_domain;
extern crate eigensolver;
extern crate integration;

pub use basis::{BasisFn, BasisFnSampler, KOLShapeFn, MaxOrthoShapeFn, ShapeFn};
pub use fem_domain::{Domain, Mesh, Point, M2D, V2D, DoF, HRef, PRef, HRefError, PRefError};
pub use integration::{CurlProduct, Integral, IntegralResult, L2InnerProduct, fill_matrices, fill_matrices_parallel};
pub use eigensolver::{GEP, SparseMatrix, solve_gep, EigenPair};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn basic_problem() {
        let mut domain = Domain::from_mesh_file("./test_input/test_mesh_b.json").unwrap();

        domain.mesh.global_p_refinement(PRef::from(1, 1)).unwrap();
        // domain.mesh.global_h_refinement(HRef::T).unwrap();
        domain.gen_dofs();

        println!("Num DoFs: {}", domain.dofs.len());
        for dof in domain.dofs.iter() {
            println!("{}", dof);
        }

        let eigenproblem = fill_matrices::<CurlProduct, L2InnerProduct, MaxOrthoShapeFn>(&domain);
        let eigen_pair = solve_gep(eigenproblem, 1.475).unwrap();

        println!("Solution: {:.10}", eigen_pair.value);
    }
}
