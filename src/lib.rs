#![cfg_attr(
    feature = "max_ortho_basis",
    feature(const_fn_floating_point_arithmetic)
)]

extern crate bytes;
extern crate nalgebra;
extern crate num_complex;
extern crate rayon;

#[macro_use]
extern crate smallvec;

#[macro_use]
extern crate json;

/// Structures and Traits for Basis Function Evaluation
pub mod basis;
/// Structures to define the geometric characteristics and refinement state of an FEM Domain
pub mod domain;
/// Structures and functions to assist in the integration of Basis Functions
pub mod integration;
/// Structures and functions to solve Generalized Eigenvalue Problems
pub mod linalg;

/// Convenient Re-Exports
pub mod prelude {
    pub use crate::basis::shape_fns::kol::KOLShapeFn;
    #[cfg(feature = "max_ortho_basis")]
    pub use crate::basis::shape_fns::max_ortho::MaxOrthoShapeFn;
    pub use crate::domain::{
        dof::{
            basis_spec::{BSAddress, BasisSpec},
            DoF,
        },
        mesh::{
            elem::Elem,
            h_refinement::{HRef, HRefError},
            p_refinement::{PRef, PRefError},
            Mesh,
        },
        Domain,
    };
    pub use crate::integration::integrals::{curl_curl::CurlCurl, inner::L2Inner};
    pub use crate::linalg::{
        slepc_gep_link::{slepc_solve_gep, EigenSolverError},
        EigenPair, GEP,
    };
}

#[cfg(test)]
mod tests {
    use super::prelude::*;

    #[test]
    fn basic_problem() {
        let mut mesh = Mesh::from_file("./test_input/test_mesh_b.json").unwrap();
        mesh.global_p_refinement(PRef::from(3, 3)).unwrap();
        mesh.global_h_refinement(HRef::T).unwrap();
        mesh.h_refine_elems(vec![6, 9, 12], HRef::T).unwrap();

        let domain = Domain::from_mesh(mesh);
        let ndofs = domain.dofs.len();
        println!("Domain constructed with {} Degrees of Freedom", ndofs);

        let eigenproblem =
            domain.galerkin_sample_gep_parallel::<KOLShapeFn, CurlCurl, L2Inner>(None);
        let solution = slepc_solve_gep(eigenproblem, 1.475).unwrap();
        println!("Found eigenvalue: {:.15}", solution.value);

        assert!((solution.value - 1.4745880937_f64).abs() < 1e-9);
        assert_eq!(solution.vector.len(), ndofs);
    }
}
