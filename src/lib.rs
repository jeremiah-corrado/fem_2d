#![cfg_attr(
    feature = "max_ortho_basis",
    feature(const_fn_floating_point_arithmetic)
)]
#![doc = include_str!("../README.md")]

/// Structures defining the FEM Domain, Mesh, and Basis Space
pub mod fem_domain;

/// Structures and Functions to define and solve the FEM Problem
pub mod fem_problem;

/// Convenient Re-Exports
pub mod prelude {
    pub use crate::fem_domain::basis::hierarchical_basis_fns::poly::HierPoly;
    #[cfg(feature = "max_ortho_basis")]
    pub use crate::fem_domain::basis::shape_fns::max_ortho::MaxOrthoShapeFn;
    pub use crate::fem_domain::domain::{
        dof::{
            basis_spec::{BSAddress, BasisSpec},
            DoF,
        },
        fields::UniformFieldSpace,
        mesh::{
            elem::Elem,
            h_refinement::{HRef, HRefError},
            p_refinement::{PRef, PRefError},
            Mesh,
        },
        ContinuityCondition, Domain,
    };
    pub use crate::fem_problem::galerkin::{galerkin_sample_gep_hcurl, GalerkinSamplingError};
    pub use crate::fem_problem::integration::integrals::{curl_curl::CurlCurl, inner::L2Inner};
    pub use crate::fem_problem::linalg::{
        nalgebra_solve::{nalgebra_solve_gep, NalgebraGEPError},
        slepc_solve::{slepc_solve_gep, SlepcGEPError},
        EigenPair, GEP,
    };
}

#[cfg(test)]
mod tests {
    use super::prelude::*;

    #[test]
    fn nalg_problem() {
        // Define Mesh
        let mut mesh = Mesh::from_file("./test_input/test_mesh_a.json").unwrap();
        mesh.global_p_refinement(PRef::from(2, 2));

        // Construct Domain
        let domain = Domain::from_mesh(mesh, ContinuityCondition::HCurl);
        let ndofs = domain.dofs.len();
        println!("Domain constructed with {} Degrees of Freedom", ndofs);

        // Fill Matrices
        let eigenproblem =
            galerkin_sample_gep_hcurl::<HierPoly, CurlCurl, L2Inner>(&domain, Some([8, 8]))
                .unwrap();

        // Solve Eigenvalue Problem
        let solution = nalgebra_solve_gep(eigenproblem, 2.64).unwrap();
        println!("Found eigenvalue: {:.15}", solution.value);

        assert!((solution.value - 2.6479657_f64).abs() < 1e-6);
        assert_eq!(solution.vector.len(), ndofs);

        let mut field_space = UniformFieldSpace::new(&domain, [8, 8]);
        let e_field_names = field_space
            .xy_fields::<HierPoly>("E", solution.normalized_eigenvector())
            .unwrap();
        field_space
            .expression_2arg(e_field_names, "E_mag", |ex, ey| {
                (ex.powi(2) + ey.powi(2)).sqrt()
            })
            .unwrap();
        field_space
            .print_all_to_vtk("./test_output/mesh_b_fields.vtk")
            .unwrap();
    }

    #[test]
    fn slepc_problem() {
        // Define Mesh
        let mut mesh = Mesh::from_file("./test_input/test_mesh_b.json").unwrap();
        mesh.global_p_refinement(PRef::from(3, 3));
        mesh.global_h_refinement(HRef::T);
        mesh.h_refine_elems(vec![6, 9, 12], HRef::T).unwrap();

        // Construct Domain
        let domain = Domain::from_mesh(mesh, ContinuityCondition::HCurl);
        let ndofs = domain.dofs.len();
        println!("Domain constructed with {} Degrees of Freedom", ndofs);

        // Fill Matrices
        let eigenproblem =
            galerkin_sample_gep_hcurl::<HierPoly, CurlCurl, L2Inner>(&domain, Some([8, 8]))
                .unwrap();

        // Solve Eigenvalue Problem
        let solution = slepc_solve_gep(eigenproblem, 1.475).unwrap();
        println!("Found eigenvalue: {:.15}", solution.value);

        assert!((solution.value - 1.4745880937_f64).abs() < 1e-9);
        assert_eq!(solution.vector.len(), ndofs);

        let mut field_space = UniformFieldSpace::new(&domain, [8, 8]);
        let e_field_names = field_space
            .xy_fields::<HierPoly>("E", solution.normalized_eigenvector())
            .unwrap();
        field_space
            .expression_2arg(e_field_names, "E_mag", |ex, ey| {
                (ex.powi(2) + ey.powi(2)).sqrt()
            })
            .unwrap();
        field_space
            .print_all_to_vtk("./test_output/mesh_b_fields.vtk")
            .unwrap();
    }
}
