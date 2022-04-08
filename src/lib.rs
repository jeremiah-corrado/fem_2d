// ! A 2D FEM implementation which supports aniosotropic hp-refinements
// !
// ! This Library is primarily focused on Solving the Maxwell Eigenvalue Problem on a PEC-terminated Mesh;
// ! however, its functionality can be easily extended by creating custom implementations of
// ! the `Integral` and `ShapeFn` Traits
// !
// ! # Example
// ! This example shows how to solve the Maxwell Eigenvalue Problem on standard waveguide and print the electric fields to a VTK file
// ! ```rust
// ! use fem_2d::prelude::*;
// !
// ! // Load a standard air-filled waveguide mesh from a JSON file
// ! let mut mesh = Mesh::from_file("./test_input/test_mesh_a.json").unwrap();
// !
// ! // Set the polynomial expansion order to 3 in both directions on all Elems
// ! mesh.set_global_expansion_orders([3, 3]);
// !
// ! // Anisotropically refine all Elems along the x-direction
// ! mesh.h_refine_elems(vec![0, 1, 2, 3], HRef::U(None));
// !
// ! // Construct a domain with Dirichlet boundary conditions
// ! let domain = Domain::from_mesh(mesh);
// ! println!("Constructed Domain with {} DoFs", domain.dofs.len());
// !
// ! // Construct a generalized eigenvalue problem for the Electric Field
// !     // (in parallel using the Rayon Global ThreadPool)
// ! let gep = domain.galerkin_sample_gep_parallel::<KOLShapeFn, CurlCurl, L2Inner>(None);
// !
// ! // Solve the generalized eigenvalue problem using Nalgebra's Eigen-Decomposition
// !     // look for an eigenvalue close to 10.0
// ! let solution = nalgebra_solve_gep(gep, 10.0).unwrap();
// ! println!("Found Eigenvalue: {:.15}", solution.value);
// !
// ! // Construct a solution-field-space over the Domain with 64 samples on each "leaf" Elem
// ! let mut field_space = UniformFieldSpace::new(&domain, [8, 8]);
// !
// ! // Compute the Electric Field in the X- and Y-directions (using the same ShapeFns as above)
// ! let e_field_names = field_space.xy_fields::<KOLShapeFn>("E", solution.normalized_eigenvector()).unwrap();
// !
// ! // Compute the magnitude of the Electric Field
// ! field_space.expression_2arg(e_field_names, "E_mag", |ex, ey| (ex.powi(2) + ey.powi(2)).sqrt()).unwrap();
// !
// ! // Print "E_x", "E_y" and "E_mag" to a VTK file
// ! field_space.print_all_to_vtk("./test_output/electric_field_solution.vtk").unwrap();
// !
// ! ````
// !
// ! # Features
// !
// ! * `json_export` (default): Functionality to export Mesh files for visualization with [this](https://github.com/jeremiah-corrado/fem_2d_mesh_plot) tool
// ! * `max_ortho_basis`: A maximally orthogonal basis (requires Nightly toolchain!)

// #![cfg_attr(
//     feature = "max_ortho_basis",
//     feature(const_fn_floating_point_arithmetic)
// )]
// #![doc = include_str!("../README.md")]

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
        fields::UniformFieldSpace,
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
        nalgebra_solve::nalgebra_solve_gep,
        slepc_solve::{slepc_solve_gep, SlepcGEPError},
        EigenPair, GEP,
    };
}

#[cfg(test)]
mod tests {
    use super::prelude::*;

    #[test]
    fn basic_problem() {
        // Define Mesh
        let mut mesh = Mesh::from_file("./test_input/test_mesh_b.json").unwrap();
        mesh.global_p_refinement(PRef::from(3, 3));
        mesh.global_h_refinement(HRef::T);
        mesh.h_refine_elems(vec![6, 9, 12], HRef::T).unwrap();

        // Construct Domain
        let domain = Domain::from_mesh(mesh);
        let ndofs = domain.dofs.len();
        println!("Domain constructed with {} Degrees of Freedom", ndofs);

        // Fill Matrices
        let eigenproblem =
            domain.galerkin_sample_gep_parallel::<KOLShapeFn, CurlCurl, L2Inner>(None);

        // Solve Eigenvalue Problem
        let solution = slepc_solve_gep(eigenproblem, 1.475).unwrap();
        println!("Found eigenvalue: {:.15}", solution.value);

        assert!((solution.value - 1.4745880937_f64).abs() < 1e-9);
        assert_eq!(solution.vector.len(), ndofs);

        let mut field_space = UniformFieldSpace::new(&domain, [8, 8]);
        let e_field_names = field_space
            .xy_fields::<KOLShapeFn>("E", solution.normalized_eigenvector())
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
