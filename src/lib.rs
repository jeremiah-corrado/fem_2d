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

#[cfg(test)]
mod tests {
    // use super::*;
    // use rayon::prelude::*;

    // #[test]
    // fn integration_correctness() {
    //     let i_max: usize = 2;
    //     let j_max: usize = 2;
    //     let p_elem_id = 1;
    //     let q_elem_id = 8;

    //     let mut se_mesh = Mesh::from_file("./test_input/test_mesh_c.json").unwrap();
    //     se_mesh.set_expansion_on_elems(vec![0], [i_max as u8, j_max as u8]).unwrap();
    //     se_mesh.h_refine_elems(vec![0], HRef::T).unwrap();
    //     se_mesh.h_refine_elems(vec![1], HRef::U(None)).unwrap();
    //     se_mesh.h_refine_elems(vec![6], HRef::V(None)).unwrap();

    //     se_mesh.export_to_json("./test_output/integration_test_mesh.json").unwrap();

    //     let (mut bs_sampler, [u_weights, v_weights]): (BasisFnSampler<KOLShapeFn>, _) =
    //         BasisFnSampler::with(i_max, j_max, None, None, false);

    //     let a_integ = CurlProduct::with_weights(&u_weights, &v_weights);
    //     let b_integ = L2InnerProduct::with_weights(&u_weights, &v_weights);

    //     println!("P:");
    //     let bs_p = bs_sampler.sample_basis_fn(&se_mesh.elems[p_elem_id], match q_elem_id == p_elem_id {
    //         false => Some(&se_mesh.elems[q_elem_id]),
    //         true => None,
    //     });
    //     println!("Q:");
    //     let bs_q = bs_sampler.sample_basis_fn(&se_mesh.elems[q_elem_id], None);

    //     for p_dir in [BasisDir::U, BasisDir::V] {
    //         for q_dir in [BasisDir::U, BasisDir::V] {
    //             for p_i in 0..=i_max {
    //                 for p_j in 0..=j_max {
    //                     for q_i in 0..=i_max {
    //                         for q_j in 0..=j_max {

    //                             let a = a_integ
    //                                 .integrate(p_dir, q_dir, [p_i, p_j], [q_i, q_j], &bs_p, &bs_q)
    //                                 .full_solution();
    //                             let b = b_integ
    //                                 .integrate(p_dir, q_dir, [p_i, p_j], [q_i, q_j], &bs_p, &bs_q)
    //                                 .full_solution();

    //                             println!(
    //                                 "[{}, ({}, {})] \t [{}, ({}, {})] \t {:.10} \t {:.10}",
    //                                 p_dir,
    //                                 p_i,
    //                                 p_j,
    //                                 q_dir,
    //                                 q_i,
    //                                 q_j,
    //                                 a,
    //                                 b
    //                             );
    //                         }
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }
    // #[test]
    // fn basic_problem_seq() {
    //     let mut mesh = Mesh::from_file("./test_input/test_mesh_b.json").unwrap();
    //     mesh.global_p_refinement(PRef::from(3, 3)).unwrap();
    //     mesh.global_h_refinement(HRef::T).unwrap();
    //     mesh.h_refine_elems(vec![6, 9, 12], HRef::T).unwrap();

    //     let domain = Domain::from_mesh(mesh);

    //     let eigenproblem = fill_matrices::<CurlProduct, L2InnerProduct, KOLShapeFn>(&domain);
    //     let eigen_pair = solve_eigenproblem(eigenproblem, 1.475).unwrap();

    //     assert!((eigen_pair.value - 1.4745880937_f64).abs() < 1e-9);
    // }

    // #[test]
    // fn basic_problem_par() {
    //     rayon::ThreadPoolBuilder::new()
    //         .num_threads(4)
    //         .build_global()
    //         .unwrap();

    //     let mut mesh = Mesh::from_file("./test_input/test_mesh_b.json").unwrap();
    //     mesh.global_p_refinement(PRef::from(3, 3)).unwrap();
    //     mesh.global_h_refinement(HRef::T).unwrap();
    //     mesh.h_refine_elems(vec![6, 9, 12], HRef::T).unwrap();

    //     let domain = Domain::from_mesh(mesh);

    //     println!("NDOFS: {}", domain.dofs.len());

    //     let eigenproblem =
    //         fill_matrices_parallel::<CurlProduct, L2InnerProduct, KOLShapeFn>(&domain);
    //     let eigen_pair = solve_eigenproblem(eigenproblem, 1.475).unwrap();

    //     assert!((eigen_pair.value - 1.4745880937_f64).abs() < 1e-9);
    // }
}
