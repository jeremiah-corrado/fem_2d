use super::{
    integration::Integral,
    linalg::{sparse_matrix::SparseMatrix, GEP},
};
use crate::fem_domain::{
    basis::{BasisFnSampler, HierBasisFnSpace},
    domain::Domain,
};
use rayon::prelude::*;

/// Fill two system matrices using a [Domain]'s Basis Space as the Testing Space. Return a Generalized Eigenproblem ([GEP])
///
/// All pairs of overlapping Shape Functions will be integrated and stored in the matrices by their associated DoF IDs
///
/// * Two [Integral]s: `AI` and `BI` must be specified. These are used to populate the A and B matrices respectively
/// * A [ShapeFn] `SF` must also be specified. This is used to evaluate the Domains [BasisSpec]s
///
/// Computations are parallelized over the Rayon Global Threadpool
pub fn galerkin_sample_gep<SF: HierBasisFnSpace, AI: Integral, BI: Integral>(
    domain: &Domain,
    [num_u_glq, num_v_glq]: [Option<usize>; 2],
) -> GEP {
    // construct an eigenproblem with a and b matrices
    let mut gep = GEP::new(domain.dofs.len());

    // construct basis sampler
    let [i_max, j_max] = domain.mesh.max_expansion_orders();
    let (bs_sampler, [u_weights, v_weights]): (BasisFnSampler<SF>, _) =
        BasisFnSampler::with(i_max as usize, j_max as usize, num_u_glq, num_v_glq, false);

    // setup integration
    let a_integrator = AI::with_weights(&u_weights, &v_weights);
    let b_integrator = BI::with_weights(&u_weights, &v_weights);

    gep.par_extend(domain.mesh.elems.par_iter().map(|elem| {
        let mut local_a = SparseMatrix::new(domain.dofs.len());
        let mut local_b = SparseMatrix::new(domain.dofs.len());

        let mut bf_sampler_elem = bs_sampler.clone();
        let elem_materials = elem.get_materials();

        // get relevant data for this Elem
        let bs_local = bf_sampler_elem.sample_basis_fn(elem, None);
        let local_basis_specs = domain.local_basis_specs(elem.id).unwrap();
        let desc_basis_specs = domain.descendant_basis_specs(elem.id).unwrap();

        let mut local_a_entries: Vec<([usize; 2], f64)> =
            Vec::with_capacity(local_basis_specs.len() * local_basis_specs.len() / 2);
        let mut local_b_entries: Vec<([usize; 2], f64)> =
            Vec::with_capacity(local_basis_specs.len() * local_basis_specs.len() / 2);

        // local - local
        for (i, (p_orders, p_dir, p_dof_id)) in local_basis_specs
            .iter()
            .map(|bs_p| bs_p.integration_data())
            .enumerate()
        {
            for (q_orders, q_dir, q_dof_id) in local_basis_specs
                .iter()
                .skip(i)
                .map(|bs_q| bs_q.integration_data())
            {
                let a = a_integrator
                    .integrate(
                        p_dir,
                        q_dir,
                        p_orders,
                        q_orders,
                        &bs_local,
                        &bs_local,
                        elem_materials,
                    )
                    .full_solution();
                let b = b_integrator
                    .integrate(
                        p_dir,
                        q_dir,
                        p_orders,
                        q_orders,
                        &bs_local,
                        &bs_local,
                        elem_materials,
                    )
                    .full_solution();

                local_a_entries.push(([p_dof_id, q_dof_id], a));
                local_b_entries.push(([p_dof_id, q_dof_id], b));
            }
        }

        local_a.insert_group(local_a_entries);
        local_b.insert_group(local_b_entries);

        let mut desc_a_entries: Vec<([usize; 2], f64)> =
            Vec::with_capacity(local_basis_specs.len() * desc_basis_specs.len());
        let mut desc_b_entries: Vec<([usize; 2], f64)> =
            Vec::with_capacity(local_basis_specs.len() * desc_basis_specs.len());

        // local - desc
        for (p_orders, p_dir, p_dof_id) in
            local_basis_specs.iter().map(|bs_p| bs_p.integration_data())
        {
            for &(q_elem_id, q_elem_basis_specs) in desc_basis_specs.iter() {
                let bs_p_sampled =
                    bf_sampler_elem.sample_basis_fn(elem, Some(&domain.mesh.elems[q_elem_id]));
                let bs_q_local =
                    bf_sampler_elem.sample_basis_fn(&domain.mesh.elems[q_elem_id], None);

                for (q_orders, q_dir, q_dof_id) in q_elem_basis_specs
                    .iter()
                    .map(|bs_q| bs_q.integration_data())
                {
                    let a = a_integrator
                        .integrate(
                            p_dir,
                            q_dir,
                            p_orders,
                            q_orders,
                            &bs_p_sampled,
                            &bs_q_local,
                            elem_materials,
                        )
                        .full_solution();
                    let b = b_integrator
                        .integrate(
                            p_dir,
                            q_dir,
                            p_orders,
                            q_orders,
                            &bs_p_sampled,
                            &bs_q_local,
                            elem_materials,
                        )
                        .full_solution();

                    desc_a_entries.push(([p_dof_id, q_dof_id], a));
                    desc_b_entries.push(([p_dof_id, q_dof_id], b));
                }
            }
        }

        local_a.insert_group(desc_a_entries);
        local_b.insert_group(desc_b_entries);

        [local_a, local_b]
    }));

    gep
}
