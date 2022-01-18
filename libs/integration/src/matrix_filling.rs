use basis::{BasisFnSampler, ShapeFn};
use fem_domain::Domain;
use sparse_matrix::SparseMatrix;

use super::Integral;

pub fn fill_matrices<AI, BI, SF>(domain: &Domain, ngq: usize) -> [SparseMatrix; 2]
where
    AI: Integral,
    BI: Integral,
    SF: ShapeFn,
{
    let mut a_mat = SparseMatrix::new(domain.dofs.len());
    let mut b_mat = SparseMatrix::new(domain.dofs.len());

    let [i_max, j_max] = domain.mesh.max_expansion_orders();

    let (mut bs_sampler, [u_weights, v_weights]): (BasisFnSampler<SF>, _) =
        BasisFnSampler::with(ngq, ngq, i_max as usize, j_max as usize, false);

    let a_integrator = AI::with_weights(&u_weights, &v_weights);
    let b_integrator = BI::with_weights(&u_weights, &v_weights);

    for elem in domain.elems() {
        let bs_local = bs_sampler.sample_basis_fn(elem, None);

        let local_basis_specs = domain.local_basis_specs(elem.id).unwrap();

        // local - local
        for (i, (p_orders, p_dir, p_dof_id)) in local_basis_specs
            .iter()
            .map(|bs_p| bs_p.get_integration_data())
            .enumerate()
        {
            for (q_orders, q_dir, q_dof_id) in local_basis_specs
                .iter()
                .skip(i)
                .map(|bs_q| bs_q.get_integration_data())
            {
                let a = a_integrator
                    .integrate(p_dir, q_dir, p_orders, q_orders, &bs_local, &bs_local)
                    .surface();
                let b = b_integrator
                    .integrate(p_dir, q_dir, p_orders, q_orders, &bs_local, &bs_local)
                    .surface();

                a_mat.insert([p_dof_id, q_dof_id], a);
                b_mat.insert([p_dof_id, q_dof_id], b);
            }
        }

        let desc_basis_specs = domain.descendant_basis_specs(elem.id).unwrap();

        // local - desc
        for (p_orders, p_dir, p_dof_id) in local_basis_specs
            .iter()
            .map(|bs_p| bs_p.get_integration_data())
        {
            for &(q_elem_id, q_elem_basis_specs) in desc_basis_specs.iter() {
                let bs_p_sampled =
                    bs_sampler.sample_basis_fn(elem, Some(domain.mesh.elem_diag_points(q_elem_id)));
                let bs_q_local = bs_sampler.sample_basis_fn(&domain.mesh.elems[q_elem_id], None);

                for (q_orders, q_dir, q_dof_id) in q_elem_basis_specs
                    .iter()
                    .map(|bs_q| bs_q.get_integration_data())
                {
                    let a = a_integrator
                        .integrate(p_dir, q_dir, p_orders, q_orders, &bs_p_sampled, &bs_q_local)
                        .surface();
                    let b = b_integrator
                        .integrate(p_dir, q_dir, p_orders, q_orders, &bs_p_sampled, &bs_q_local)
                        .surface();

                    a_mat.insert([p_dof_id, q_dof_id], a);
                    b_mat.insert([p_dof_id, q_dof_id], b);
                }
            }
        }
    }

    [a_mat, b_mat]
}
