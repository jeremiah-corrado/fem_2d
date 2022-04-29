use super::{
    integration::HierCurlIntegral,
    linalg::{sparse_matrix::SparseMatrix, GEP},
};
use crate::fem_domain::{
    basis::{BasisFnSampler, HierCurlBasisFn, HierCurlBasisFnSpace},
    domain::{ContinuityCondition, Domain},
};
use rayon::prelude::*;
use std::fmt;

/// Minimum number of Gauss Legendre Quadrature Points Allowed for Galerkin Sampling
pub const MIN_GLQ_ORDER: usize = 4;

/// Fill two system matrices using a [Domain]'s Basis Space as the Testing Space. Return a Generalized Eigenproblem ([GEP])
///
/// All pairs of overlapping Shape Functions will be integrated and stored in the matrices by their associated DoF IDs
///
/// Computations are parallelized over the Rayon Global Threadpool
///
/// # Arguments
/// * `domain`: The [Domain] over which the Galerkin Sampling is to be performed
/// * `glq_grid_dim`: The number of Gauss Legendre Quadrature Points in to use for integration along each direction. If `None`, the default values are used.
/// * Two [HierCurlIntegral]s: `AI` and `BI` must be specified as Generic Arguments. These are used to populate the A and B matrices respectively
/// * A [HierCurlBasisFnSpace] `BSpace` must also be specified as a Generic Argument. This is used to instantiate the Domains `BasisSpec`s as [HierCurlBasisFn]s
///
/// # Returns
/// * An `Err` if the `Domain` was not constructed with an `H(Curl)` [ContinuityCondition]
/// * An `Err` if the `Domain` doesn't have any Degrees of Freedom
/// * An `Err` if the specified number of Gauss Legendre Points is too small
/// * A [GEP], otherwise
///
pub fn galerkin_sample_gep_hcurl<
    BSpace: HierCurlBasisFnSpace,
    AI: HierCurlIntegral,
    BI: HierCurlIntegral,
>(
    domain: &Domain,
    glq_grid_dim: Option<[usize; 2]>,
) -> Result<GEP, GalerkinSamplingError> {
    // check for errors
    if domain.cc != ContinuityCondition::HCurl {
        return Err(GalerkinSamplingError::WrongContinuityCondition(
            ContinuityCondition::HCurl,
            domain.cc,
        ));
    }
    if domain.dofs.is_empty() {
        return Err(GalerkinSamplingError::EmptyDOFSet);
    }
    let [num_glq_u, num_glq_v] = match glq_grid_dim {
        Some([u_dim, v_dim]) => {
            if u_dim < MIN_GLQ_ORDER || v_dim < MIN_GLQ_ORDER {
                return Err(GalerkinSamplingError::InvalidGLQSettings);
            }
            [Some(u_dim), Some(v_dim)]
        }
        None => [None; 2],
    };

    // construct an eigenproblem with a and b matrices
    let mut gep = GEP::new(domain.dofs.len());

    // construct basis sampler
    let [i_max, j_max] = domain.mesh.max_expansion_orders();
    let (bs_sampler, [u_weights, v_weights]): (BasisFnSampler<HierCurlBasisFn<BSpace>>, _) =
        BasisFnSampler::with(i_max as usize, j_max as usize, num_glq_u, num_glq_v, false);

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

    Ok(gep)
}

/// Error Type for Galerkin Sampling Functions
#[derive(Debug)]
pub enum GalerkinSamplingError {
    WrongContinuityCondition(ContinuityCondition, ContinuityCondition),
    EmptyDOFSet,
    InvalidGLQSettings,
}

impl std::error::Error for GalerkinSamplingError {}

impl fmt::Display for GalerkinSamplingError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::WrongContinuityCondition(required, received) => {
                write!(f, "Wrong Continuity Condition on Domain (required: {}, received: {}); Cannot execute Galerkin Sampling!", required, received)
            }
            Self::EmptyDOFSet => write!(
                f,
                "No Degrees-of-Freedom Defined over Domain; Cannot execute Galerkin Sampling!"
            ),
            Self::InvalidGLQSettings => {
                write!(f, "Invalid GLQ Settings (the number of GLQ points must be at least {}); Cannot execute Galerkin Sampling!", MIN_GLQ_ORDER)
            }
        }
    }
}
