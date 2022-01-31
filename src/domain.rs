/// Degrees of Freedom
pub mod dof;
/// Structures used to compute solution fields over a Domain
pub mod fields;
/// The internal geometric structure of a Domain
pub mod mesh;

use crate::basis::{BasisFnSampler, ParBasisFnSampler, ShapeFn};
use crate::integration::Integral;
use crate::linalg::{sparse_matrix::SparseMatrix, GEP};
use dof::{
    basis_spec::{BSAddress, BasisDir, BasisLoc, BasisSpec},
    DoF,
};
use mesh::*;

use rayon::prelude::*;
use smallvec::smallvec;
use std::collections::BTreeMap;
use std::sync::{Arc, Mutex};

/// High Level Description of an FEM Domain
pub struct Domain {
    pub mesh: Mesh,
    /// Degrees of Freedom (collections of BasisSpecs matched via Tangential Continuity)
    pub dofs: Vec<DoF>,
    /// Individual BasisSpecs sorted by their associated Elem
    pub basis_specs: Vec<Vec<BasisSpec>>,
}

impl Domain {
    pub fn blank() -> Self {
        Self {
            mesh: Mesh::blank(),
            dofs: Vec::new(),
            basis_specs: Vec::new(),
        }
    }

    pub fn from_mesh(mesh: Mesh) -> Self {
        let mut dom = Self {
            mesh,
            dofs: Vec::new(),
            basis_specs: Vec::new(),
        };

        dom.gen_dofs();

        dom
    }

    /// Iterate over all `Elem`s in the mesh
    pub fn elems<'a>(&'a self) -> impl Iterator<Item = &'a mesh::elem::Elem> + '_ {
        self.mesh.elems.iter()
    }

    /// Iterate over all `Edge`s in the mesh
    pub fn edges<'a>(&'a self) -> impl Iterator<Item = &'a mesh::edge::Edge> + '_ {
        self.mesh.edges.iter()
    }

    /// Iterate over all `Node`s in the mesh
    pub fn nodes<'a>(&'a self) -> impl Iterator<Item = &'a mesh::node::Node> + '_ {
        self.mesh.nodes.iter()
    }

    // Generate Degrees of Freedom over the mesh according to the Polynomial Expansion orders on each Elem
    fn gen_dofs(&mut self) {
        // prepare for fresh set of DoFs and BasisSpecs
        self.basis_specs = vec![Vec::new(); self.mesh.elems.len()];
        self.dofs.clear();
        self.mesh.set_edge_activation();
        let mut dof_id_tracker = IdTracker::new(0);

        // Generate lists of BasisSpecs associated with Elems, Edges, and Nodes, sorted by their IDs
        let [elem_bs, edge_bs, _] = self.gen_basis_specs();

        // Designate all elem-type BasisSpecs located on shell Elems as DoFs
        for (elem_id, mut elem_bs_list) in elem_bs {
            if !self.mesh.elems[elem_id].has_children() {
                self.basis_specs[elem_id] = Vec::with_capacity(elem_bs_list.len());

                for elem_bs in elem_bs_list
                    .drain(0..)
                    .filter(|bs| bs.dir == BasisDir::U || bs.dir == BasisDir::V)
                {
                    let dof_id = dof_id_tracker.next_id();
                    let address = self.push_basis_spec(elem_bs, dof_id);
                    self.dofs.push(DoF::new(dof_id, smallvec![address]));
                }
            }
        }

        // Create DoFs from pairs of matched BasisSpecs on the active Elems associated with each Edge
        for (edge_id, mut edge_bs_list) in edge_bs {
            if let Some(active_elem_ids) = self.mesh.edges[edge_id].active_elem_pair() {
                // only basis specs associated with the active pair of Elems need to be considered here
                let rel_basis_specs: Vec<BasisSpec> = edge_bs_list
                    .drain(0..)
                    .filter(|bs| {
                        (bs.dir == BasisDir::U || bs.dir == BasisDir::V)
                            && active_elem_ids.contains(&bs.elem_id)
                    })
                    .collect();

                // allocate space for the new basis specs
                let num_expected = rel_basis_specs.len() / 2;
                for elem_id in active_elem_ids {
                    if self.basis_specs[elem_id].is_empty() {
                        self.basis_specs[elem_id] = Vec::with_capacity(num_expected);
                    } else {
                        self.basis_specs[elem_id].reserve(num_expected);
                    }
                }

                // iterate over each pair of BasisSpecs (once) and look for matches
                let mut active_pairs: Vec<[usize; 2]> = Vec::with_capacity(num_expected);
                for (a, bs_0) in rel_basis_specs.iter().enumerate() {
                    for (b, bs_1) in rel_basis_specs.iter().enumerate().skip(a + 1) {
                        if bs_0.matches_with_edge(bs_1) {
                            active_pairs.push([a, b]);
                            break;
                        }
                    }
                }

                // Store the matched BasisSpecs and create new DoFs
                for pair in active_pairs {
                    let dof_id = dof_id_tracker.next_id();
                    let addresses = pair
                        .iter()
                        .map(|rel_idx| {
                            // TODO: should use MaybeUninit in BasisSpec (or some other method) to avoid expensive Clone  here!
                            self.push_basis_spec(rel_basis_specs[*rel_idx].clone(), dof_id)
                        })
                        .collect();

                    self.dofs.push(DoF::new(dof_id, addresses));
                }
            }
        }

        // TODO: implement node-type BasisSpec Matching!
    }

    fn gen_basis_specs(&self) -> [BTreeMap<usize, Vec<BasisSpec>>; 3] {
        let mut elem_bs: BTreeMap<usize, Vec<BasisSpec>> = BTreeMap::new();
        let mut edge_bs: BTreeMap<usize, Vec<BasisSpec>> = BTreeMap::new();
        let mut node_bs: BTreeMap<usize, Vec<BasisSpec>> = BTreeMap::new();

        let mut bs_id_tracker = IdTracker::new(0);

        for elem in self.elems() {
            for dir in [BasisDir::U, BasisDir::V, BasisDir::W] {
                for poly_ij in elem.poly_orders.permutations(dir) {
                    let bs = BasisSpec::new(bs_id_tracker.next_id(), poly_ij, dir, elem);

                    match bs.loc {
                        BasisLoc::ElemBs => elem_bs
                            .entry(elem.id)
                            .and_modify(|bs_list| bs_list.push(bs.clone()))
                            .or_insert_with(|| vec![bs]),
                        BasisLoc::EdgeBs(_, edge_id) => edge_bs
                            .entry(edge_id)
                            .and_modify(|bs_list| bs_list.push(bs.clone()))
                            .or_insert_with(|| vec![bs]),
                        BasisLoc::NodeBs(_, node_id) => node_bs
                            .entry(node_id)
                            .and_modify(|bs_list| bs_list.push(bs.clone()))
                            .or_insert_with(|| vec![bs]),
                    };
                }
            }
        }

        [elem_bs, edge_bs, node_bs]
    }

    /// Retrieve a list of [BasisSpec]s on an `Elem` by ID
    pub fn local_basis_specs(&self, elem_id: usize) -> Result<&Vec<BasisSpec>, String> {
        if elem_id >= self.mesh.elems.len() {
            Err(format!(
                "Elem {} doesn't exist; Cannot retrieve BasisSpecs!",
                elem_id
            ))
        } else {
            Ok(&self.basis_specs[elem_id])
        }
    }

    /// Retrieve a list of an `Elem`s descendant [BasisSpec]s (All the [`BasisSpec`]s on its descendant `Elem`s)
    pub fn descendant_basis_specs(
        &self,
        elem_id: usize,
    ) -> Result<Vec<(usize, &Vec<BasisSpec>)>, String> {
        if elem_id >= self.mesh.elems.len() {
            Err(format!(
                "Elem {} doesn't exist; Cannot retrieve Descendant BasisSpecs!",
                elem_id
            ))
        } else {
            let desc_elem_ids = self.mesh.descendant_elems(elem_id, false)?;
            Ok(desc_elem_ids
                .iter()
                .map(|elem_id| (*elem_id, &self.basis_specs[*elem_id]))
                .collect())
        }
    }

    /// Retrieve a list of an `Elem`s ancestor [BasisSpec]s (All the [`BasisSpec`]s on its ancestor `Elem`s)
    pub fn ancestor_basis_specs(&self, elem_id: usize) -> Result<Vec<(usize, &Vec<BasisSpec>)>, String> {
        if elem_id >= self.mesh.elems.len() {
            Err(format!(
                "Elem {} doesn't exist; Cannot retrieve Ancestor BasisSpecs!",
                elem_id
            ))
        } else {
            let anc_elem_ids = self.mesh.ancestor_elems(elem_id, false)?;
            Ok(anc_elem_ids
                .iter()
                .map(|elem_id| (*elem_id, &self.basis_specs[*elem_id]))
                .collect())
        }
    }

    // push a new `BasisSpec` onto the list, updating its ID to match its position in its elem's list
    // return its [BSAddress] composed of its element id and index
    fn push_basis_spec(&mut self, mut bs: BasisSpec, dof_id: usize) -> BSAddress {
        let elem_id = bs.elem_id;
        let elem_idx = self.basis_specs[elem_id].len();

        bs.set_dof_and_idx(dof_id, elem_idx);
        self.basis_specs[elem_id].push(bs);

        BSAddress::new(elem_id, elem_idx)
    }

    /// Fill two system matrices using the domains it's Basis Space as the Testing Space. Return a Generalized Eigenproblem ([GEP])
    ///
    /// All pairs of overlapping Shape Functions will be integrated and stored in the matrices by their associated DoF IDs
    ///
    /// * Two [Integral]s: `AI` and `BI` must be specified. These are used to populate the A and B matrices respectively
    /// * A [ShapeFn] `SF` must also be specified. This is used to evaluate the Domains [BasisSpec]s
    pub fn galerkin_sample_gep<SF, AI, BI>(&self, num_gauss_quad: Option<usize>) -> GEP
    where
        SF: ShapeFn,
        AI: Integral,
        BI: Integral,
    {
        // construct an eigenproblem with a and b matrices
        let mut gep = GEP::new(self.dofs.len());

        // construct basis sampler
        let [i_max, j_max] = self.mesh.max_expansion_orders();
        let (mut bf_sampler, [u_weights, v_weights]): (BasisFnSampler<SF>, _) =
            BasisFnSampler::with(
                i_max as usize,
                j_max as usize,
                num_gauss_quad,
                num_gauss_quad,
                false,
            );

        // setup integration
        let a_integrator = AI::with_weights(&u_weights, &v_weights);
        let b_integrator = BI::with_weights(&u_weights, &v_weights);

        for elem in self.elems() {
            // get relevant data for this Elem
            let local_basis_specs = self.local_basis_specs(elem.id).unwrap();
            let desc_basis_specs = self.descendant_basis_specs(elem.id).unwrap();
            let bs_local = bf_sampler.sample_basis_fn(elem, None);

            let elem_materials = elem.get_materials();

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
                        .integrate(p_dir, q_dir, p_orders, q_orders, &bs_local, &bs_local)
                        .full_solution();
                    let b = b_integrator
                        .integrate(p_dir, q_dir, p_orders, q_orders, &bs_local, &bs_local)
                        .full_solution();

                    local_a_entries.push(([p_dof_id, q_dof_id], a / elem_materials.mu_rel.re));
                    local_b_entries.push(([p_dof_id, q_dof_id], b * elem_materials.eps_rel.re));
                }
            }

            gep.a.insert_group(local_a_entries);
            gep.b.insert_group(local_b_entries);

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
                        bf_sampler.sample_basis_fn(elem, Some(&self.mesh.elems[q_elem_id]));

                    let bs_q_local = bf_sampler.sample_basis_fn(&self.mesh.elems[q_elem_id], None);

                    for (q_orders, q_dir, q_dof_id) in q_elem_basis_specs
                        .iter()
                        .map(|bs_q| bs_q.integration_data())
                    {
                        let a = a_integrator
                            .integrate(p_dir, q_dir, p_orders, q_orders, &bs_p_sampled, &bs_q_local)
                            .full_solution();
                        let b = b_integrator
                            .integrate(p_dir, q_dir, p_orders, q_orders, &bs_p_sampled, &bs_q_local)
                            .full_solution();

                        desc_a_entries.push(([p_dof_id, q_dof_id], a / elem_materials.mu_rel.re));
                        desc_b_entries.push(([p_dof_id, q_dof_id], b * elem_materials.eps_rel.re));
                    }
                }
            }

            gep.a.insert_group(desc_a_entries);
            gep.b.insert_group(desc_b_entries);
        }

        gep
    }

    /// Same as `galerkin_sample_gep`, except matrices are filled in parallel using the Rayon Global ThreadPool
    pub fn galerkin_sample_gep_parallel<SF, AI, BI>(&self, num_gauss_quad: Option<usize>) -> GEP
    where
        SF: ShapeFn,
        AI: Integral,
        BI: Integral,
    {
        // construct an eigenproblem with a and b matrices
        let mut gep = GEP::new(self.dofs.len());

        // construct basis sampler
        let [i_max, j_max] = self.mesh.max_expansion_orders();
        let (bs_sampler, [u_weights, v_weights]): (ParBasisFnSampler<SF>, _) =
            ParBasisFnSampler::with(
                i_max as usize,
                j_max as usize,
                num_gauss_quad,
                num_gauss_quad,
                false,
            );
        let bf_sampler_send = Arc::new(Mutex::new(bs_sampler));

        // setup integration
        let a_integrator = AI::with_weights(&u_weights, &v_weights);
        let b_integrator = BI::with_weights(&u_weights, &v_weights);

        gep.par_extend(self.mesh.elems.par_iter().map(|elem| {
            let mut local_a = SparseMatrix::new(self.dofs.len());
            let mut local_b = SparseMatrix::new(self.dofs.len());

            let bf_sampler_elem = bf_sampler_send.clone();

            let elem_materials = elem.get_materials();

            // get relevant data for this Elem
            let bs_local = bf_sampler_elem.lock().unwrap().sample_basis_fn(elem, None);
            let local_basis_specs = self.local_basis_specs(elem.id).unwrap();
            let desc_basis_specs = self.descendant_basis_specs(elem.id).unwrap();

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
                        .integrate(p_dir, q_dir, p_orders, q_orders, &bs_local, &bs_local)
                        .full_solution();
                    let b = b_integrator
                        .integrate(p_dir, q_dir, p_orders, q_orders, &bs_local, &bs_local)
                        .full_solution();

                    local_a_entries.push(([p_dof_id, q_dof_id], a / elem_materials.mu_rel.re));
                    local_b_entries.push(([p_dof_id, q_dof_id], b * elem_materials.eps_rel.re));
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
                    let bs_p_sampled = bf_sampler_elem
                        .lock()
                        .unwrap()
                        .sample_basis_fn(elem, Some(&self.mesh.elems[q_elem_id]));
                    let bs_q_local = bf_sampler_elem
                        .lock()
                        .unwrap()
                        .sample_basis_fn(&self.mesh.elems[q_elem_id], None);

                    for (q_orders, q_dir, q_dof_id) in q_elem_basis_specs
                        .iter()
                        .map(|bs_q| bs_q.integration_data())
                    {
                        let a = a_integrator
                            .integrate(p_dir, q_dir, p_orders, q_orders, &bs_p_sampled, &bs_q_local)
                            .full_solution();
                        let b = b_integrator
                            .integrate(p_dir, q_dir, p_orders, q_orders, &bs_p_sampled, &bs_q_local)
                            .full_solution();

                        desc_a_entries.push(([p_dof_id, q_dof_id], a / elem_materials.mu_rel.re));
                        desc_b_entries.push(([p_dof_id, q_dof_id], b * elem_materials.eps_rel.re));
                    }
                }
            }

            local_a.insert_group(desc_a_entries);
            local_b.insert_group(desc_b_entries);

            [local_a, local_b]
        }));

        gep
    }

    /// Retrieve a [BasisSpec] at a particular [BSAddress]
    ///
    /// Returns an error if the designated `Elem` does not exist, or does not have that [BasisSpec]
    pub fn get_basis_spec(&self, address: BSAddress) -> Result<&BasisSpec, String> {
        if address.elem_id < self.basis_specs.len() {
            if address.elem_idx < self.basis_specs[address.elem_id].len() {
                Ok(&self.basis_specs[address.elem_id][address.elem_idx])
            } else {
                Err(format!("Could not retrieve {}", address))
            }
        } else {
            Err(format!("Could not retrieve {}", address))
        }
    }
}

struct IdTracker {
    next_id: usize,
}

impl IdTracker {
    pub fn new(start: usize) -> Self {
        Self { next_id: start }
    }

    pub fn next_id(&mut self) -> usize {
        self.next_id += 1;
        self.next_id - 1
    }

    pub fn next_two_ids(&mut self) -> [usize; 2] {
        let ids = [self.next_id, self.next_id + 1];
        self.next_id += 2;
        ids
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use mesh::h_refinement::HRef;
    use mesh::p_refinement::PRef;

    #[test]
    fn create_domain() {
        let mut mesh = Mesh::from_file("./test_input/test_mesh_a.json").unwrap();
        mesh.set_global_expansion_orders([5, 5]).unwrap();
        mesh.global_h_refinement(HRef::T).unwrap();
        mesh.h_refine_elems(vec![4, 5], HRef::T).unwrap();
        mesh.h_refine_elems(vec![6, 7], HRef::u()).unwrap();
        mesh.h_refine_elems(vec![8, 9], HRef::v()).unwrap();
        mesh.p_refine_elems(vec![10, 11, 12, 13], PRef::from(2, -1))
            .unwrap();

        let dom = Domain::from_mesh(mesh);
        dom.local_basis_specs(0).unwrap();
        dom.descendant_basis_specs(0).unwrap();
    }
}
