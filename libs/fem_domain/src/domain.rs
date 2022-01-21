mod dof;
mod mesh;

pub use dof::{BasisDir, BasisLoc, BasisSpec, DoF};
pub use mesh::*;

use dof::BSAddress;

use std::collections::BTreeMap;

/// Description of an FEM Domain with Degrees of Freedom
pub struct Domain {
    /// Geometric discretization
    pub mesh: Mesh,
    /// Degrees of Freedom (collections of BasisSpecs matched via Tangential Continuity)
    pub dofs: Vec<DoF>,
    /// Individual BasisSpecs sorted by Elem
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

    pub fn from_mesh_file(path: impl AsRef<str>) -> std::io::Result<Self> {
        Ok(Self {
            mesh: Mesh::from_file(path)?,
            dofs: Vec::new(),
            basis_specs: Vec::new(),
        })
    }

    /// Iterate over all [Elem]s in the mesh
    pub fn elems<'a>(&'a self) -> impl Iterator<Item = &'a Elem> + '_ {
        self.mesh.elems.iter()
    }

    /// Iterate over all [Edge]s in the mesh
    pub fn edges<'a>(&'a self) -> impl Iterator<Item = &'a Edge> + '_ {
        self.mesh.edges.iter()
    }

    /// Iterate over all [Node]s in the mesh
    pub fn nodes<'a>(&'a self) -> impl Iterator<Item = &'a Node> + '_ {
        self.mesh.nodes.iter()
    }

    /// Generate Degrees of Freedom over the mesh according to the Polynomial Expansion orders on each [Elem]
    pub fn gen_dofs(&mut self) {
        // prepare for fresh set of DoFs and BasisSpecs
        self.basis_specs = vec![Vec::new(); self.mesh.elems.len()];
        self.dofs.clear();
        self.mesh.set_edge_activation();
        let mut dof_id_tracker = IdTracker::new(0);

        // Generate lists of BasisSpecs associated with Elems, Edges, and Nodes, sorted by their IDs
        let [elem_bs, edge_bs, _] = self.gen_basis_specs();

        // Designate all elem-type BasisSpecs located on childless Elems as DoFs
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
                    if self.basis_specs[elem_id].len() == 0 {
                        self.basis_specs[elem_id] = Vec::with_capacity(num_expected);
                    } else {
                        self.basis_specs[elem_id].reserve(num_expected);
                    }
                }

                // iterate over each pair of BasisSpecs (once) and look for matches
                let mut active_pairs: Vec<[usize; 2]> = Vec::with_capacity(num_expected);
                for (a, bs_0) in rel_basis_specs.iter().enumerate() {
                    for (b, bs_1) in rel_basis_specs.iter().enumerate().skip(a + 1) {
                        if bs_0.matches_with(bs_1) {
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
                        BasisLoc::ELEM => elem_bs
                            .entry(elem.id)
                            .and_modify(|bs_list| bs_list.push(bs.clone()))
                            .or_insert(vec![bs]),
                        BasisLoc::EDGE(_, edge_id) => edge_bs
                            .entry(edge_id)
                            .and_modify(|bs_list| bs_list.push(bs.clone()))
                            .or_insert(vec![bs]),
                        BasisLoc::NODE(_, node_id) => node_bs
                            .entry(node_id)
                            .and_modify(|bs_list| bs_list.push(bs.clone()))
                            .or_insert(vec![bs]),
                    };
                }
            }
        }

        [elem_bs, edge_bs, node_bs]
    }

    /// Retrieve a list of [BasisSpec]s on an [Elem] by ID
    pub fn local_basis_specs<'a>(&'a self, elem_id: usize) -> Result<&'a Vec<BasisSpec>, String> {
        if elem_id >= self.mesh.elems.len() {
            Err(format!(
                "Elem {} doesn't exist; Cannot retrieve BasisSpecs!",
                elem_id
            ))
        } else {
            Ok(&self.basis_specs[elem_id])
        }
    }

    /// Retrieve a list of an [Elem]s descendant [BasisSpec]s (All the [`BasisSpec`]s on its descendant [`Elem`]s)
    pub fn descendant_basis_specs<'a>(
        &'a self,
        elem_id: usize,
    ) -> Result<Vec<(usize, &'a Vec<BasisSpec>)>, String> {
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

    /// Retrieve a list of an [Elem]s ancestor [BasisSpec]s (All the [`BasisSpec`]s on its ancestor [`Elem`]s)
    pub fn ancestor_basis_specs<'a>(
        &'a self,
        elem_id: usize,
    ) -> Result<Vec<&'a BasisSpec>, String> {
        if elem_id >= self.mesh.elems.len() {
            Err(format!(
                "Elem {} doesn't exist; Cannot retrieve Ancestor BasisSpecs!",
                elem_id
            ))
        } else {
            let anc_elem_ids = self.mesh.ancestor_elems(elem_id, false)?;
            Ok(anc_elem_ids
                .iter()
                .flat_map(|elem_id| self.basis_specs[*elem_id].iter())
                .collect())
        }
    }

    // // push a new BasisSpec onto the list, updating its ID to match its position in its elem's list
    // // return its [elem_id, elem_list_position]
    // fn push_basis_spec(&mut self, mut bs: BasisSpec, dof_id: usize) -> BSAddress {
    //     let new_id = self.basis_specs[bs.elem_id].len();
    //     let elem_id = bs.elem_id;

    //     bs.update_ids(new_id, dof_id);
    //     self.basis_specs[elem_id].push(bs);

    //     BSAddress::new(elem_id, new_id)
    // }

    // push a new BasisSpec onto the list, updating its ID to match its position in its elem's list
    // return its [elem_id, elem_list_position]
    fn push_basis_spec(&mut self, mut bs: BasisSpec, dof_id: usize) -> BSAddress {
        let bs_id = bs.id;
        let elem_id = bs.elem_id;

        // bs.update_ids(new_id, dof_id);
        bs.assign_dof_id(dof_id);
        self.basis_specs[elem_id].push(bs);

        BSAddress::new(elem_id, bs_id)
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

    #[test]
    fn create_dofs() {
        let mut dom_a = Domain::from_mesh_file("../../test_input/test_mesh_a.json").unwrap();

        dom_a.mesh.set_global_expansion_orders([5, 5]).unwrap();
        dom_a.gen_dofs();

        dom_a.mesh.global_h_refinement(HRef::T).unwrap();
        dom_a.gen_dofs();

        dom_a.mesh.h_refine_elems(vec![4, 5], HRef::T).unwrap();
        dom_a.mesh.h_refine_elems(vec![6, 7], HRef::u()).unwrap();
        dom_a.mesh.h_refine_elems(vec![8, 9], HRef::v()).unwrap();
        dom_a
            .mesh
            .p_refine_elems(vec![10, 11, 12, 13], PRef::from(2, -1))
            .unwrap();
        dom_a.gen_dofs();
    }
}
