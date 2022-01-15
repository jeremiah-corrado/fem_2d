mod dof;
mod mesh;

pub use dof::{BasisDir, DoF, BasisLoc, BasisSpec};
pub use mesh::*;

use std::collections::BTreeMap;

/// Description of an FEM Domain with Degrees of Freedom
pub struct Domain {
    /// Geometric discretization
    pub mesh: Mesh,
    /// Degrees of Freedom (collections of BasisSpecs matched via Tangential Continuity)
    pub dofs: Vec<DoF>,
    /// Individual BasisSpecs sorted by Elem
    basis_specs: Vec<Vec<BasisSpec>>,
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
        self.basis_specs.clear();
        self.dofs.clear();

        let [elem_bs, edge_bs, _] = self.gen_basis_specs();

        let mut dof_id_tracker = IdTracker::new(0);

        // Designate all elem-type BasisSpecs located on childless Elems as DoFs
        for (elem_id, elem_bs_list) in elem_bs {
            if !self.mesh.elems[elem_id].has_children() {
                for elem_bs in elem_bs_list.iter() {
                    self.dofs
                        .push(DoF::new(dof_id_tracker.next_id(), &[elem_bs]));
                }
                self.basis_specs[elem_id] = elem_bs_list;
            }
        }

        // Create Dofs from pairs of matched BasisSpecs on the active Elems associated with each Edge
        for (edge_id, mut edge_bs_list) in edge_bs {
            if let Some(active_elem_ids) = self.mesh.edges[edge_id].active_elem_pair() {
                let relevant_basis_specs: Vec<BasisSpec> = edge_bs_list
                    .drain(0..)
                    .filter(|bs| active_elem_ids.contains(&bs.elem_id))
                    .collect();

                let mut active_pairs: Vec<[usize; 2]> =
                    Vec::with_capacity(relevant_basis_specs.len() / 2);
                self.basis_specs[active_elem_ids[0]] =
                    Vec::with_capacity(relevant_basis_specs.len() / 2);
                self.basis_specs[active_elem_ids[1]] =
                    Vec::with_capacity(relevant_basis_specs.len() / 2);

                // iterate over each pair of BasisSpecs (once) and look for matches
                for (a, bs_0) in relevant_basis_specs.iter().enumerate() {
                    for (b, bs_1) in relevant_basis_specs.iter().enumerate().skip(a + 1) {
                        if bs_0.matches_with(bs_1) {
                            active_pairs.push([a, b]);
                        }
                    }
                }

                // create new DoFs from the matched pairs
                // and add the individual BasisSpecs to their Elem's list
                for ap in active_pairs {
                    self.dofs.push(DoF::new(
                        dof_id_tracker.next_id(),
                        &[&relevant_basis_specs[ap[0]], &relevant_basis_specs[ap[1]]],
                    ));

                    self.basis_specs[relevant_basis_specs[ap[0]].elem_id]
                        .push(relevant_basis_specs[ap[0]].clone());
                    self.basis_specs[relevant_basis_specs[ap[1]].elem_id]
                        .push(relevant_basis_specs[ap[1]].clone());
                }
            }
        }

        // ignoring Node-Type DoFs for now. Needs to be implemented!
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
    ) -> Result<Vec<&'a BasisSpec>, String> {
        if elem_id >= self.mesh.elems.len() {
            Err(format!(
                "Elem {} doesn't exist; Cannot retrieve Descendant BasisSpecs!",
                elem_id
            ))
        } else {
            let desc_elem_ids = self.mesh.descendant_elems(elem_id, false)?;
            Ok(desc_elem_ids
                .iter()
                .flat_map(|elem_id| self.basis_specs[*elem_id].iter())
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
