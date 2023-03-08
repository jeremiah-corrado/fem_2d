/// Degrees of Freedom (Groups of one or more Basis Function that can take a scalar value in the FEM solution)
pub mod dof;
/// Structures used to render a solution over a Domain
pub mod fields;
/// The internal geometric structure of a Domain. This is modified by hp-refinements.
pub mod mesh;

use dof::{
    basis_spec::{BSAddress, BasisDir, BasisLoc, BasisSpec},
    DoF,
};
use mesh::*;
use smallvec::smallvec;
use std::collections::BTreeMap;
use std::fmt;

/// The Continuity Condition to be enforced by the Domain. Only H(Curl) is currently supported!!!
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ContinuityCondition {
    HCurl,
    HDiv,
    Discontinuous,
}

impl fmt::Display for ContinuityCondition {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::HCurl => write!(f, "H(Curl)"),
            Self::HDiv => write!(f, "H(Div)"),
            Self::Discontinuous => write!(f, "Discontinuous"),
        }
    }
}

/// High Level Description of an FEM Domain
///
/// This struct contains:
/// * The [Mesh] - A discretization of the Domain of interest
/// * A list of [DoF]s (Degrees of Freedom) - Sets of Basis Functions that can take a scalar value in the FEM solution
/// * A list of [BasisSpec]s - The specification of each Basis Function in the Domain
///
pub struct Domain {
    pub mesh: Mesh,
    /// Degrees of Freedom (collections of BasisSpecs matched via Tangential Continuity)
    pub dofs: Vec<DoF>,
    /// Individual BasisSpecs sorted by their associated Elem
    pub basis_specs: Vec<Vec<BasisSpec>>,
    /// The continuity condition enforced over the Basis Space
    pub cc: ContinuityCondition,
}

impl Domain {
    /// Create a Domain around an empty mesh
    pub fn blank(cc: ContinuityCondition) -> Self {
        Self {
            mesh: Mesh::blank(),
            dofs: Vec::new(),
            basis_specs: Vec::new(),
            cc,
        }
    }

    /// Create a Domain from a unit Mesh
    pub fn unit(cc: ContinuityCondition) -> Self {
        Self::from_mesh(Mesh::unit(), cc)
    }

    /// Create a Domain from a Mesh
    pub fn from_mesh(mut mesh: Mesh, cc: ContinuityCondition) -> Self {
        // prepare for basis function matching
        mesh.set_edge_activation();
        // mesh.set_node_activation(); (TODO: this is not needed until node-type basis functions are implemented)

        // create dof and basis_spec collections
        let mut dof_id_tracker = IdTracker::new(0);
        let mut basis_specs = vec![Vec::new(); mesh.elems.len()];
        let mut dofs = Vec::new();

        // Generate lists of BasisSpecs associated with Elems, Edges, and (Nodes), sorted by their IDs
        let [elem_bs, edge_bs, _] = Self::gen_basis_specs(&mesh, cc);

        // Designate all elem-type BasisSpecs located on shell Elems as DoFs
        for (elem_id, mut elem_bs_list) in elem_bs {
            if !mesh.elems[elem_id].has_children() {
                basis_specs[elem_id] = Vec::with_capacity(elem_bs_list.len());

                for elem_bs in elem_bs_list
                    .drain(0..)
                    .filter(|bs| bs.dir == BasisDir::U || bs.dir == BasisDir::V)
                {
                    let dof_id = dof_id_tracker.next_id();
                    let address = Self::push_basis_spec(&mut basis_specs, elem_bs, dof_id);
                    dofs.push(DoF::new(dof_id, smallvec![address]));
                }
            }
        }

        // Create DoFs from pairs of matched BasisSpecs on the active Elems associated with each Edge
        for (edge_id, mut edge_bs_list) in edge_bs {
            if let Some(active_elem_ids) = mesh.edges[edge_id].active_elem_pair() {
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
                    if basis_specs[elem_id].is_empty() {
                        basis_specs[elem_id] = Vec::with_capacity(num_expected);
                    } else {
                        basis_specs[elem_id].reserve(num_expected);
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
                            Self::push_basis_spec(
                                &mut basis_specs,
                                rel_basis_specs[*rel_idx].clone(),
                                dof_id,
                            )
                        })
                        .collect();

                    dofs.push(DoF::new(dof_id, addresses));
                }
            }
        }

        // TODO: implement node-type BasisSpec Matching!

        Self {
            mesh,
            dofs,
            basis_specs,
            cc,
        }
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

    /// Retrieve a [BasisSpec] at a particular [BSAddress]
    ///
    /// # Returns
    /// * A reference to the [BasisSpec] if it exists
    /// * An `Err` if the `elem` does not exist
    /// * An `Err` if the `elem` does not have a [BasisSpec] at the given `address`
    ///
    pub fn get_basis_spec(&self, bs_address: BSAddress) -> Result<&BasisSpec, String> {
        if bs_address.elem_id > self.mesh.elems.len() {
            Err(format!(
                "Elem {} does not exist; cannot retrieve BasisSpec!",
                bs_address.elem_id
            ))
        } else if self.basis_specs[bs_address.elem_id].len() <= bs_address.elem_idx {
            Err(format!(
                "Elem {} only has {} BasisSpecs; cannot retrive BasisSpec({})",
                bs_address.elem_id,
                self.basis_specs[bs_address.elem_id].len(),
                bs_address.elem_idx
            ))
        } else {
            Ok(&self.basis_specs[bs_address.elem_id][bs_address.elem_idx])
        }
    }

    fn gen_basis_specs(
        mesh: &Mesh,
        cc: ContinuityCondition,
    ) -> [BTreeMap<usize, Vec<BasisSpec>>; 3] {
        let mut elem_bs: BTreeMap<usize, Vec<BasisSpec>> = BTreeMap::new();
        let mut edge_bs: BTreeMap<usize, Vec<BasisSpec>> = BTreeMap::new();
        let mut node_bs: BTreeMap<usize, Vec<BasisSpec>> = BTreeMap::new();

        let mut bs_id_tracker = IdTracker::new(0);

        for elem in mesh.elems.iter() {
            for dir in [BasisDir::U, BasisDir::V, BasisDir::W] {
                for poly_ij in elem.poly_orders.permutations(dir) {
                    let bs = BasisSpec::new(bs_id_tracker.next_id(), poly_ij, dir, elem, cc);

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
    ///
    /// # Example
    /// ```
    /// use fem_2d::prelude::*;
    /// let mut mesh = Mesh::unit();
    /// mesh.set_global_expansion_orders(Orders::new(2, 2));
    ///
    /// let dom = Domain::from_mesh(mesh, ContinuityCondition::HCurl);
    ///
    /// // get Elem 0's basis specs
    /// let basis_specs = dom.local_basis_specs(0).unwrap();
    ///
    /// // 2 elem-type basis specs in each direction
    /// assert_eq!(basis_specs.len(), 4);
    /// ```
    pub fn local_basis_specs(&self, elem_id: usize) -> Result<&Vec<BasisSpec>, MeshAccessError> {
        if elem_id >= self.mesh.elems.len() {
            Err(MeshAccessError::ElemDoesNotExist(elem_id))
        } else {
            Ok(&self.basis_specs[elem_id])
        }
    }

    /// Retrieve a list of an `Elem`s descendant [BasisSpec]s (All the [`BasisSpec`]s on its descendant `Elem`s)
    ///
    /// # Returns
    /// * A vector of tuples in the form: `(desc_elem_id: usize, desc_elem_basis_specs: Vec<BasisSpec>)`. Where each tuple corresponds to one descendant `Elem`
    /// * An error if the given `elem_id` does not exist
    ///
    /// The [BasisSpec]s on `elem_id` itself are not included in the returned vector.
    ///
    /// # Example
    /// ```
    /// use fem_2d::prelude::*;
    ///
    /// let mut mesh = Mesh::unit();
    /// mesh.set_global_expansion_orders(Orders::new(2, 2));
    /// mesh.global_h_refinement(HRef::T);
    ///
    /// let dom = Domain::from_mesh(mesh, ContinuityCondition::HCurl);
    ///
    /// // get the basis specs from all of Elem 0's descendants
    /// let basis_specs = dom.descendant_basis_specs(0).unwrap();
    ///
    /// // Elem 0 has 4 descendant `Elem`s due to the T-Type h-Refinement
    /// assert_eq!(basis_specs.len(), 4);
    ///
    /// // there are 4 element-type and 4 edge-type basis specs in each direction on each descendant Elem
    /// assert_eq!(basis_specs[0].1.len(), 8);
    /// assert_eq!(basis_specs[1].1.len(), 8);
    /// assert_eq!(basis_specs[2].1.len(), 8);
    /// assert_eq!(basis_specs[3].1.len(), 8);
    /// ```
    pub fn descendant_basis_specs(
        &self,
        elem_id: usize,
    ) -> Result<Vec<(usize, &Vec<BasisSpec>)>, MeshAccessError> {
        if elem_id >= self.mesh.elems.len() {
            Err(MeshAccessError::ElemDoesNotExist(elem_id))
        } else {
            let desc_elem_ids = self.mesh.descendant_elems(elem_id, false)?;
            Ok(desc_elem_ids
                .iter()
                .map(|elem_id| (*elem_id, &self.basis_specs[*elem_id]))
                .collect())
        }
    }

    /// Retrieve a list of an `Elem`s ancestor [BasisSpec]s (All the [`BasisSpec`]s on its ancestor `Elem`s)
    ///
    /// # Returns
    /// * A vector of tuples in the form: `(anc_elem_id: usize, anc_elem_basis_specs: Vec<BasisSpec>)`. Where each tuple corresponds to one ancestor `Elem`
    /// * An error if the given `elem_id` does not exist
    ///
    /// The [BasisSpec]s on `elem_id` itself are not included in the returned vector.
    ///
    /// # Example
    /// ```
    /// use fem_2d::prelude::*;
    ///
    /// let mut mesh = Mesh::unit();
    /// mesh.set_global_expansion_orders(Orders::new(2, 2));
    /// mesh.global_h_refinement(HRef::T);
    /// mesh.h_refine_elems(vec![1], HRef::T).unwrap();
    ///
    /// let dom = Domain::from_mesh(mesh, ContinuityCondition::HCurl);
    ///
    /// // get the basis specs from elem 5's ancestors
    /// let basis_specs = dom.ancestor_basis_specs(5).unwrap();
    ///
    /// // there are two ancestor Elems: 0 and 1
    /// assert_eq!(basis_specs.len(), 2);
    ///
    /// // elem 0 is inactive and has no basis specs
    /// let elem_0_bs = basis_specs.iter().find(|(elem_id, _)| *elem_id == 0).unwrap().1;
    /// assert!(elem_0_bs.is_empty());
    ///
    /// // elem 1 has 2 edge-type basis specs on two of its edges
    /// let elem_1_bs = basis_specs.iter().find(|(elem_id, _)| *elem_id == 1).unwrap().1;
    /// assert_eq!(elem_1_bs.len(), 4);
    /// ```
    pub fn ancestor_basis_specs(
        &self,
        elem_id: usize,
    ) -> Result<Vec<(usize, &Vec<BasisSpec>)>, MeshAccessError> {
        if elem_id >= self.mesh.elems.len() {
            Err(MeshAccessError::ElemDoesNotExist(elem_id))
        } else {
            let anc_elem_ids = self.mesh.ancestor_elems(elem_id, false)?;
            Ok(anc_elem_ids
                .iter()
                .map(|elem_id| (*elem_id, &self.basis_specs[*elem_id]))
                .collect())
        }
    }

    // Push a new `BasisSpec` onto the list, updating its ID to match its position in its elem's list
    // return its [BSAddress] composed of its element id and index
    fn push_basis_spec(
        basis_specs: &mut Vec<Vec<BasisSpec>>,
        mut bs: BasisSpec,
        dof_id: usize,
    ) -> BSAddress {
        let elem_id = bs.elem_id;
        let elem_idx = basis_specs[elem_id].len();

        bs.set_dof_and_idx(dof_id, elem_idx);
        basis_specs[elem_id].push(bs);

        BSAddress::new(elem_id, elem_idx)
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
    use mesh::p_refinement::{PRef, Orders};

    #[test]
    fn create_domain() {
        let mut mesh = Mesh::from_file("./test_input/test_mesh_a.json").unwrap();
        mesh.set_global_expansion_orders(Orders::try_new(5, 5).unwrap());
        mesh.global_h_refinement(HRef::T);
        mesh.h_refine_elems(vec![4, 5], HRef::T).unwrap();
        mesh.h_refine_elems(vec![6, 7], HRef::u()).unwrap();
        mesh.h_refine_elems(vec![8, 9], HRef::v()).unwrap();
        mesh.p_refine_elems(vec![10, 11, 12, 13], PRef::from(2, -1))
            .unwrap();

        let dom = Domain::from_mesh(mesh, ContinuityCondition::HCurl);
        dom.local_basis_specs(0).unwrap();
        dom.descendant_basis_specs(0).unwrap();
    }
}
