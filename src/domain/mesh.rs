/// A line between two Nodes
pub mod edge;
/// A Finite Element in Parametric Space
pub mod elem;
/// A Finite Element in Real Space
pub mod element;
/// Structures and Functions to facilitate RBS based anisotropic h-refinement
pub mod h_refinement;
/// A Point in Real Space
pub mod node;
/// Structures and Functions to facilitate anisotropic p-refinement
pub mod p_refinement;
/// Structures to describe the 2D real and parametric spaces defining a Mesh
pub mod space;

use edge::Edge;
use elem::{Elem, ElemUninit};
use element::{Element, Materials};
use h_refinement::{HRef, HRefError};
use node::Node;
use p_refinement::{PRef, PRefError};
use space::{ParaDir, Point};

use super::IdTracker;

use json::{object, JsonValue};
use smallvec::SmallVec;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::fs::{read_to_string, File};
use std::io::BufWriter;
use std::sync::Arc;

/// Minimum Edge length in parametric space. h-Refinements will fail after edges are smaller than this value.
pub const MIN_EDGE_LENGTH: f64 = 3.0518e-5; // 15ish refinement layers with unit sized cells

/// The expected "h-Refinement" depth. This determines the stack allocation size of some `SmallVec`s related to h-Refinement
pub const EXPECTED_NUM_H_REFINEMENTS: usize = 8;

/// Maximum Polynomial expansion. p-Refinements will fail when Elem's expansion orders exceed this value.
pub const MAX_POLYNOMIAL_ORDER: u8 = 20;

/// Information used to Define the geometric structure and refinement state of a Domain.
#[derive(Debug, Clone)]
pub struct Mesh {
    pub elements: Vec<Arc<Element>>,
    pub elems: Vec<Elem>,
    pub nodes: Vec<Node>,
    pub edges: Vec<Edge>,
}

impl Mesh {
    /// Construct a completely empty Mesh
    pub fn blank() -> Self {
        Self {
            elements: Vec::new(),
            elems: Vec::new(),
            nodes: Vec::new(),
            edges: Vec::new(),
        }
    }

    /// Construct a Mesh from a JSON file with the following format
    ///
    /// The first "Element" and "Node" describe the meaning of each variable
    ///
    /// The following entries in each array describe this two element mesh:
    /// ```text
    ///     3               4               5
    /// 0.5 *---------------*---------------*
    ///     |               |               |
    ///     |      air      |    teflon     |
    ///     |               |               |
    /// 0.0 *---------------*---------------*
    ///  y  0               1               2
    ///  x 0.0             1.0             2.0
    /// ```
    ///
    /// mesh.json
    /// ```JSON
    /// {
    ///     "Elements": [
    ///         {
    ///             "materials": [eps_rel_re, eps_rel_im, mu_rel_re, mu_rel_im],
    ///             "node_ids": [node_0_id, node_1_id, node_2_id, node_3_id],
    ///         },
    ///         {
    ///             "materials": [1.0, 0.0, 1.0, 0.0],
    ///             "node_ids": [1, 2, 4, 5],
    ///         },
    ///         {
    ///             "materials": [1.2, 0.0, 0.9999, 0.0],
    ///             "node_ids": [2, 3, 5, 6],
    ///         }
    ///     ],
    ///     "Nodes": [
    ///         [x_coordinate, y_coordinate],
    ///         [0.0, 0.0],
    ///         [1.0, 0.0],
    ///         [2.0, 0.0],
    ///         [0.0, 0.5],
    ///         [1.0, 0.5],
    ///         [2.0, 0.5],
    ///     ]
    /// }
    /// ```
    pub fn from_file(path: impl AsRef<str>) -> std::io::Result<Self> {
        // parse mesh file as JSON
        let mesh_file_contents = read_to_string(path.as_ref())?;
        let mesh_file_json =
            json::parse(&mesh_file_contents).expect("Unable to parse Mesh File as JSON!");

        // extract element material parameters and node_id sets (panicking if JSON format is not correct)
        let (mut element_materials, mut element_node_ids) =
            parse_element_information(&mesh_file_json);

        // extract node locations (panicking if JSON format is not correct)
        let points = parse_node_information(&mesh_file_json);

        // build a vector of elements with the specified nodes and material properties
        let elements: Vec<Arc<Element>> = element_materials
            .drain(0..)
            .zip(element_node_ids.iter())
            .enumerate()
            .map(|(element_id, (materials, node_ids))| {
                Arc::new(Element::new(
                    element_id,
                    [
                        points[node_ids[0]],
                        points[node_ids[1]],
                        points[node_ids[2]],
                        points[node_ids[3]],
                    ],
                    materials,
                ))
            })
            .collect();

        // count the number of times each point/node is referenced by an element
        let mut node_connection_counts = vec![0; points.len()];
        for node_ids in element_node_ids.iter() {
            for node_id in node_ids.iter() {
                node_connection_counts[*node_id] += 1;
            }
        }
        // mark nodes with fewer than 4 references as boundary nodes
        let boundary_nodes: Vec<bool> = node_connection_counts
            .iter()
            .map(|count| {
                assert!(
                    *count <= 4,
                    "Nodes can only be shared by a maximum of 4 Elements; Cannot construct Mesh from file!"
                );
                *count < 4
            })
            .collect();

        // build a vector of nodes from the above information
        let mut nodes: Vec<Node> = points
            .iter()
            .enumerate()
            .map(|(node_id, point)| Node::new(node_id, *point, boundary_nodes[node_id]))
            .collect();

        // build a map which describes all the edges and which elements/elems they are adjacent to on each side
        // {[node_id_0, node_id_1] => [LB element_id, TR element_id]}
        let mut edge_node_pairs: BTreeMap<[usize; 2], [Option<usize>; 2]> = BTreeMap::new();
        for (element_id, element_node_ids) in element_node_ids.iter().enumerate() {
            for (edge_index_pair, element_side_index) in EDGE_IDX_DEFS {
                edge_node_pairs
                    .entry([
                        element_node_ids[edge_index_pair[0]],
                        element_node_ids[edge_index_pair[1]],
                    ])
                    .and_modify(|edges_element_ids| {
                        assert!(
                            edges_element_ids[element_side_index].is_none(),
                            "'Edge': {:?}'s, Elem ({}); has already been set to {:?}; cannot set to {}",
                            edge_index_pair,
                            element_side_index,
                            edges_element_ids[element_side_index],
                            element_id,
                        );
                        edges_element_ids[element_side_index] = Some(element_id);
                    })
                    .or_insert(match element_side_index {
                        0 => [Some(element_id), None],
                        1 => [None, Some(element_id)],
                        _ => unreachable!(),
                    });
            }
        }

        // decide which edges are on the boundary
        //  - first, if there is only one element attached to the edge
        let boundary_edges: Vec<bool> = edge_node_pairs
            .values()
            .enumerate()
            .map(
                |(edge_id, adj_elem_ids)| match adj_elem_ids.iter().flatten().count() {
                    0 => panic!(
                        "Edge {} has no adjacent Elements; cannot construct Mesh from file!",
                        edge_id
                    ),
                    1 => true,
                    2 => false,
                    _ => unreachable!(),
                },
            )
            .collect();

        // build a vector of edges defined by the above sets of two nodes. Mark them as boundary edges if they have only one adjacent element
        let mut edges: Vec<Edge> = edge_node_pairs
            .keys()
            .enumerate()
            .map(|(edge_id, node_ids)| {
                Edge::new(
                    edge_id,
                    [&nodes[node_ids[0]], &nodes[node_ids[1]]],
                    boundary_edges[edge_id],
                )
            })
            .collect();

        // invert 'edge_node_pairs' S.T. we have a list of edges associated with each element in the correct order
        let mut elem_edges: Vec<[Option<usize>; 4]> = vec![[None; 4]; elements.len()];
        for (edge_id, adj_element_ids) in edge_node_pairs.values().enumerate() {
            for (elem_side_idx, elem_id) in adj_element_ids
                .iter()
                .enumerate()
                .filter(|(_, adj_elem_id)| adj_elem_id.is_some())
                .map(|(side_idx, ajd_elem_id)| (side_idx, ajd_elem_id.unwrap()))
            {
                let edge_idx = match (elem_side_idx, edges[edge_id].dir) {
                    (0, ParaDir::U) => 1,
                    (1, ParaDir::U) => 0,
                    (0, ParaDir::V) => 3,
                    (1, ParaDir::V) => 2,
                    _ => unreachable!(),
                };

                assert!(
                    elem_edges[elem_id][edge_idx].is_none(),
                    "Edge ({}) of Elem {} has already been set to {:?}; cannot set to {}!",
                    edge_idx,
                    elem_id,
                    elem_edges[elem_id][edge_idx],
                    edge_id,
                );

                elem_edges[elem_id][edge_idx] = Some(edge_id);
            }
        }

        // create a vector of Elems from the above information and connect them to the relevant Edges
        let elems: Vec<Elem> = element_node_ids
            .drain(0..)
            .enumerate()
            .map(|(elem_id, node_ids)| {
                let edge_ids: [usize; 4] = elem_edges[elem_id]
                    .iter()
                    .flatten()
                    .copied()
                    .collect::<Vec<usize>>()
                    .try_into()
                    .unwrap();

                let elem = Elem::new(elem_id, node_ids, edge_ids, elements[elem_id].clone());

                for edge_id in elem.edges.iter() {
                    edges[*edge_id].connect_elem(&elem)
                }

                for node_id in elem.nodes.iter() {
                    nodes[*node_id].connect_elem(&elem);
                }

                elem
            })
            .collect();

        Ok(Self {
            elements,
            elems,
            nodes,
            edges,
        })
    }

    /// Print the mesh to a JSON file specified by path.
    #[cfg(feature = "json_export")]
    pub fn export_to_json(&self, path: impl AsRef<str>) -> std::io::Result<()> {
        let f = File::create(path.as_ref())?;
        let mut w = BufWriter::new(&f);

        let mesh_object = object! {
            "Elements": JsonValue::from(self.elements.iter().map(|element| element.to_json()).collect::<Vec<_>>()),
            "Elems": JsonValue::from(self.elems.iter().map(|elem| elem.to_json()).collect::<Vec<_>>()),
            "Nodes": JsonValue::from(self.nodes.iter().map(|node| node.to_json()).collect::<Vec<_>>()),
            "Edges": JsonValue::from(self.edges.iter().map(|edge| edge.to_json()).collect::<Vec<_>>()),
        };

        mesh_object.write_pretty(&mut w, 4)?;

        Ok(())
    }

    // ----------------------------------------------------------------------------------------------------
    // General Data Retrieval
    // ----------------------------------------------------------------------------------------------------

    /// Get the four [Point]s composing an [`Elem`]
    pub fn elem_points(&self, elem_id: usize) -> [&Point; 4] {
        assert!(elem_id < self.elems.len());
        self.elems[elem_id]
            .nodes
            .map(|node_id| &self.nodes[node_id].coords)
    }

    /// Get the "smallest" and "largest" [Point]s composing an [`Elem`]
    pub fn elem_diag_points(&self, elem_id: usize) -> [&Point; 2] {
        assert!(elem_id < self.elems.len());
        [
            &self.nodes[self.elems[elem_id].nodes[0]].coords,
            &self.nodes[self.elems[elem_id].nodes[3]].coords,
        ]
    }

    /// Get the two [Point]s composing an [`Edge`]
    pub fn edge_points(&self, edge_id: usize) -> [&Point; 2] {
        assert!(edge_id < self.edges.len());
        [
            &self.nodes[self.edges[edge_id].nodes[0]].coords,
            &self.nodes[self.edges[edge_id].nodes[1]].coords,
        ]
    }

    /// Get a list of an [`Elem`]s descendant's IDs
    pub fn descendant_elems(
        &self,
        elem_id: usize,
        include_starting_elem: bool,
    ) -> Result<Vec<usize>, String> {
        if elem_id >= self.elems.len() {
            Err(format!(
                "Elem {} doesn't exist; Cannot retrieve descendant Elems!",
                elem_id
            ))
        } else {
            let mut descendants = Vec::new();
            self.rec_descendant_elems(elem_id, include_starting_elem, &mut descendants);
            Ok(descendants)
        }
    }

    fn rec_descendant_elems(&self, elem_id: usize, include: bool, desc: &mut Vec<usize>) {
        if include {
            desc.push(elem_id);
        }
        if let Some(child_elem_ids) = self.elems[elem_id].child_ids() {
            for cei in child_elem_ids {
                self.rec_descendant_elems(cei, true, desc);
            }
        }
    }

    /// Get a list of an [`Elem`]s ancestors's IDs
    pub fn ancestor_elems(
        &self,
        elem_id: usize,
        include_starting_elem: bool,
    ) -> Result<Vec<usize>, String> {
        if elem_id >= self.elems.len() {
            Err(format!(
                "Elem {} doesn't exist; Cannot retrieve ancestor Elems!",
                elem_id
            ))
        } else {
            let mut ancestors = Vec::new();
            self.rec_ancestor_elems(elem_id, include_starting_elem, &mut ancestors);
            Ok(ancestors)
        }
    }

    fn rec_ancestor_elems(&self, elem_id: usize, include: bool, anc: &mut Vec<usize>) {
        if include {
            anc.push(elem_id);
        }
        if let Some(parent_elem_id) = self.elems[elem_id].parent_id() {
            self.rec_ancestor_elems(parent_elem_id, true, anc);
        }
    }

    /// Get a list of an [`Edge`]s descendant's IDs
    pub fn descendant_edges(
        &self,
        edge_id: usize,
        include_starting_edge: bool,
    ) -> Result<Vec<usize>, String> {
        if edge_id > self.edges.len() {
            Err(format!(
                "Edge {} does not exist; Cannot retrieve descendant Edges!",
                edge_id
            ))
        } else {
            let mut descendants = Vec::new();
            self.rec_descendant_edges(edge_id, include_starting_edge, &mut descendants);
            Ok(descendants)
        }
    }

    fn rec_descendant_edges(&self, edge_id: usize, include: bool, desc: &mut Vec<usize>) {
        if include {
            desc.push(edge_id);
        }
        if let Some(child_edge_ids) = self.edges[edge_id].child_ids() {
            for cei in child_edge_ids {
                self.rec_descendant_edges(cei, true, desc);
            }
        }
    }

    /// Maximum polynomial expansion orders represented among all `Elem`s in the `Mesh`
    pub fn max_expansion_orders(&self) -> [u8; 2] {
        self.elems
            .iter()
            .fold([0; 2], |acc, elem| elem.poly_orders.max_with(acc))
    }

    /// Determine if this Elem can be h-refined
    /// * returns false if the Elem already has children
    /// * returns false if any of the Elem's Edges are shorter than [MIN_EDGE_LENGTH]
    /// * returns an `Err` if the Mesh doesn't have `elem_id`
    pub fn elem_is_h_refineable(&self, elem_id: usize) -> Result<bool, HRefError> {
        if elem_id >= self.elems.len() {
            Err(HRefError::ElemDoesntExist(elem_id))
        } else {
            let elem = &self.elems[elem_id];
            Ok(!elem.has_children()
                && elem
                    .edges
                    .iter()
                    .all(|edge_id| self.edges[*edge_id].length > MIN_EDGE_LENGTH))
        }
    }

    /// Determine if this Elem can be p-refined (in the positive direction)
    /// * returns false if the Elem's expansion orders have exceeded [MAX_POLYNOMIAL_ORDER] in either direction
    /// * returns an `Err` if the Mesh doesn't have `elem_id`
    pub fn elem_id_p_refineable(&self, elem_id: usize) -> Result<bool, PRefError> {
        if elem_id >= self.elems.len() {
            Err(PRefError::ElemDoesntExist(elem_id))
        } else {
            Ok(self.elems[elem_id].poly_orders.ni < MAX_POLYNOMIAL_ORDER
                && self.elems[elem_id].poly_orders.nj < MAX_POLYNOMIAL_ORDER)
        }
    }

    // ----------------------------------------------------------------------------------------------------
    // h-refinement methods
    // ----------------------------------------------------------------------------------------------------

    /// Apply an [HRef] to all [Elem]s in the Mesh that are eligible for h-refinement
    pub fn global_h_refinement(&mut self, refinement: HRef) -> Result<(), HRefError> {
        self.execute_h_refinements(
            self.elems
                .iter()
                .filter(|elem| self.elem_is_h_refineable(elem.id).unwrap())
                .map(|shell_elem| (shell_elem.id, refinement))
                .collect(),
        )
    }

    /// Apply an [HRef] to a list of [Elem]s by their ID
    pub fn h_refine_elems(&mut self, elems: Vec<usize>, refinement: HRef) -> Result<(), HRefError> {
        self.execute_h_refinements(elems.iter().map(|elem_id| (*elem_id, refinement)).collect())
    }

    /// h-refine [Elem]s according to an external filter function
    pub fn h_refine_with_filter<F>(&mut self, filt: F) -> Result<(), HRefError>
    where
        F: Fn(&Elem) -> Option<HRef>,
    {
        self.execute_h_refinements(
            self.elems
                .iter()
                .filter(|elem| self.elem_is_h_refineable(elem.id).unwrap())
                .map(|elem| (elem.id, filt(elem)))
                .filter(|(_, refinement)| refinement.is_some())
                .map(|(id, refinement)| (id, refinement.unwrap()))
                .collect(),
        )
    }

    /// Execute a series of [HRef]s on [Elem]s specified by their id
    pub fn execute_h_refinements(
        &mut self,
        refinements: Vec<(usize, HRef)>,
    ) -> Result<(), HRefError> {
        let mut refinements_map: BTreeMap<usize, HRef> = BTreeMap::new();
        for (elem_id, h_ref) in refinements {
            if elem_id >= self.elems.len() {
                return Err(HRefError::ElemDoesntExist(elem_id));
            }
            if refinements_map.insert(elem_id, h_ref).is_some() {
                return Err(HRefError::DoubleRefinement(elem_id));
            }
        }

        let mut refinement_extensions: Vec<(usize, HRef)> = Vec::new();
        let mut elem_id_tracker = self.elems.len();
        let mut node_id_tracker = IdTracker::new(self.nodes.len());
        let mut edge_id_tracker = IdTracker::new(self.edges.len());

        for (elem_id, refinement) in refinements_map {
            let new_uninitialized_elems =
                self.elems[elem_id].h_refine(refinement, &mut elem_id_tracker)?;

            let mut new_elems = match refinement {
                HRef::T => self.execute_t_refinement(
                    new_uninitialized_elems,
                    elem_id,
                    &mut node_id_tracker,
                    &mut edge_id_tracker,
                )?,
                HRef::U(extension) => {
                    let new_u_elems = self.execute_u_refinement(
                        new_uninitialized_elems,
                        elem_id,
                        &mut node_id_tracker,
                        &mut edge_id_tracker,
                    )?;
                    if let Some(child_idx) = extension {
                        refinement_extensions.push((new_u_elems[child_idx].id, HRef::V(None)))
                    }
                    new_u_elems
                }
                HRef::V(extension) => {
                    let new_v_elems = self.execute_v_refinement(
                        new_uninitialized_elems,
                        elem_id,
                        &mut node_id_tracker,
                        &mut edge_id_tracker,
                    )?;
                    if let Some(child_idx) = extension {
                        refinement_extensions.push((new_v_elems[child_idx].id, HRef::U(None)))
                    }
                    new_v_elems
                }
            };

            self.elems.extend(new_elems.drain(0..));
        }

        if !refinement_extensions.is_empty() {
            self.execute_h_refinements(refinement_extensions)?;
        }

        Ok(())
    }

    fn execute_t_refinement(
        &mut self,
        mut new_elems: Vec<ElemUninit>,
        parent_elem_id: usize,
        node_id_tracker: &mut IdTracker,
        edge_id_tracker: &mut IdTracker,
    ) -> Result<Vec<Elem>, HRefError> {
        assert_eq!(new_elems.len(), 4);

        // create a new node in the center of the parent Elem
        let parent_elem_points = self.elem_points(parent_elem_id);
        let center_node_id = node_id_tracker.next_id();
        let center_point = Point::between(parent_elem_points[0], parent_elem_points[3]);

        assert_eq!(center_node_id, self.nodes.len());
        self.nodes
            .push(Node::new(center_node_id, center_point, false));

        // connect the child Elems to the center node and the parents nodes
        for (elem_idx, elem_uninit) in new_elems.iter_mut().enumerate() {
            elem_uninit.set_node(3 - elem_idx, center_node_id);
            elem_uninit.set_node(elem_idx, self.elems[parent_elem_id].nodes[elem_idx]);
        }

        // Iterate over each parent Edge and refine it if necessary
        //      connect the child Elems to the edges descendant edges and node
        // Create a new Edge between the center_node and the new descendant Node
        //      connect the child Elems to the new Edge
        for (edge_index, adj_child_elem_indices, shared_node_indices, internal_edge_idx) in [
            (0, [0, 1], [1, 0], [3, 2]),
            (1, [2, 3], [3, 2], [3, 2]),
            (2, [0, 2], [2, 0], [1, 0]),
            (3, [1, 3], [3, 1], [1, 0]),
        ] {
            // get ids of child edges and node. Refine the parent edge if it hasn't been refined already
            let (child_edge_ids, shared_node_id) = self.h_refine_edge_if_needed(
                self.elems[parent_elem_id].edges[edge_index],
                node_id_tracker,
                edge_id_tracker,
            )?;

            // connect the uninitialized elements to the child Edges and Node
            new_elems[adj_child_elem_indices[0]].set_edge(edge_index, child_edge_ids[0]);
            new_elems[adj_child_elem_indices[1]].set_edge(edge_index, child_edge_ids[1]);

            new_elems[adj_child_elem_indices[0]].set_node(shared_node_indices[0], shared_node_id);
            new_elems[adj_child_elem_indices[1]].set_node(shared_node_indices[1], shared_node_id);

            // create a new edge between the shared node and central node
            let new_edge_id = self.new_edge_between_nodes(
                [shared_node_id, center_node_id],
                edge_id_tracker,
                parent_elem_id,
            )?;

            // connect it to the child Elems
            new_elems[adj_child_elem_indices[0]].set_edge(internal_edge_idx[0], new_edge_id);
            new_elems[adj_child_elem_indices[1]].set_edge(internal_edge_idx[1], new_edge_id);
        }

        // upgrade the ElemUninits to Elems (They should each have 4 node_ids and 4 edge_ids by this point)
        // connect the Elems to their relevant edges in the process
        self.upgrade_uninit_elems(new_elems)
    }

    fn execute_u_refinement(
        &mut self,
        mut new_elems: Vec<ElemUninit>,
        parent_elem_id: usize,
        node_id_tracker: &mut IdTracker,
        edge_id_tracker: &mut IdTracker,
    ) -> Result<Vec<Elem>, HRefError> {
        assert_eq!(new_elems.len(), 2);
        let mut outer_node_ids = [0; 2];

        // h_refine Edges 0 and 1 if necessary and connect child Elems to the relevant Nodes and Edges
        for (edge_index, shared_node_indices, outer_node_indices) in
            [(0, [1, 0], [0, 1]), (1, [3, 2], [2, 3])]
        {
            // get ids of child Edges and Node. Create them if they haven't been already
            let (child_edge_ids, shared_node_id) = self.h_refine_edge_if_needed(
                self.elems[parent_elem_id].edges[edge_index],
                node_id_tracker,
                edge_id_tracker,
            )?;

            outer_node_ids[edge_index] = shared_node_id;

            // connect the new Edges and Node to the child Elems
            new_elems[0].set_edge(edge_index, child_edge_ids[0]);
            new_elems[1].set_edge(edge_index, child_edge_ids[1]);
            new_elems[0].set_node(shared_node_indices[0], shared_node_id);
            new_elems[1].set_node(shared_node_indices[1], shared_node_id);

            // connect the parents outer Nodes to the child Elems
            new_elems[0].set_node(
                outer_node_indices[0],
                self.elems[parent_elem_id].nodes[outer_node_indices[0]],
            );
            new_elems[1].set_node(
                outer_node_indices[1],
                self.elems[parent_elem_id].nodes[outer_node_indices[1]],
            );
        }

        // create a new Edge between the two new Nodes
        let new_edge_id =
            self.new_edge_between_nodes(outer_node_ids, edge_id_tracker, parent_elem_id)?;

        // connect the new Edge to both child Elems
        new_elems[0].set_edge(3, new_edge_id);
        new_elems[1].set_edge(2, new_edge_id);

        // connect the parents outer unrefined Edges to the child Elems
        new_elems[0].set_edge(2, self.elems[parent_elem_id].edges[2]);
        new_elems[1].set_edge(3, self.elems[parent_elem_id].edges[3]);

        // upgrade the ElemUninits to Elems (They should each have 4 node_ids and 4 edge_ids by this point)
        // connect the Elems to their relevant edges in the process
        self.upgrade_uninit_elems(new_elems)
    }

    fn execute_v_refinement(
        &mut self,
        mut new_elems: Vec<ElemUninit>,
        parent_elem_id: usize,
        node_id_tracker: &mut IdTracker,
        edge_id_tracker: &mut IdTracker,
    ) -> Result<Vec<Elem>, HRefError> {
        assert_eq!(new_elems.len(), 2);
        let mut outer_node_ids = [0; 2];

        // h_refine Edges 0 and 1 if necessary and connect child Elems to the relevant Nodes and Edges
        for (edge_index, shared_node_indices, outer_node_indices) in
            [(2, [2, 0], [0, 2]), (3, [3, 1], [1, 3])]
        {
            // get ids of child Edges and Node. Create them if they haven't been already
            let (child_edge_ids, shared_node_id) = self.h_refine_edge_if_needed(
                self.elems[parent_elem_id].edges[edge_index],
                node_id_tracker,
                edge_id_tracker,
            )?;

            outer_node_ids[edge_index - 2] = shared_node_id;

            // connect the new Edges and Node to the child Elems
            new_elems[0].set_edge(edge_index, child_edge_ids[0]);
            new_elems[1].set_edge(edge_index, child_edge_ids[1]);
            new_elems[0].set_node(shared_node_indices[0], shared_node_id);
            new_elems[1].set_node(shared_node_indices[1], shared_node_id);

            // connect the parents outer Nodes to the child Elems
            new_elems[0].set_node(
                outer_node_indices[0],
                self.elems[parent_elem_id].nodes[outer_node_indices[0]],
            );
            new_elems[1].set_node(
                outer_node_indices[1],
                self.elems[parent_elem_id].nodes[outer_node_indices[1]],
            );
        }

        // create a new Edge between the two new Nodes
        let new_edge_id =
            self.new_edge_between_nodes(outer_node_ids, edge_id_tracker, parent_elem_id)?;

        // connect the new Edge to both child Elems
        new_elems[0].set_edge(1, new_edge_id);
        new_elems[1].set_edge(0, new_edge_id);

        // connect the parents outer unrefined Edges to the child Elems
        new_elems[0].set_edge(0, self.elems[parent_elem_id].edges[0]);
        new_elems[1].set_edge(1, self.elems[parent_elem_id].edges[1]);

        // upgrade the ElemUninits to Elems (They should each have 4 node_ids and 4 edge_ids by this point)
        // connect the Elems to their relevant edges in the process
        self.upgrade_uninit_elems(new_elems)
    }

    fn h_refine_edge_if_needed(
        &mut self,
        parent_edge_id: usize,
        node_id_tracker: &mut IdTracker,
        edge_id_tracker: &mut IdTracker,
    ) -> Result<(SmallVec<[usize; 2]>, usize), HRefError> {
        Ok(if self.edges[parent_edge_id].has_children() {
            (
                self.edges[parent_edge_id].child_ids().unwrap(),
                self.edges[parent_edge_id].child_node_id().unwrap(),
            )
        } else {
            let new_edge_ids = edge_id_tracker.next_two_ids();
            let new_node_id = node_id_tracker.next_id();

            let mut new_edges = self.edges[parent_edge_id].h_refine(new_edge_ids, new_node_id)?;
            self.edges.extend(new_edges.drain(0..));

            let parent_edge_points = self.edge_points(parent_edge_id);
            let node_coords = Point::between(parent_edge_points[0], parent_edge_points[1]);

            assert_eq!(new_node_id, self.nodes.len());

            self.nodes.push(Node::new(
                new_node_id,
                node_coords,
                self.edges[parent_edge_id].boundary,
            ));

            (SmallVec::from(new_edge_ids), new_node_id)
        })
    }

    fn new_edge_between_nodes(
        &mut self,
        node_ids: [usize; 2],
        edge_id_tracker: &mut IdTracker,
        parent_elem_id: usize,
    ) -> Result<usize, HRefError> {
        assert_ne!(node_ids[0], node_ids[1]);

        let new_edge_id = edge_id_tracker.next_id();
        let node_0 = &self.nodes[node_ids[0]];
        let node_1 = &self.nodes[node_ids[1]];

        let ordered_nodes = match self.elems[parent_elem_id]
            .element
            .order_points(&node_0.coords, &node_1.coords)
        {
            Ordering::Equal => return Err(HRefError::EdgeOnEqualPoints(parent_elem_id)),
            Ordering::Less => [node_0, node_1],
            Ordering::Greater => [node_1, node_0],
        };

        // print!("Edge {} \t [{} --- {}] \t", new_edge_id, ordered_nodes[0].coords, ordered_nodes[1].coords);
        let new_edge = Edge::new(new_edge_id, ordered_nodes, false);
        // println!("Dir: {:?}", new_edge.dir);

        self.edges.push(new_edge);

        Ok(new_edge_id)
    }

    fn upgrade_uninit_elems(
        &mut self,
        elems_uninit: Vec<ElemUninit>,
    ) -> Result<Vec<Elem>, HRefError> {
        let mut elems = Vec::with_capacity(4);
        for elem_uninit in elems_uninit {
            elems.push(elem_uninit.into_elem()?);
        }

        for elem in elems.iter() {
            for edge_id in elem.edges {
                self.edges[edge_id].connect_elem(elem);
            }

            for node_id in elem.nodes {
                self.nodes[node_id].connect_elem(elem);
            }
        }

        Ok(elems)
    }

    pub(crate) fn set_edge_activation(&mut self) {
        for edge in self.edges.iter_mut() {
            edge.reset_activation();
        }

        let mut base_edge_ids: Vec<usize> = self
            .edges
            .iter()
            .filter(|edge| edge.parent_id().is_none() && !edge.boundary)
            .map(|edge| edge.id)
            .collect();

        for base_edge_id in base_edge_ids.drain(0..) {
            if !self.rec_set_edge_activation_in_tree(base_edge_id) {
                panic!("Unable to find active Edge pair over Edge {}; Something must be wrong with the mesh!", base_edge_id);
            }
        }
    }

    fn rec_set_edge_activation_in_tree(&mut self, edge_id: usize) -> bool {
        if self.edges[edge_id].set_activation() {
            if let Some(child_edge_ids) = self.edges[edge_id].child_ids() {
                match (
                    self.rec_set_edge_activation_in_tree(child_edge_ids[0]), 
                    self.rec_set_edge_activation_in_tree(child_edge_ids[1]),
                ) { 
                    (true, true) => self.edges[edge_id].reset_activation(),
                    (false, false) => (),
                    _ => panic!("Children of Edge {} do not have consistent support for Basis Functions; Cannot set activation states!", edge_id),
                };
            }
            true
        } else {
            false
        }
    }

    // ----------------------------------------------------------------------------------------------------
    // p-refinement methods
    // ----------------------------------------------------------------------------------------------------

    /// Apply a [PRef] to all [Elem]s
    pub fn global_p_refinement(&mut self, refinement: PRef) -> Result<(), PRefError> {
        self.execute_p_refinements(
            self.elems
                .iter()
                .map(|elem| (elem.id, refinement))
                .collect(),
        )
    }

    /// Apply a [PRef] to a list of [Elem]s by their ID
    pub fn p_refine_elems(&mut self, elems: Vec<usize>, refinement: PRef) -> Result<(), PRefError> {
        self.execute_p_refinements(elems.iter().map(|elem_id| (*elem_id, refinement)).collect())
    }

    /// p-refine [Elem]s according to an external filter function
    pub fn p_refine_with_filter<F>(&mut self, filt: F) -> Result<(), PRefError>
    where
        F: Fn(&Elem) -> Option<PRef>,
    {
        self.execute_p_refinements(
            self.elems
                .iter()
                .map(|elem| (elem.id, filt(elem)))
                .filter(|(_, refinement)| refinement.is_some())
                .map(|(id, refinement)| (id, refinement.unwrap()))
                .collect(),
        )
    }

    /// Execute a series of [PRef]s on [Elem]s specified by their id
    pub fn execute_p_refinements(
        &mut self,
        refinements: Vec<(usize, PRef)>,
    ) -> Result<(), PRefError> {
        let mut refinements_map: BTreeMap<usize, PRef> = BTreeMap::new();
        for (elem_id, p_ref) in refinements {
            if elem_id >= self.elems.len() {
                return Err(PRefError::ElemDoesntExist(elem_id));
            }
            if refinements_map.insert(elem_id, p_ref).is_some() {
                return Err(PRefError::DoubleRefinement(elem_id));
            }
        }

        for (elem_id, refinement) in refinements_map {
            self.elems[elem_id].poly_orders.refine(refinement)?;
        }

        Ok(())
    }

    /// Set the expansion orders on all [Elem]s
    pub fn set_global_expansion_orders(&mut self, orders: [u8; 2]) -> Result<(), PRefError> {
        self.set_expansion_orders(self.elems.iter().map(|elem| (elem.id, orders)).collect())
    }

    /// Set the expansion orders on a list of [Elem]s by their ID
    pub fn set_expansion_on_elems(
        &mut self,
        elems: Vec<usize>,
        orders: [u8; 2],
    ) -> Result<(), PRefError> {
        self.set_expansion_orders(elems.iter().map(|elem_id| (*elem_id, orders)).collect())
    }

    /// set expansion orders on [Elem]s according to an external filter function
    pub fn set_expansions_with_filter<F>(&mut self, filt: F) -> Result<(), PRefError>
    where
        F: Fn(&Elem) -> Option<[u8; 2]>,
    {
        self.set_expansion_orders(
            self.elems
                .iter()
                .map(|elem| (elem.id, filt(elem)))
                .filter(|(_, orders)| orders.is_some())
                .map(|(id, orders)| (id, orders.unwrap()))
                .collect(),
        )
    }

    pub fn set_expansion_orders(
        &mut self,
        poly_orders: Vec<(usize, [u8; 2])>,
    ) -> Result<(), PRefError> {
        let mut poly_orders_map: BTreeMap<usize, [u8; 2]> = BTreeMap::new();
        for (elem_id, orders) in poly_orders {
            if elem_id >= self.elems.len() {
                return Err(PRefError::ElemDoesntExist(elem_id));
            }
            if poly_orders_map.insert(elem_id, orders).is_some() {
                return Err(PRefError::DoubleRefinement(elem_id));
            }
        }

        for (elem_id, orders) in poly_orders_map {
            self.elems[elem_id].poly_orders.set(orders)?;
        }

        Ok(())
    }
}

// ----------------------------------------------------------------------------------------------------
// Mesh construction from JSON Utility functions
// ----------------------------------------------------------------------------------------------------

/*
    edge - node_pair - side relationships

    edge 0 : [node_0, node_1], top
    edge 1 : [node_2, node_3], bottom,
    edge 2 : [node_0, node_2], right,
    edge 3 : [node_1, node_3], left,
*/
const EDGE_IDX_DEFS: [([usize; 2], usize); 4] =
    [([0, 1], 1), ([2, 3], 0), ([0, 2], 1), ([1, 3], 0)];

fn parse_element_information(mesh_file_json: &JsonValue) -> (Vec<Materials>, Vec<[usize; 4]>) {
    assert!(
        mesh_file_json["Elements"].is_array(),
        "Elements must be an Array!"
    );

    let num_nodes = mesh_file_json["Nodes"].members().count();

    mesh_file_json["Elements"]
        .members()
        .map(|json_element| {
            assert!(
                json_element["node_ids"].is_array(),
                "Elements must have an Array of node_ids!"
            );
            assert_eq!(
                json_element["node_ids"].members().count(),
                4,
                "Elements Array of node_ids must have a length of 4!"
            );

            assert!(
                json_element["materials"].is_array(),
                "Elements must have an Array of materials!"
            );
            assert_eq!(
                json_element["materials"].members().count(),
                4,
                "Elements Array of materials must have a length of 4!"
            );

            let node_ids: [usize; 4] = json_element["node_ids"]
                .members()
                .map(|node_id_json| {
                    let node_id = node_id_json
                        .as_usize()
                        .expect("node_ids must be positive integers!");
                    assert!(
                        node_id < num_nodes,
                        "node_ids must be smaller than the total number of nodes!"
                    );
                    node_id
                })
                .collect::<Vec<usize>>()
                .try_into()
                .unwrap();
            assert!(
                !has_duplicates(&node_ids),
                "Element's node_ids should have 4 unique values!"
            );

            let material_props: [f64; 4] = json_element["materials"]
                .members()
                .map(|mp_json| {
                    mp_json
                        .as_f64()
                        .expect("Element materials must be numerical values")
                })
                .collect::<Vec<f64>>()
                .try_into()
                .unwrap();

            (Materials::from_array(material_props), node_ids)
        })
        .unzip()
}

fn parse_node_information(mesh_file_json: &JsonValue) -> Vec<Point> {
    assert!(
        mesh_file_json["Nodes"].is_array(),
        "Nodes must be an Array!"
    );

    let node_points: Vec<Point> = mesh_file_json["Nodes"]
        .members()
        .map(|json_node_point| {
            assert!(json_node_point.is_array(), "nodes must be arrays!");
            assert_eq!(
                json_node_point.members().count(),
                2,
                "nodes must be arrays of length 2!"
            );

            let x = json_node_point[0]
                .as_f64()
                .expect("nodes must be composed of numerical values!");
            let y = json_node_point[1]
                .as_f64()
                .expect("nodes must be composed of numerical values!");

            Point::new(x, y)
        })
        .collect();

    assert!(
        !has_duplicates(&node_points),
        "All Nodes must be at unique locations!"
    );

    node_points
}

fn has_duplicates<T>(values: &[T]) -> bool
where
    T: PartialEq,
{
    for (i, val) in values.iter().enumerate() {
        for val_cmp in values.iter().skip(i + 1) {
            if val == val_cmp {
                return true;
            }
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    const MESH_A_POINTS_X: [[f64; 4]; 4] = [
        [0.0, 1.0, 0.0, 1.0],
        [1.0, 2.0, 1.0, 2.0],
        [0.0, 1.0, 0.0, 1.0],
        [1.0, 2.0, 1.0, 2.0],
    ];

    const MESH_A_POINTS_Y: [[f64; 4]; 4] = [
        [0.0, 0.0, 0.5, 0.5],
        [0.0, 0.0, 0.5, 0.5],
        [0.5, 0.5, 1.0, 1.0],
        [0.5, 0.5, 1.0, 1.0],
    ];

    const MESH_A_NEIGHBORS: [[Option<usize>; 4]; 4] = [
        [None, Some(2), None, Some(1)],
        [None, Some(3), Some(0), None],
        [Some(0), None, None, Some(3)],
        [Some(1), None, Some(2), None],
    ];

    #[test]
    fn mesh_from_file() {
        let mut mesh_a = Mesh::from_file("./test_input/test_mesh_a.json").unwrap();
        mesh_a.set_edge_activation();

        for (element, mats_cmp) in mesh_a
            .elements
            .iter()
            .zip([[1.0, 1.0], [1.0, 2.0], [2.0, 1.0], [2.0, 2.0]].iter())
        {
            assert!((element.materials.eps_rel.re - mats_cmp[0]).abs() < 1e-14);
        }

        for (elem_id, elem) in mesh_a.elems.iter().enumerate() {
            for n_idx in 0..4 {
                let p = &mesh_a.nodes[elem.nodes[n_idx]].coords;

                assert!((p.x - MESH_A_POINTS_X[elem_id][n_idx]).abs() < 1e-14);
                assert!((p.y - MESH_A_POINTS_Y[elem_id][n_idx]).abs() < 1e-14);

                if let Some(expected_neighbor) = MESH_A_NEIGHBORS[elem_id][n_idx] {
                    assert_eq!(
                        mesh_a.edges[elem.edges[n_idx]]
                            .other_active_elem_id(elem_id)
                            .unwrap(),
                        expected_neighbor
                    );
                }
            }
        }

        for element in mesh_a.elements.iter() {
            assert!(element.points[3].x > element.points[0].x);
            assert!(element.points[3].x >= element.points[1].x);
            assert!(element.points[3].x > element.points[2].x);

            assert!(element.points[3].y > element.points[0].y);
            assert!(element.points[3].y > element.points[1].y);
            assert!(element.points[3].y >= element.points[2].y);

            assert!(element.points[0].x < element.points[1].x);
            assert!(element.points[0].y < element.points[2].y);
        }
    }

    #[test]
    fn basic_h_refinements() {
        let mut mesh_c = Mesh::from_file("./test_input/test_mesh_c.json").unwrap();
        mesh_c.h_refine_elems(vec![0], HRef::T).unwrap();
        mesh_c.h_refine_elems(vec![1, 2], HRef::U(None)).unwrap();
        mesh_c.h_refine_elems(vec![3, 4], HRef::V(None)).unwrap();
    }

    #[test]
    fn basic_p_refinements() {
        let mut mesh_b = Mesh::from_file("./test_input/test_mesh_b.json").unwrap();
        mesh_b.p_refine_elems(vec![0], PRef::from(2, 2)).unwrap();
        mesh_b.p_refine_elems(vec![1], PRef::from(2, 1)).unwrap();
        mesh_b.p_refine_elems(vec![2], PRef::from(1, 2)).unwrap();

        assert_eq!(mesh_b.elems[0].poly_orders.ni, 3);
        assert_eq!(mesh_b.elems[0].poly_orders.nj, 3);
        assert_eq!(mesh_b.elems[1].poly_orders.ni, 3);
        assert_eq!(mesh_b.elems[1].poly_orders.nj, 2);
        assert_eq!(mesh_b.elems[2].poly_orders.ni, 2);
        assert_eq!(mesh_b.elems[2].poly_orders.nj, 3);
    }

    #[test]
    fn refined_mesh_to_file() {
        let mut mesh_b = Mesh::from_file("./test_input/test_mesh_b.json").unwrap();
        mesh_b
            .export_to_json("./test_output/mesh_b_copy.json")
            .unwrap();

        mesh_b
            .execute_p_refinements(vec![
                (0, PRef::from(3, 3)),
                (1, PRef::from(2, 3)),
                (2, PRef::from(3, 2)),
            ])
            .unwrap();
        mesh_b
            .execute_h_refinements(vec![
                (0, HRef::T),
                (1, HRef::u_extened(0).unwrap()),
                (2, HRef::v_extened(0).unwrap()),
            ])
            .unwrap();

        mesh_b
            .p_refine_elems(vec![3, 4, 5], PRef::from(-1, -1))
            .unwrap();
        mesh_b.h_refine_elems(vec![6, 14, 12], HRef::T).unwrap();

        mesh_b.set_edge_activation();
        mesh_b
            .export_to_json("./test_output/mesh_b_refined.json")
            .unwrap();
    }

    #[test]
    #[should_panic]
    fn h_refine_non_existent() {
        let mut mesh_c = Mesh::from_file("./test_input/test_mesh_c.json").unwrap();
        mesh_c.h_refine_elems(vec![0, 1], HRef::T).unwrap();
    }

    #[test]
    #[should_panic]
    fn h_refine_elem_with_children() {
        let mut mesh_c = Mesh::from_file("../test_input/test_mesh_c.json").unwrap();
        mesh_c.h_refine_elems(vec![0], HRef::T).unwrap();
        mesh_c.h_refine_elems(vec![0], HRef::T).unwrap();
    }

    #[test]
    #[should_panic]
    fn double_h_refinement() {
        let mut mesh_a = Mesh::from_file("./test_input/test_mesh_c.json").unwrap();
        mesh_a
            .execute_h_refinements(vec![(0, HRef::T), (1, HRef::T), (0, HRef::u())])
            .unwrap();
    }

    #[test]
    #[should_panic]
    fn bad_extended_h_refinement() {
        let mut mesh_c = Mesh::from_file("./test_input/test_mesh_c.json").unwrap();
        mesh_c
            .h_refine_elems(vec![0], HRef::u_extened(2).unwrap())
            .unwrap();
    }

    #[test]
    #[should_panic]
    fn minimum_edge_length_exceeded() {
        let mut mesh_c = Mesh::from_file("./test_input/test_mesh_c.json").unwrap();

        // repeatedly refine the bottom left cell
        for _ in 0..18 {
            mesh_c
                .h_refine_with_filter(|elem| {
                    if !elem.has_children() && elem.nodes[0] == 0 {
                        Some(HRef::T)
                    } else {
                        None
                    }
                })
                .unwrap()
        }
    }

    #[test]
    #[should_panic]
    fn p_refine_non_existent() {
        let mut mesh_c = Mesh::from_file("./test_input/test_mesh_c.json").unwrap();
        mesh_c.p_refine_elems(vec![0, 1], PRef::from(1, 1)).unwrap();
    }

    #[test]
    #[should_panic]
    fn double_p_refinement() {
        let mut mesh_c = Mesh::from_file("./test_input/test_mesh_c.json").unwrap();
        mesh_c.p_refine_elems(vec![0, 0], PRef::from(1, 1)).unwrap();
    }

    #[test]
    #[should_panic]
    fn neg_p_refinement_i() {
        let mut mesh_c = Mesh::from_file("./test_input/test_mesh_c.json").unwrap();
        mesh_c.set_global_expansion_orders([3, 3]).unwrap();
        mesh_c.p_refine_elems(vec![0], PRef::from(-3, 1)).unwrap();
    }

    #[test]
    #[should_panic]
    fn neg_p_refinement_j() {
        let mut mesh_c = Mesh::from_file("./test_input/test_mesh_c.json").unwrap();
        mesh_c.set_global_expansion_orders([3, 3]).unwrap();

        mesh_c.p_refine_elems(vec![0], PRef::from(1, -3)).unwrap();
    }

    #[test]
    #[should_panic]
    fn p_refinement_over_max_i() {
        let mut mesh_c = Mesh::from_file("./test_input/test_mesh_c.json").unwrap();
        let max_exp_as_i8 = MAX_POLYNOMIAL_ORDER.try_into().unwrap();
        mesh_c
            .p_refine_elems(vec![0], PRef::from(max_exp_as_i8, 0))
            .unwrap();
    }

    #[test]
    #[should_panic]
    fn p_refinement_over_max_j() {
        let mut mesh_c = Mesh::from_file("../../test_input/test_mesh_c.json").unwrap();
        let max_exp_as_i8 = MAX_POLYNOMIAL_ORDER.try_into().unwrap();
        mesh_c
            .p_refine_elems(vec![0], PRef::from(0, max_exp_as_i8))
            .unwrap();
    }
}
