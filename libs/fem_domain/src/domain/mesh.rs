mod edge;
mod elem;
mod element;
mod h_refinement;
mod node;
mod p_refinement;
mod space;

pub use edge::Edge;
pub use elem::{Elem, ElemUninit};
pub use element::{Element, Materials};
pub use h_refinement::{Bisection, HRef, HRefError, HRefLoc, Quadrant, HLevels};
pub use node::Node;
pub use p_refinement::{PRef, PRefError, PolyOrders};
pub use space::{ParaDir, Point, M2D, V2D};

use json::JsonValue;

use std::collections::BTreeMap;
use std::fs::read_to_string;
use std::rc::Rc;

/// Minimum Edge length in parametric space. h-Refinements will fail after edges are smaller than this value.
pub const MIN_EDGE_LENGTH: f64 = 1e-6;

/// Maximum Polynomial expansion. p-Refinements will fail when Elem's expansion orders exceed this value.
pub const MAX_POLYNOMIAL_ORDER: u8 = 20;

/// Information used to Define the geometric structure and refinement state of a Domain.
pub struct Mesh {
    pub elements: Vec<Rc<Element>>,
    pub elems: Vec<Elem>,
    pub nodes: Vec<Node>,
    pub edges: Vec<Edge>,
    refinement_generations: Vec<RefinementGeneration>,
}

struct RefinementGeneration {
    pub open: bool,
    pub elem_id_range: [usize; 2],
}

impl RefinementGeneration {
    pub fn new(starting_id: usize, ending_id: usize) -> Self {
        Self {
            open: true,
            elem_id_range: [starting_id, ending_id],
        }
    }

    pub fn next_elem_id(&mut self) -> usize {
        assert!(self.open, "Cannot increment elem_id on a closed refinement generation");
        self.elem_id_range[1] += 1;
        self.elem_id_range[1]
    }

    pub fn next_gen(&mut self) -> Self {
        self.open = false;
        Self::new(self.elem_id_range[1], self.elem_id_range[1])
    }
}


impl Mesh {
    /// Construct a completely empty Mesh
    pub fn blank() -> Self {
        Self {
            elements: Vec::new(),
            elems: Vec::new(),
            nodes: Vec::new(),
            edges: Vec::new(),
            refinement_generations: Vec::new(),
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
        let elements: Vec<Rc<Element>> = element_materials
            .drain(0..)
            .zip(element_node_ids.iter())
            .enumerate()
            .map(|(element_id, (materials, node_ids))| {
                Rc::new(Element::new(
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
        let nodes: Vec<Node> = points
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
        for (edge_id, elements) in edge_node_pairs.values().enumerate() {
            for (elem_side_idx, elem_id) in elements
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

                for edge_id in edge_ids.iter() {
                    edges[*edge_id].connect_elem(&elem)
                }

                elem
            })
            .collect();

        let ref_gen_0 = RefinementGeneration::new(0, elems.len() - 1);

        Ok(Self {
            elements,
            elems,
            nodes,
            edges,
            refinement_generations: vec![ref_gen_0]
        })
    }

    // ----------------------------------------------------------------------------------------------------
    // h-refinement methods
    // ----------------------------------------------------------------------------------------------------

    /// Execute a series of [HRef]s on [Elem]s specified by their id
    pub fn execute_h_refinements(
        &mut self,
        refinements: impl Iterator<Item = (usize, HRef)>,
    ) -> Result<(), HRefError> {
        let mut refinements_map: BTreeMap<usize, HRef> = BTreeMap::new();
        for (elem_id, h_ref) in refinements {
            if elem_id >= self.elems.len() {
                return Err(HRefError::ElemDoesntExist(elem_id));
            }
            if let Some(_) = refinements_map.insert(elem_id, h_ref) {
                return Err(HRefError::DoubleRefinement(elem_id))
            }
        }
        
        self.start_new_refinement_gen();
        let mut refinement_extensions: Vec<(usize, HRef)> = Vec::new();
        let mut elem_id_tracker = self.elems.len();

        for (elem_id, refinement) in refinements_map {
            let new_uninitialized_elems = self.elems[elem_id].h_refine(refinement, &mut elem_id_tracker)?;
            let mut new_elems = match refinement {
                HRef::T => self.execute_t_refinement(new_uninitialized_elems)?,
                HRef::U(extension) => {
                    let new_u_elems = self.execute_u_refinement(new_uninitialized_elems)?;
                    if let Some(v_ref) = extension {
                        refinement_extensions.push((new_u_elems[v_ref.index()].id, HRef::V(None)))
                    }
                    new_u_elems
                },
                HRef::V(extension) => {
                    let new_v_elems = self.execute_v_refinement(new_uninitialized_elems)?;
                    if let Some(u_ref) = extension {
                        refinement_extensions.push((new_v_elems[u_ref.index()].id, HRef::U(None)))
                    } 
                    new_v_elems
                }
            };

            self.elems.extend(new_elems.drain(0..));
        }

        if refinement_extensions.len() > 0 {
            self.execute_h_refinements(refinement_extensions.into_iter())?;
        }

        Ok(())
    }

    fn start_new_refinement_gen(&mut self)  {
        let next_gen = self.refinement_generations.last_mut().unwrap().next_gen();
        self.refinement_generations.push(next_gen);
    }

    fn execute_t_refinement(&mut self, new_elems: Vec<ElemUninit>) -> Result<Vec<Elem>, HRefError> {
        unimplemented!()
    }
    
    fn execute_u_refinement(&mut self, new_elems: Vec<ElemUninit>) -> Result<Vec<Elem>, HRefError> {
        unimplemented!()
    }
    
    fn execute_v_refinement(&mut self, new_elems: Vec<ElemUninit>) -> Result<Vec<Elem>, HRefError> {
        unimplemented!()
    }


    // ----------------------------------------------------------------------------------------------------
    // p-refinement methods
    // ----------------------------------------------------------------------------------------------------

    /// Execute a series of [PRef]s on [Elem]s specified by their id
    pub fn execute_p_refinements(&mut self, refinements: impl Iterator<Item = (usize, PRef)>) -> Result<(), PRefError> {
        let mut refinements_map: BTreeMap<usize, PRef> = BTreeMap::new();
        for (elem_id, p_ref) in refinements {
            if elem_id >= self.elems.len() {
                return Err(PRefError::ElemDoesntExist(elem_id));
            }
            if let Some(_) = refinements_map.insert(elem_id, p_ref) {
                return Err(PRefError::DoubleRefinement(elem_id));
            }
        }

        for (elem_id, refinement) in refinements_map {
            self.elems[elem_id].poly_orders.refine(refinement)?;
        }
        
        Ok(())
    }
}

// ----------------------------------------------------------------------------------------------------
// Helper methods for Mesh construction from JSON file
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
    use super::Mesh;

    #[test]
    fn mesh_from_file() {
        let mesh_a = Mesh::from_file("../../test_input/test_mesh_a.json").unwrap();
        let mesh_b = Mesh::from_file("../../test_input/test_mesh_b.json").unwrap();
    }
}
