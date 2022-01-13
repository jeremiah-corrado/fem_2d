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
pub use h_refinement::{HRef, HRefError, HRefLoc, Quadrant, Bisection};
pub use p_refinement::{PRef, PRefError};
pub use node::Node;
pub use space::{ParaDir, Point, M2D, V2D};

use std::rc::Rc;
use std::fs::read_to_string;
use json::JsonValue;

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
}

impl Mesh {
    pub fn blank() -> Self {
        Self {
            elements: Vec::new(),
            elems: Vec::new(),
            nodes: Vec::new(),
            edges: Vec::new(),
        }
    }

    pub fn from_file(path: impl AsRef<str>) -> std::io::Result<Self> {
        let mesh_file_contents = read_to_string(path.as_ref())?;
        let mesh_file_json = json::parse(&mesh_file_contents).expect("Unable to parse Mesh File as JSON!");

        let (mut materials, node_ids, points) = parse_mesh_json(mesh_file_json);

        let elements = materials.drain(0..).zip(node_ids.iter()).enumerate().map(|(element_id, (element_materials, element_node_ids))| {
            Rc::new(Element::new(element_id, element_node_ids.iter().map(|node_id| {
                points[*node_id].clone()
            }).collect::<Vec<_>>().try_into().unwrap(), element_materials))
        }).collect();

        let mut node_connection_counts = vec![0; points.len()];
        for element_node_ids in node_ids.iter() {
            for node_id in element_node_ids.iter() {
                node_connection_counts[*node_id] += 1;
            }
        }
        let node_boundary_statuses: Vec<bool> = node_connection_counts.iter().map(|count| {
            assert!(*count <= 4, "Nodes should be shared by a maximum of 4 Elements!");
            *count < 4
        }).collect();

        let nodes = points.iter().enumerate().map(|(node_id, point)| {
            Node::new(node_id, *point, node_boundary_statuses[node_id])
        }).collect();

        let mut points_by_u: Vec<(usize, Point)> = points.iter().enumerate().map(|(node_id, p)| {
            (node_id, *p)
        }).collect();
        let mut points_by_v = points_by_u.clone();

        Ok(Self{
            elements,
            elems: Vec::new(),
            nodes,
            edges: Vec::new(),
        })
    }
}



fn parse_mesh_json(mesh_file_json: JsonValue) -> (Vec<Materials>, Vec<[usize; 4]>, Vec<Point>) {
    assert!(mesh_file_json["Elements"].is_array(), "Elements must be an Array!");
    assert!(mesh_file_json["Nodes"].is_array(), "Nodes must be an Array!");

    let num_nodes = mesh_file_json["Nodes"].members().count();

    let (element_materials, element_node_ids) : (Vec<Materials>, Vec<[usize; 4]>) = mesh_file_json["Elements"].members().map(|json_element| {
        assert!(json_element["node_ids"].is_array(), "Elements must have an Array of node_ids!");
        assert_eq!(json_element["node_ids"].members().count(), 4, "Elements Array of node_ids must have a length of 4!");

        assert!(json_element["materials"].is_array(), "Elements must have an Array of materials!");
        assert_eq!(json_element["materials"].members().count(), 4, "Elements Array of materials must have a length of 4!");

        let node_ids: [usize; 4] = json_element["node_ids"].members().map(|node_id_json| {
            let node_id = node_id_json.as_usize().expect("node_ids must be positive integers!");
            assert!(node_id < num_nodes, "node_ids must be smaller than the total number of nodes!");
            node_id
        }).collect::<Vec<usize>>().try_into().unwrap();
        assert!(!has_duplicates(&node_ids), "Element's node_ids should have 4 unique values!");

        let material_props: [f64; 4] = json_element["materials"].members().map(|mp_json| {
            mp_json.as_f64().expect("Element materials must be numerical values")
        }).collect::<Vec<f64>>().try_into().unwrap();

        (Materials::from_array(material_props), node_ids)
    }).unzip();

    let node_points: Vec<Point> = mesh_file_json["Nodes"].members().map(|json_node_point| {
        assert!(json_node_point.is_array(), "nodes must be arrays!");
        assert_eq!(json_node_point.members().count(), 2, "nodes must be arrays of length 2!");

        let x = json_node_point[0].as_f64().expect("nodes must be composed of numerical values!");
        let y = json_node_point[1].as_f64().expect("nodes must be composed of numerical values!");

        Point::new(x, y)
    }).collect();

    assert!(!has_duplicates(&node_points), "All Nodes must be at unique locations!");

    (element_materials, element_node_ids, node_points)
}

fn has_duplicates<T>(values: &[T]) -> bool where T: PartialEq {
    for (i, val) in values.iter().enumerate() {
        for (val_cmp) in values.iter().skip(i) {
            if val == val_cmp {
                return true
            }
        }
    }
    false
}
