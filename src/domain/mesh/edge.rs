use super::{elem::Elem, h_refinement::HRefError, node::Node, space::ParaDir, MIN_EDGE_LENGTH};
use json::{array, object, JsonValue};
use smallvec::{smallvec, SmallVec};
use std::collections::BTreeMap;

/// Edges describe a strait line in parametric space between two `Node`s
///
/// Edges keep track of all adjacent `Elem`s and are responsible for identifying whether or not they can support edge-type Degrees of Freedom.
/// They also decide which pair of adjacent `Elem`s are to be used for Shape Function matching if any.
///
/// ## Layout
/// There are two Edge variants: U-directed and V-directed
///
/// 1. **U-Dir**:
///
/// Edges are designated as u-directed if the angle between its two points is less than 45° relative to the Real x-axis
///
/// `Node` and `Elem` Indices:
/// ```text
///            T
///            1
///    0 *-----------* 1
///            0         
///            B
/// ```
///
/// Indices of child-edges after h-refinement:
/// ```text
///         0     1
///      *-----*-----*
/// ```
///     
/// 1. **V-Dir**:
///
/// Edges are designated as v-directed if the angle between its two points is greater than 45° relative to the Real x-axis
///
/// `Node` and `Elem` Indices:
/// ```text
///         1
///         *
///         |
///         |
///   L  0  |  1  R
///         |
///         |
///         *
///         0
///
///
/// ```
///
/// Indices of child-edges after h-refinement:
/// ```text
///         *
///         |
///         |  1
///         *
///         |
///         |  0
///         *  
/// ```
///
///
#[derive(Debug, Clone)]
pub struct Edge {
    pub id: usize,
    pub nodes: [usize; 2],
    pub boundary: bool,
    pub dir: ParaDir,
    pub length: f64,
    children: Option<[usize; 2]>,
    parent: Option<usize>,
    elems: [BTreeMap<[u8; 2], usize>; 2],
    active_elems: Option<[usize; 2]>,
    child_node: Option<usize>,
}

impl Edge {
    /// Construct a new edge between two points in real space
    pub fn new(id: usize, nodes: [&Node; 2], boundary: bool) -> Self {
        let dir = nodes[0].coords.orientation_with(&nodes[1].coords);

        Self {
            id,
            nodes: [nodes[0].id, nodes[1].id],
            boundary,
            dir,
            length: nodes[0].coords.dist(&nodes[1].coords),
            children: None,
            parent: None,
            elems: [BTreeMap::new(), BTreeMap::new()],
            active_elems: None,
            child_node: None,
        }
    }

    pub(crate) fn connect_elem(&mut self, elem: &Elem) {
        if let Some(index_of_self) = elem.edges.iter().position(|edge_id| edge_id == &self.id) {
            let address = elem.h_levels.edge_ranking(self.dir);
            let side_index = match index_of_self {
                0 | 2 => 1,
                1 | 3 => 0,
                _ => unreachable!(),
            };

            if let Some(prev_elem_id) = self.elems[side_index].insert(address, elem.id) {
                assert_eq!(
                    prev_elem_id, elem.id,
                    "Edge {} is already connected to Elem {} at {:?} (on side {}); cannot connect to Elem {}",
                    self.id,
                    prev_elem_id,
                    address,
                    side_index,
                    elem.id,
                );
            }
        } else {
            panic!(
                "Elem {} is not connected to Edge {}; cannot reciprocate connection!",
                elem.id, self.id
            );
        }
    }

    /// Produce two child Edges from this edge and connect them to a new Node along its center
    pub(crate) fn h_refine(
        &mut self,
        new_ids: [usize; 2],
        new_node_id: usize,
    ) -> Result<SmallVec<[Self; 2]>, HRefError> {
        match self.children {
            Some(_) => Err(HRefError::EdgeHasChildren(self.id)),
            None => {
                let child_edge_length = self.length / 2.0;

                if child_edge_length < MIN_EDGE_LENGTH {
                    Err(HRefError::MinEdgeLength(self.id))
                } else {
                    self.children = Some(new_ids);
                    self.child_node = Some(new_node_id);
                    Ok(smallvec![
                        Self {
                            id: new_ids[0],
                            nodes: [self.nodes[0], new_node_id],
                            boundary: self.boundary,
                            dir: self.dir,
                            length: child_edge_length,
                            children: None,
                            parent: Some(self.id),
                            elems: [BTreeMap::new(), BTreeMap::new()],
                            active_elems: None,
                            child_node: None,
                        },
                        Self {
                            id: new_ids[1],
                            nodes: [new_node_id, self.nodes[1]],
                            boundary: self.boundary,
                            dir: self.dir,
                            length: child_edge_length,
                            children: None,
                            parent: Some(self.id),
                            elems: [BTreeMap::new(), BTreeMap::new()],
                            active_elems: None,
                            child_node: None,
                        },
                    ])
                }
            }
        }
    }

    /// Id of the Parent Edge if this Edge has a parent
    pub fn parent_id(&self) -> Option<usize> {
        self.parent
    }

    /// Returns a vector of child Edge ids. Will return an empty vector if this Edge has no children.
    pub fn child_ids(&self) -> Option<SmallVec<[usize; 2]>> {
        self.children.map(SmallVec::from)
    }

    /// Has been h-Refined?
    pub fn has_children(&self) -> bool {
        self.children.is_some()
    }

    /// Returns the Node sitting on the center of this Edge if it has been h-Refined; otherwise, `None` is returned.
    pub fn child_node_id(&self) -> Option<usize> {
        self.child_node
    }

    /// Which two Elem's should support edge-type Shape Functions (if any)
    pub fn active_elem_pair(&self) -> Option<[usize; 2]> {
        self.active_elems
    }

    /// Is this edge on the border of the Mesh (does it have elements only on one side)
    pub fn is_boundary(&self) -> bool {
        self.elems[0].is_empty() || self.elems[1].is_empty()
    }

    /// Attempts to establish an active pair of Elems. Returns false if none can be established
    pub(crate) fn set_activation(&mut self) -> bool {
        match (self.last_entry(0), self.last_entry(1)) {
            (Some(bl_elem_id), Some(tr_elem_id)) => {
                self.active_elems = Some([bl_elem_id, tr_elem_id]);
                true
            }
            (_, _) => {
                self.active_elems = None;
                false
            }
        }
    }

    fn last_entry(&self, side_idx: usize) -> Option<usize> {
        if let Some((_, elem_id)) = self.elems[side_idx].iter().rev().take(1).next() {
            Some(*elem_id)
        } else {
            None
        }
    }

    pub(crate) fn reset_activation(&mut self) {
        self.active_elems = None;
    }

    /// Get the ID of the other active `Elem` connected to this Edge
    ///
    /// *Returns None if this Edge doesn't have a pair of Active `Elem`s
    /// *Panics if the edge has an active pair, but the `elem_id` argument isn't one of the active `Elem`s
    pub fn other_active_elem_id(&self, elem_id: usize) -> Option<usize> {
        match self.active_elems {
            Some(active_elem_ids) => {
                match active_elem_ids
                    .iter()
                    .position(|eid_cmp| eid_cmp == &elem_id)
                {
                    Some(this_pos) => match this_pos {
                        0 => Some(active_elem_ids[1]),
                        1 => Some(active_elem_ids[0]),
                        _ => unreachable!(),
                    },
                    None => panic!(
                        "{} isn't an active Elem on Edge {}; Cannot retrieve other active Elem!",
                        elem_id, self.id
                    ),
                }
            }
            None => None,
        }
    }

    /// Produce a Json Object that describes this Edge
    #[cfg(feature = "json_export")]
    pub fn to_json(&self) -> JsonValue {
        let mut edge_json = object! {
            "id": self.id,
            "boundary": self.boundary,
            "direction": self.dir,
            "nodes": array![self.nodes[0], self.nodes[1]],
            "parent": self.parent_id(),
            "children": match self.children {
                Some(child_ids) => array![child_ids[0], child_ids[1]],
                None => array![],
            },
            "elems": array![array![], array![]],
        };

        for side_idx in 0..2 {
            for (elem_key, elem_id) in self.elems[side_idx].iter() {
                edge_json["elems"][side_idx]
                    .push(object! {
                        "level_key": array![elem_key[0], elem_key[1]],
                        "cell_id": *elem_id,
                    })
                    .unwrap();
            }
        }

        if let Some([active_elem_bl, active_elem_tr]) = self.active_elems {
            edge_json["active_elems"] = array![active_elem_bl, active_elem_tr];
        } else {
            edge_json["active_elems"] = array![]
        }

        edge_json
    }
}
