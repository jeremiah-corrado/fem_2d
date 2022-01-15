use super::{
    h_refinement::{Bisection, HRefError},
    Elem, Node, ParaDir, MIN_EDGE_LENGTH,
};
use json::JsonValue;
use smallvec::SmallVec;
use std::collections::BTreeMap;

#[derive(Debug)]
/// The connection between two neighboring [Elem]s. Defined by two points in real space.
pub struct Edge {
    pub id: usize,
    pub nodes: [usize; 2],
    pub boundary: bool,
    pub dir: ParaDir,
    pub length: f64,
    children: Option<[usize; 2]>,
    parent: Option<(usize, Bisection)>,
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

    /// Connect an [Elem] to this node based on their relative orientation in space and the h-refinement state of the Elem
    pub fn connect_elem(&mut self, elem: &Elem) {
        let index_of_self = match elem.edges.iter().position(|edge_id| edge_id == &self.id) {
            Some(index) => index,
            None => panic!(
                "Elem {} does not have Edge {}; Cannot form connection!",
                elem.id, self.id
            ),
        };

        let side_index = match index_of_self {
            0 | 2 => 1,
            1 | 3 => 0,
            _ => unreachable!(),
        };
        let address = elem.h_levels.edge_ranking(self.dir);
        assert!(
            self.elems[side_index].get(&address).is_none(),
            "Edge {} is already connected another Elem at {:?}; Cannot connect to Elem {}!",
            self.id,
            address,
            elem.id
        );

        self.elems[side_index].insert(address, elem.id);
    }

    /// Produce two child Edges from this edge and a new Node along its center
    pub fn h_refine(
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
                    self.children = Some(new_ids.clone());
                    self.child_node = Some(new_node_id);
                    Ok(smallvec![
                        Self {
                            id: new_ids[0],
                            nodes: [self.nodes[0], new_node_id],
                            boundary: self.boundary,
                            dir: self.dir,
                            length: child_edge_length,
                            children: None,
                            parent: Some((self.id, Bisection::BL)),
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
                            parent: Some((self.id, Bisection::TR)),
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
        match self.parent {
            Some((parent_id, _)) => Some(parent_id),
            None => None,
        }
    }

    /// Returns a vector of child Edge ids. Will return an empty vector if this Edge has no children.
    pub fn child_ids(&self) -> Vec<usize> {
        match self.children {
            Some(child_ids) => child_ids.to_vec(),
            None => Vec::new(),
        }
    }

    pub fn has_children(&self) -> bool {
        self.children.is_some()
    }

    pub fn child_node_id(&self) -> Option<usize> {
        self.child_node
    }

    /// Produce a Json Object that describes this Elem
    pub fn to_json(&self) -> JsonValue {
        let mut edge_json = object! {
            "id": self.id,
            "boundary": self.boundary,
            "direction": self.dir,
            "nodes": array![self.nodes[0], self.nodes[1]],
            "parent": self.parent_id(),
            "children": JsonValue::from(self.child_ids()),
            "cells": array![array![], array![]],
        };

        for side_idx in 0..2 {
            for (elem_key, elem_id) in self.elems[side_idx].iter() {
                edge_json["cells"][side_idx]
                    .push(object! {
                        "level_key": array![elem_key[0], elem_key[1]],
                        "cell_id": *elem_id,
                    })
                    .unwrap();
            }
        }

        if let Some([active_elem_bl, active_elem_tr]) = self.active_elems {
            edge_json["active_cells"] = array![active_elem_bl, active_elem_tr];
        }

        edge_json
    }
}
