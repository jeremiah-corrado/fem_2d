use super::{Point, Elem};
use json::JsonValue;
use std::collections::BTreeMap;

/// A point in 2D space.
/// Pairs of points describe Edges and Groups of 4 points describe Elems.
#[derive(Debug, Clone)]
pub struct Node {
    pub id: usize,
    pub coords: Point,
    pub boundary: bool,
    elems: [BTreeMap<[u8; 2], usize>; 4],
    active_elems: Option<[usize; 4]>,
}

impl Node {
    /// construct a new Node defined by its location in real space
    pub fn new(id: usize, coords: Point, boundary: bool) -> Self {
        Self {
            id,
            coords,
            boundary,
            elems: [
                BTreeMap::new(),
                BTreeMap::new(), 
                BTreeMap::new(), 
                BTreeMap::new(), 
            ],
            active_elems: None,
        }
    }

    pub(crate) fn connect_elem(&mut self, elem: &Elem) {
        if let Some(index_of_self) = elem.nodes.iter().position(|node_id| node_id == &self.id) {
            let address = elem.h_levels.node_ranking();

            if let Some(prev_elem_id) = self.elems[index_of_self].insert(address, elem.id) {
                assert_eq!(
                    prev_elem_id, elem.id, 
                    "Node {} is already connected to Elem {} at {:?} ({}); cannot connect to Elem {}",
                    self.id,
                    prev_elem_id,
                    address,
                    index_of_self,
                    elem.id, 

                );
            }
        } else {
            panic!("Elem {} is not connected to Node {}; cannot reciprocate connection!", elem.id, self.id);
        }
    }

    pub fn active_elems(&self) -> Option<[usize; 4]> {
        self.active_elems
    }

    pub(crate) fn reset_activation(&mut self) {
        self.active_elems = None;
    }

    pub(crate) fn set_activation(&mut self) -> bool {
        let mut combination_sets: [BTreeMap<[usize; 2], u8>; 4] = [
            BTreeMap::new(), // "edge" 0 : [0, 1] -> v level
            BTreeMap::new(), // "edge" 1 : [2, 3] -> v level
            BTreeMap::new(), // "edge" 2 : [0, 2] -> u level
            BTreeMap::new(), // "edge" 3 : [1, 3] -> u level
        ];

        for (edge_idx, elem_indices, check_index) in [
            (0, [0, 1], 1),
            (1, [2, 3], 1),
            (2, [0, 2], 0),
            (3, [1, 3], 0),
        ] {
            for (key_0, elem_id_0) in self.elems[elem_indices[0]].iter() {
                for (key_1, elem_id_1) in self.elems[elem_indices[1]].iter() {
                    if key_0[check_index] == key_1[check_index] {
                        combination_sets[edge_idx].insert([*elem_id_0, *elem_id_1], key_0[check_index]);
                    }
                }
            }
        }

        unimplemented!()
    }

    /// Produce a Json Object that describes this Node
    pub fn to_json(&self) -> JsonValue {
        object! {
            "id": self.id,
            "boundary": self.boundary,
            "point": self.coords,
            "elems": JsonValue::from(self.elems.iter().map(|elem_list| {
                elem_list.values().copied().collect()
            }).collect::<Vec<Vec<usize>>>())
        }
    }
}
