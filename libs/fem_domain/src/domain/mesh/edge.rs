use super::{
    h_refinement::{Bisection, HRefError},
    Node, ParaDir, MIN_EDGE_LENGTH,
};
use std::collections::BTreeMap;

/// A flat line in parametric space defined as the line between two Points in real space.
#[derive(Debug)]
pub struct Edge {
    pub id: usize,
    pub nodes: [usize; 2],
    pub boundary: bool,
    pub dir: ParaDir,
    pub length: f64,
    children: Option<[usize; 2]>,
    parent: Option<(usize, Bisection)>,
    elems: [BTreeMap<[u8; 2], usize>; 2],
}

impl Edge {
    pub fn new(id: usize, nodes: [&Node; 2], boundary: bool) -> Self {
        let dir = nodes[0].coords.orientation_with(&nodes[1].coords).expect("Nodes must share either an x or y coordinate; Cannot construct Edge!");

        Self {
            id,
            nodes: [nodes[0].id, nodes[1].id],
            boundary,
            dir,
            length: match dir {
                ParaDir::U => nodes[1].coords.x - nodes[0].coords.x,
                ParaDir::V => nodes[1].coords.y - nodes[0].coords.y,
            },
            children: None,
            parent: None,
            elems: [BTreeMap::new(), BTreeMap::new()],
        }
    }

    pub fn h_refine(
        &mut self,
        new_ids: [usize; 2],
        new_node_id: usize,
    ) -> Result<[Self; 2], HRefError> {
        match self.children {
            Some(_) => Err(HRefError::EdgeHasChildren(self.id)),
            None => {
                let length_child_edges = self.length / 2.0;

                if length_child_edges < MIN_EDGE_LENGTH {
                    Err(HRefError::MinEdgeLength(self.id))
                } else {
                    self.children = Some(new_ids.clone());
                    Ok([
                        Self {
                            id: new_ids[0],
                            nodes: [self.nodes[0], new_node_id],
                            boundary: self.boundary,
                            dir: self.dir,
                            length: length_child_edges,
                            children: None,
                            parent: Some((self.id, Bisection::BL)),
                            elems: [BTreeMap::new(), BTreeMap::new()],
                        },
                        Self {
                            id: new_ids[1],
                            nodes: [new_node_id, self.nodes[1]],
                            boundary: self.boundary,
                            dir: self.dir,
                            length: length_child_edges,
                            children: None,
                            parent: Some((self.id, Bisection::TR)),
                            elems: [BTreeMap::new(), BTreeMap::new()],
                        },
                    ])
                }
            }
        }
    }
}
