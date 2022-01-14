use super::{
    h_refinement::{Bisection, HRefError},
    Elem, Node, ParaDir, MIN_EDGE_LENGTH,
};
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
}

impl Edge {
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
        }
    }

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

    pub fn h_refine(
        &mut self,
        new_ids: [usize; 2],
        new_node_id: usize,
    ) -> Result<[Self; 2], HRefError> {
        match self.children {
            Some(_) => Err(HRefError::EdgeHasChildren(self.id)),
            None => {
                let child_edge_length = self.length / 2.0;

                if child_edge_length < MIN_EDGE_LENGTH {
                    Err(HRefError::MinEdgeLength(self.id))
                } else {
                    self.children = Some(new_ids.clone());
                    Ok([
                        Self {
                            id: new_ids[0],
                            nodes: [self.nodes[0], new_node_id],
                            boundary: self.boundary,
                            dir: self.dir,
                            length: child_edge_length,
                            children: None,
                            parent: Some((self.id, Bisection::BL)),
                            elems: [BTreeMap::new(), BTreeMap::new()],
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
                        },
                    ])
                }
            }
        }
    }
}
