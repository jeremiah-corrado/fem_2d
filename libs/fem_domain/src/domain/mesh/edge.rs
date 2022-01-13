use super::{ParaDir, h_refinement::{Bisection, HRefError}};
use std::collections::BTreeMap;

#[derive(Debug)]
pub struct Edge {
    pub id: usize,
    pub nodes: [usize; 2],
    pub boundary: bool,
    pub dir: ParaDir,
    children: Option<[usize; 2]>,
    parent: Option<(usize, Bisection)>,
    elems: [BTreeMap<[u8; 2], usize>; 2],
}

impl Edge {
    pub fn new(id: usize, nodes: [usize; 2], boundary: bool, dir: ParaDir) -> Self {
        Self {
            id,
            nodes,
            boundary,
            dir,
            children: None,
            parent: None,
            elems: [BTreeMap::new(), BTreeMap::new()],
        }
    }

    pub fn h_refine(&self, new_ids: [usize; 2], new_node_id: usize) -> Result<[Self; 2], HRefError> {
        Ok([Self::new(), Self::new()])
    }
}