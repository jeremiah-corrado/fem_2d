use crate::domain::mesh::Elem;

#[derive(Clone, Debug)]
pub struct BasisSpec {
    pub id: usize,
    pub i: u8,
    pub j: u8,
    pub dir: BasisDir,
    pub(crate) loc: BasisLoc,
    pub elem_id: usize,
}

impl BasisSpec {
    pub fn new(id: usize, [i, j]: [u8; 2], dir: BasisDir, elem: &Elem) -> Self {
        let loc = match (i, j, dir) {
            (_, 0..=1, BasisDir::U) => BasisLoc::edge_bs(elem, j),
            (0..=1, _, BasisDir::V) => BasisLoc::edge_bs(elem, i + 2),
            (_, _, BasisDir::W) => match (i < 2, j < 2) {
                (true, false) => BasisLoc::edge_bs(elem, i + 2),
                (false, true) => BasisLoc::edge_bs(elem, j),
                (true, true) => BasisLoc::node_bs(elem, i + 2 * j),
                (false, false) => BasisLoc::ELEM,
            },
            (2..=u8::MAX, 2..=u8::MAX, _) => BasisLoc::ELEM,
            (_, _, _) => BasisLoc::ELEM,
        };

        Self {
            id,
            i,
            j,
            dir,
            loc,
            elem_id: elem.id,
        }
    }

    /// Assumes the two basis specs are associated with the same geometric components
    pub fn matches_with(&self, other: &Self) -> bool {
        match (self.loc, other.loc) {
            (BasisLoc::EDGE(idx_0, edge_id_0), BasisLoc::EDGE(idx_1, edge_id_1)) => {
                assert_eq!(
                    edge_id_0, edge_id_1,
                    "Cannot attempt to match Edge-Type BasisSpecs associated with different Edges!"
                );
                match (self.dir, other.dir) {
                    (BasisDir::U, BasisDir::U) => {
                        self.i == other.i && self.j + other.j == 1 && idx_0 + idx_1 == 1
                    }
                    (BasisDir::V, BasisDir::V) => {
                        self.j == other.j && self.i + other.i == 1 && idx_0 + idx_1 == 5
                    }
                    (BasisDir::W, BasisDir::W) => {
                        if self.i >= 2 && other.i >= 2 {
                            self.i == other.i && self.j + other.j == 1 && idx_0 + idx_1 == 1
                        } else if self.j >= 2 && other.j >= 2 {
                            self.j == other.j && self.i + other.i == 1 && idx_0 + idx_1 == 5
                        } else {
                            false
                        }
                    }
                    (_, _) => false,
                }
            }
            (BasisLoc::NODE(idx_0, node_id_0), BasisLoc::NODE(idx_1, node_id_1)) => {
                assert_eq!(
                    node_id_0, node_id_1,
                    "Cannot attempt to match Node-Type BasisSpecs associated with different Nodes!"
                );
                unimplemented!()
            }
            (_, _) => false,
        }
    }

    pub fn update_id(&mut self, new_id: usize) {
        self.id = new_id;
    }
}

#[derive(Clone, Copy, Debug)]
pub enum BasisDir {
    U,
    V,
    W,
}

#[derive(Clone, Copy, Debug)]
pub enum BasisLoc {
    ELEM,
    EDGE(u8, usize),
    NODE(u8, usize),
}

impl BasisLoc {
    pub fn edge_bs(elem: &Elem, idx: u8) -> Self {
        Self::EDGE(idx, elem.edges[idx as usize])
    }

    pub fn node_bs(elem: &Elem, idx: u8) -> Self {
        Self::NODE(idx, elem.nodes[idx as usize])
    }
}

pub struct BSAddress {
    elem_id: usize,
    bs_id: usize,
}

impl BSAddress {
    pub fn new(elem_id: usize, bs_id: usize) -> Self {
        Self {
            elem_id, bs_id,
        }
    }
}
