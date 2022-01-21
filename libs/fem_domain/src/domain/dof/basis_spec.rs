use crate::domain::mesh::Elem;
use std::fmt;

#[derive(Clone, Debug)]
pub struct BasisSpec {
    pub id: usize,
    pub i: u8,
    pub j: u8,
    pub dir: BasisDir,
    pub elem_id: usize,
    pub dof_id: Option<usize>,
    pub(crate) loc: BasisLoc,
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
            dof_id: None,
            elem_id: elem.id,
            loc,
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
            (BasisLoc::NODE(_, node_id_0), BasisLoc::NODE(_, node_id_1)) => {
                assert_eq!(
                    node_id_0, node_id_1,
                    "Cannot attempt to match Node-Type BasisSpecs associated with different Nodes!"
                );
                unimplemented!()
            }
            (_, _) => false,
        }
    }

    pub fn assign_dof_id(&mut self, dof_id: usize) {
        assert!(
            self.dof_id.is_none(),
            "BasisSpec {} is already connected to DoF {}; cannot connect to {}!",
            self.id,
            self.dof_id.unwrap(),
            dof_id
        );
        self.dof_id = Some(dof_id);
    }

    pub fn update_ids(&mut self, new_id: usize, dof_id: usize) {
        self.id = new_id;
        self.dof_id = Some(dof_id);
    }

    pub fn integration_data(&self) -> ([usize; 2], BasisDir, usize) {
        debug_assert!(
            self.dof_id.is_some(),
            "Cannot get integration data for unmatched BasisSpec: {}",
            self.id
        );
        (
            [self.i as usize, self.j as usize],
            self.dir,
            self.dof_id.unwrap(),
        )
    }
}

impl fmt::Display for BasisSpec {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "id:{} \t [{}, {}] - {} dir \t elem id: {} \t dof_id: {}",
            self.id,
            self.i,
            self.j,
            self.dir,
            self.elem_id,
            match self.dof_id {
                Some(id) => id.to_string(),
                None => String::from("_"),
            }
        )
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BasisDir {
    U,
    V,
    W,
}

impl fmt::Display for BasisDir {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            BasisDir::U => write!(f, "U"),
            BasisDir::V => write!(f, "V"),
            BasisDir::W => write!(f, "W"),
        }
    }
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

#[derive(Clone, Copy, Debug)]
pub struct BSAddress {
    pub elem_id: usize,
    pub bs_id: usize,
}

impl BSAddress {
    pub fn new(elem_id: usize, bs_id: usize) -> Self {
        Self { elem_id, bs_id }
    }
}
