mod basis_spec;

pub use basis_spec::{BasisDir, BasisLoc, BasisSpec};
use smallvec::SmallVec;
use std::fmt;

/// A single degree of freedom
pub struct DoF {
    pub id: usize,
    basis_specs: BasisSpecGroup,
}

impl DoF {
    pub fn new(id: usize, basis_specs: &[&BasisSpec]) -> Self {
        Self {
            id,
            basis_specs: match basis_specs.len() {
                1 => BasisSpecGroup::ELEM(basis_specs[0].id),
                2 => BasisSpecGroup::EDGE([basis_specs[0].id, basis_specs[1].id]),
                4 => BasisSpecGroup::NODE([
                    basis_specs[0].id,
                    basis_specs[1].id,
                    basis_specs[2].id,
                    basis_specs[3].id,
                ]),
                _ => panic!(
                    "BasisSpec groups must contain 1, 2, or 4 BasisSpecs; cannot construct DoF {}!",
                    id
                ),
            },
        }
    }

    pub fn basis_spec_ids(&self) -> SmallVec<[usize; 4]> {
        match self.basis_specs {
            BasisSpecGroup::ELEM(elem_bs_id) => smallvec![elem_bs_id],
            BasisSpecGroup::EDGE(edge_bs_ids) => smallvec![edge_bs_ids[0], edge_bs_ids[1]],
            BasisSpecGroup::NODE(node_bs_ids) => smallvec![
                node_bs_ids[0],
                node_bs_ids[1],
                node_bs_ids[2],
                node_bs_ids[3]
            ],
        }
    }
}

enum BasisSpecGroup {
    ELEM(usize),
    EDGE([usize; 2]),
    NODE([usize; 4]),
}
