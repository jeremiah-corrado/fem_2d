mod basis_spec;

pub use basis_spec::{BSAddress, BasisDir, BasisLoc, BasisSpec};
use smallvec::SmallVec;
use std::fmt;

/// A single degree of freedom
pub struct DoF {
    pub id: usize,
    basis_specs: BasisSpecGroup,
}

impl DoF {
    pub fn new(id: usize, bs_addresses: SmallVec<[BSAddress; 4]>) -> Self {
        Self {
            id,
            basis_specs: match bs_addresses.len() {
                1 => BasisSpecGroup::ELEM(bs_addresses[0]),
                2 => BasisSpecGroup::EDGE([bs_addresses[0], bs_addresses[1]]),
                4 => BasisSpecGroup::NODE([
                    bs_addresses[0],
                    bs_addresses[1],
                    bs_addresses[2],
                    bs_addresses[3],
                ]),
                _ => panic!(
                    "BasisSpec groups must contain 1, 2, or 4 Addresses; cannot construct DoF {}!",
                    id
                ),
            },
        }
    }

    pub fn basis_spec_addresses(&self) -> SmallVec<[BSAddress; 4]> {
        match self.basis_specs {
            BasisSpecGroup::ELEM(elem_bs_address) => smallvec![elem_bs_address],
            BasisSpecGroup::EDGE(edge_bs_addresses) => {
                smallvec![edge_bs_addresses[0], edge_bs_addresses[1]]
            }
            BasisSpecGroup::NODE(node_bs_addresses) => smallvec![
                node_bs_addresses[0],
                node_bs_addresses[1],
                node_bs_addresses[2],
                node_bs_addresses[3]
            ],
        }
    }
}

impl fmt::Display for DoF {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} \t {}", self.id, self.basis_specs,)
    }
}

enum BasisSpecGroup {
    ELEM(BSAddress),
    EDGE([BSAddress; 2]),
    NODE([BSAddress; 4]),
}

impl fmt::Display for BasisSpecGroup {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::ELEM(address) => write!(f, "Elem[{}, ({})]", address.elem_id, address.bs_id),
            Self::EDGE([address_0, address_1]) => write!(
                f,
                "Edge[elem_ids: {}, {}, bs_ids: ({}, {})]",
                address_0.elem_id, address_1.elem_id, address_0.bs_id, address_1.bs_id
            ),
            Self::NODE([address_0, address_1, address_2, address_3]) => write!(
                f,
                "Elem[elem_ids: {}, {}, {}, {}, bs_ids: ({}, {}, {}, {})]",
                address_0.elem_id,
                address_1.elem_id,
                address_2.elem_id,
                address_3.elem_id,
                address_0.bs_id,
                address_1.bs_id,
                address_2.bs_id,
                address_3.bs_id,
            ),
        }
    }
}
