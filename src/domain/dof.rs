//! Specification for a BasisFn defined on Some Elem.
//! Keeps track of the orders, direction, and associated DoF.
pub mod basis_spec;

use basis_spec::BSAddress;
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

    /// Get the list of addresses for the 1, 2 or 4 BasisSpecs associated with this DoF.
    pub fn get_basis_specs(&self) -> SmallVec<[BSAddress; 4]> {
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
        write!(f, "DoF {} \t {}", self.id, self.basis_specs)
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
            Self::ELEM(address) => write!(f, "ElemType[ {} ]", address),
            Self::EDGE([address_0, address_1]) => {
                write!(f, "EdgeType[ {} {} ]", address_0, address_1)
            }
            Self::NODE(bs_addresses) => {
                write!(f, "NodeType[")?;
                for bsa in bs_addresses.iter() {
                    write!(f, " {}", bsa)?;
                }
                write!(f, " ]")
            }
        }
    }
}
