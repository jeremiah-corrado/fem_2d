/// Structures to specify Basis Functions over `Elem`s
pub mod basis_spec;

use basis_spec::BSAddress;
use smallvec::{smallvec, SmallVec};
use std::fmt;

/// A single Degree-of-Freedom in the connected system
///
/// A DoF specifies a basis function or group of basis functions who must share the same solution values
pub struct DoF {
    pub id: usize,
    basis_specs: BasisSpecGroup,
}

impl DoF {
    /// Create a new DoF with the given id and the associated basis functions
    pub fn new(id: usize, bs_addresses: SmallVec<[BSAddress; 4]>) -> Self {
        Self {
            id,
            basis_specs: match bs_addresses.len() {
                1 => BasisSpecGroup::ElemGroup(bs_addresses[0]),
                2 => BasisSpecGroup::EdgeGroup([bs_addresses[0], bs_addresses[1]]),
                4 => BasisSpecGroup::NodeGroup([
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
            BasisSpecGroup::ElemGroup(elem_bs_address) => smallvec![elem_bs_address],
            BasisSpecGroup::EdgeGroup(edge_bs_addresses) => {
                smallvec![edge_bs_addresses[0], edge_bs_addresses[1]]
            }
            BasisSpecGroup::NodeGroup(node_bs_addresses) => smallvec![
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
    ElemGroup(BSAddress),
    EdgeGroup([BSAddress; 2]),
    NodeGroup([BSAddress; 4]),
}

impl fmt::Display for BasisSpecGroup {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::ElemGroup(address) => write!(f, "ElemType[ {} ]", address),
            Self::EdgeGroup([address_0, address_1]) => {
                write!(f, "EdgeType[ {} {} ]", address_0, address_1)
            }
            Self::NodeGroup(bs_addresses) => {
                write!(f, "NodeType[")?;
                for bsa in bs_addresses.iter() {
                    write!(f, " {}", bsa)?;
                }
                write!(f, " ]")
            }
        }
    }
}
