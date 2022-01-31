use super::MAX_POLYNOMIAL_ORDER;
use crate::domain::dof::basis_spec::BasisDir;
use json::{object, JsonValue};
use std::fmt;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PolyOrders {
    /// Maximum u-directed polynomial expansion order
    pub ni: u8,
    /// Maximum v-directed polynomial expansion order
    pub nj: u8,
}

impl PolyOrders {
    pub const fn from(i: u8, j: u8) -> Self {
        Self { ni: i, nj: j }
    }

    /// Update the u- and v-directed expansion orders according to a [PRef]
    pub fn refine(&mut self, refinement: PRef) -> Result<(), PRefError> {
        self.ni = refinement.refine_i(self.ni)?;
        self.nj = refinement.refine_j(self.nj)?;

        Ok(())
    }

    pub fn set(&mut self, [ni, nj]: [u8; 2]) -> Result<(), PRefError> {
        if ni > MAX_POLYNOMIAL_ORDER || nj > MAX_POLYNOMIAL_ORDER {
            return Err(PRefError::ExceededMaxExpansion);
        }

        self.ni = ni;
        self.nj = nj;

        Ok(())
    }

    /// Get the permutations of [i, j] for the u- or v-directed shape functions.
    ///
    /// * For u-directed: i ∈ [0, Ni) and j ∈ [0, Nj]
    /// * For v-directed: i ∈ [0, Ni] and j ∈ [0, Nj)
    /// * For w-directed: i ∈ [0, Ni] and j ∈ [0, Nj]
    pub fn permutations(&self, dir: BasisDir) -> Box<dyn Iterator<Item = [u8; 2]> + '_> {
        match dir {
            BasisDir::U => Box::new(
                (0..self.ni)
                    .flat_map(move |i_order| (0..=self.nj).map(move |j_order| [i_order, j_order])),
            ),
            BasisDir::V => Box::new(
                (0..=self.ni)
                    .flat_map(move |i_order| (0..self.nj).map(move |j_order| [i_order, j_order])),
            ),
            BasisDir::W => Box::new(
                (0..=self.ni)
                    .flat_map(move |i_order| (0..=self.nj).map(move |j_order| [i_order, j_order])),
            ),
        }
    }

    pub fn max_with(&self, orders: [u8; 2]) -> [u8; 2] {
        [
            std::cmp::max(self.ni, orders[0]),
            std::cmp::max(self.nj, orders[1]),
        ]
    }
}

impl Default for PolyOrders {
    fn default() -> Self {
        Self { ni: 1, nj: 1 }
    }
}

#[cfg(feature = "json_export")]
impl From<PolyOrders> for JsonValue {
    fn from(orders: PolyOrders) -> Self {
        object! {
            "u": orders.ni,
            "v": orders.nj,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PRefInt {
    Increment(u8),
    Decrement(u8),
    None,
}

impl PRefInt {
    fn refine(&self, n: u8) -> Result<u8, PRefError> {
        match self {
            Self::Increment(delta) => {
                if n + *delta > MAX_POLYNOMIAL_ORDER {
                    Err(PRefError::ExceededMaxExpansion)
                } else {
                    Ok(n + *delta)
                }
            }
            Self::Decrement(delta) => {
                if *delta >= n {
                    Err(PRefError::NegExpansion)
                } else {
                    Ok(n - *delta)
                }
            }
            Self::None => Ok(n),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
/// Description of a p-Refinement
pub struct PRef {
    di: PRefInt,
    dj: PRefInt,
}

impl PRef {
    pub const fn from(i: i8, j: i8) -> Self {
        Self {
            di: match i {
                0 => PRefInt::None,
                d if d > 0 => PRefInt::Increment(d as u8),
                d if d < 0 => PRefInt::Decrement(-d as u8),
                _ => unreachable!(),
            },
            dj: match j {
                0 => PRefInt::None,
                d if d > 0 => PRefInt::Increment(d as u8),
                d if d < 0 => PRefInt::Decrement(-d as u8),
                _ => unreachable!(),
            },
        }
    }

    fn refine_i(&self, i_current: u8) -> Result<u8, PRefError> {
        self.di.refine(i_current)
    }

    fn refine_j(&self, j_current: u8) -> Result<u8, PRefError> {
        self.dj.refine(j_current)
    }
}

#[derive(Debug)]
pub enum PRefError {
    NegExpansion,
    ExceededMaxExpansion,
    ElemDoesntExist(usize),
    DoubleRefinement(usize),
}

impl fmt::Display for PRefError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::NegExpansion => write!(f, "Negative p-Refinement will result in 0 or negative expansion; Cannot p-Refine!"),
            Self::ExceededMaxExpansion => write!(f, "Positive p-Refinement will result in expansion order over maximum; Cannot p-Refine!"),
            Self::ElemDoesntExist(elem_id) => write!(f, "Elem {} does not exist; Cannot apply p-Refinement!", elem_id),
            Self::DoubleRefinement(elem_id) => write!(f, "Multiple p-refinements were specified for Elem {}; Cannot apply p-Refinements", elem_id),
        }
    }
}
