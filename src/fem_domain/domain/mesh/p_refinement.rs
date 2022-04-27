use super::MeshAccessError;
use super::MAX_POLYNOMIAL_ORDER;
use crate::fem_domain::domain::{dof::basis_spec::BasisDir, mesh::space::ParaDir};
use json::{object, JsonValue};
use std::{cmp::Ordering, fmt, ops::AddAssign};

/// The record of polynomial expansion orders associated with an Elem
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PolyOrders {
    /// Maximum u-directed polynomial expansion order
    pub ni: u8,
    /// Maximum v-directed polynomial expansion order
    pub nj: u8,
}

impl PolyOrders {
    /// Update the u- and v-directed expansion orders according to a [PRef]
    ///
    /// Return an `Err` if the refinement causes the expansion orders to fall outside of the valid range
    pub fn refine(&mut self, refinement: PRef) -> Result<(), PRefError> {
        self.ni = refinement.refine_i(self.ni)?;
        self.nj = refinement.refine_j(self.nj)?;

        Ok(())
    }

    /// Directly update the u- and v-directed expansion orders to the given values
    ///
    /// Return an `Err` if the given values are out of the valid range
    pub fn set(&mut self, [ni, nj]: [u8; 2]) -> Result<(), PRefError> {
        if ni > MAX_POLYNOMIAL_ORDER || nj > MAX_POLYNOMIAL_ORDER {
            return Err(PRefError::ExceededMaxExpansion);
        }
        if ni < 1 || nj < 1 {
            return Err(PRefError::NegExpansion);
        }

        self.ni = ni;
        self.nj = nj;

        Ok(())
    }

    /// Get the permutations of [i, j] for the u-, v- or w-directed basis functions
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

    /// The maximum orders from self and the given orders
    pub fn max_with(&self, orders: [u8; 2]) -> [u8; 2] {
        [
            std::cmp::max(self.ni, orders[0]),
            std::cmp::max(self.nj, orders[1]),
        ]
    }

    /// Return the expansion orders as an array of `usize`
    pub fn as_array(&self) -> [usize; 2] {
        [self.ni as usize, self.nj as usize]
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

// the internal p-refinement type
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

    pub const fn from_i8(delta: i8) -> Self {
        match delta {
            0 => PRefInt::None,
            d if d > 0 => PRefInt::Increment(d as u8),
            d if d < 0 => PRefInt::Decrement(-d as u8),
            _ => unreachable!(),
        }
    }

    pub fn as_i8(&self) -> i8 {
        match self {
            Self::None => 0,
            Self::Increment(delta) => *delta as i8,
            Self::Decrement(delta) => -1 * (*delta as i8),
        }
    }

    pub fn constrain(&mut self, bounds: [i8; 2]) {
        let si = self.as_i8();

        if si < bounds[0] {
            *self = Self::from_i8(bounds[0]);
        } else if si > bounds[1] {
            *self = Self::from_i8(bounds[1]);
        }
    }

    pub fn within(&self, bounds: [i8; 2]) -> bool {
        let si = self.as_i8();
        si >= bounds[0] && si <= bounds[1]
    }
}

impl AddAssign for PRefInt {
    fn add_assign(&mut self, rhs: Self) {
        let sum = match *self {
            Self::Increment(s_delta) => match rhs {
                Self::Increment(r_delta) => Self::Increment(s_delta + r_delta),
                Self::Decrement(r_delta) => match s_delta.cmp(&r_delta) {
                    Ordering::Equal => Self::None,
                    Ordering::Greater => Self::Increment(s_delta - r_delta),
                    Ordering::Less => Self::Decrement(r_delta - s_delta),
                },
                Self::None => self.clone(),
            },
            Self::Decrement(s_delta) => match rhs {
                Self::Decrement(r_delta) => Self::Decrement(s_delta + r_delta),
                Self::Increment(r_delta) => match s_delta.cmp(&r_delta) {
                    Ordering::Equal => Self::None,
                    Ordering::Greater => Self::Decrement(s_delta - r_delta),
                    Ordering::Less => Self::Increment(r_delta - s_delta),
                },
                Self::None => self.clone(),
            },
            Self::None => rhs.clone(),
        };

        *self = sum;
    }
}

impl fmt::Display for PRefInt {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            PRefInt::None => write!(f, "0"),
            PRefInt::Decrement(delta) => write!(f, "-{}", delta),
            PRefInt::Increment(delta) => write!(f, "+{}", delta),
        }
    }
}

/// The p-Refinement Type
///
/// p-refinements are used to modify the expansion orders associated with an Elem
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PRef {
    di: PRefInt,
    dj: PRefInt,
}

impl PRef {
    /// Create a new p-refinement with the given deltas for the i and j expansion orders
    pub const fn from(delta_i: i8, delta_j: i8) -> Self {
        Self {
            di: PRefInt::from_i8(delta_i),
            dj: PRefInt::from_i8(delta_j),
        }
    }

    /// Create a new p-refinement with the given deltas in a specific parametric direction
    pub fn on_dir(dir: ParaDir, delta: i8) -> Self {
        match dir {
            ParaDir::U => Self::from(delta, 0),
            ParaDir::V => Self::from(0, delta),
        }
    }

    /// Get the p-refinement deltas as an array of `i8`s
    pub fn as_array(&self) -> [i8; 2] {
        [self.di.as_i8(), self.dj.as_i8()]
    }

    /// Constrict the p-refinement to sit within a given range
    pub fn constrained_to(mut self, [u_bounds, v_bounds]: [[i8; 2]; 2]) -> Self {
        self.di.constrain(u_bounds);
        self.dj.constrain(v_bounds);

        self
    }

    /// Check if the refinement falls within the given bounds
    pub fn falls_within(&self, [u_bounds, v_bounds]: [[i8; 2]; 2]) -> bool {
        self.di.within(u_bounds) && self.dj.within(v_bounds)
    }

    fn refine_i(&self, i_current: u8) -> Result<u8, PRefError> {
        self.di.refine(i_current)
    }

    fn refine_j(&self, j_current: u8) -> Result<u8, PRefError> {
        self.dj.refine(j_current)
    }
}

impl AddAssign for PRef {
    fn add_assign(&mut self, rhs: Self) {
        self.di += rhs.di;
        self.dj += rhs.dj;
    }
}

impl fmt::Display for PRef {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "PRef (i: {}, j: {})", self.di, self.dj)
    }
}

/// The Error Type for invalid p-refinements
#[derive(Debug)]
pub enum PRefError {
    // Errors caused by internal problems with the Mesh (Should never happen)
    NegExpansion,
    ExceededMaxExpansion,
    // Public Errors
    DuplicateElemIds,
    ElemDoesNotExist(usize),
    RefinementOutOfBounds(usize),
}

impl std::error::Error for PRefError {}

impl fmt::Display for PRefError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::NegExpansion => write!(
                f,
                "Negative p-Refinement results in 0 or negative expansion order; Cannot p-Refine!"
            ),
            Self::ExceededMaxExpansion => write!(
                f,
                "Positive p-Refinement results in expansion order over maximum; Cannot p-Refine!"
            ),
            Self::DuplicateElemIds => write!(
                f,
                "Duplicate element ids in p-Refinement; Cannot apply p-Refinement!"
            ),
            Self::ElemDoesNotExist(elem_id) => write!(
                f,
                "Elem {} does not exist; Cannot apply p-Refinement!",
                elem_id
            ),
            Self::RefinementOutOfBounds(elem_id) => write!(
                f,
                "Refinement out of bounds for Elem {}; Cannot apply p-Refinement!",
                elem_id
            ),
        }
    }
}

impl From<MeshAccessError> for PRefError {
    fn from(err: MeshAccessError) -> Self {
        match err {
            MeshAccessError::ElemDoesNotExist(elem_id) => Self::ElemDoesNotExist(elem_id),
            _ => unreachable!(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn p_ref_constraints() {
        let pr = PRef::from(4, 6);
        let prc = pr.constrained_to([[0, 1], [-2, 5]]).as_array();
        assert_eq!(prc[0], 1);
        assert_eq!(prc[1], 5);

        let pr_neg = PRef::from(-3, -2);
        let prc_neg = pr_neg.constrained_to([[-2, 4], [0, 3]]).as_array();
        assert_eq!(prc_neg[0], -2);
        assert_eq!(prc_neg[1], 0);
    }

    #[test]
    fn p_ref_falls_within() {
        let pr = PRef::from(0, 10);
        assert!(pr.falls_within([[0, 1], [-2, 10]]));
        assert!(!pr.falls_within([[-2, -1], [3, 10]]));
        assert!(!pr.falls_within([[-1, 1], [11, 14]]));
    }
}
