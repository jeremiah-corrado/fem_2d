use super::MAX_POLYNOMIAL_ORDER;
use crate::domain::{mesh::space::ParaDir, dof::basis_spec::BasisDir};
use json::{object, JsonValue};
use std::{cmp::Ordering, fmt, ops::AddAssign};

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

    pub fn as_i8(&self) -> i8 {
        match self {
            Self::None => 0,
            Self::Increment(delta) => *delta as i8,
            Self::Decrement(delta) => -1 * (*delta as i8),
        }
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

    pub fn on_dir(dir: ParaDir, delta: i8) -> Self {
        match dir {
            ParaDir::U => Self::from(delta, 0),
            ParaDir::V => Self::from(0, delta),
        }
    }

    pub fn as_array(&self) -> [i8; 2] {
        [
            self.di.as_i8(),
            self.dj.as_i8(),
        ]
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
