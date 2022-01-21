use json::JsonValue;

use super::ParaDir;
use std::fmt;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
/// Description of an Elem's h-Refinement levels in the u and v directions
pub struct HLevels {
    pub u: u8,
    pub v: u8,
}

impl HLevels {
    pub fn from(u: u8, v: u8) -> Self {
        Self { u, v }
    }

    pub fn refined(&self, refinement: HRef) -> Self {
        match refinement {
            HRef::T => Self::from(self.u + 1, self.v + 1),
            HRef::U(_) => Self::from(self.u + 1, self.v),
            HRef::V(_) => Self::from(self.u, self.v + 1),
        }
    }

    pub fn edge_ranking(self, edge_dir: ParaDir) -> [u8; 2] {
        match edge_dir {
            ParaDir::U => [self.v, self.u],
            ParaDir::V => [self.u, self.v],
        }
    }

    pub fn node_ranking(self) -> [u8; 2] {
        [self.u, self.v]
    }
}

impl Default for HLevels {
    fn default() -> Self {
        Self { u: 0, v: 0 }
    }
}

impl Into<JsonValue> for HLevels {
    fn into(self) -> JsonValue {
        object! {
            "u": self.u,
            "v": self.v,
        }
    }
}

/// Description of an h-Refinement
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HRef {
    /// isotropic
    T,
    /// anisotropic about the u-direction (with the option for the subsequent v-refinement of one child Elem)
    U(Option<usize>),
    /// anisotropic about the v-direction (with the option for the subsequent u-refinement of one child Elem)
    V(Option<usize>),
}

impl HRef {
    pub const fn u() -> Self {
        Self::U(None)
    }

    pub const fn v() -> Self {
        Self::V(None)
    }

    pub const fn u_extened(child_idx: u8) -> Result<Self, HRefError> {
        match child_idx {
            0 => Ok(Self::U(Some(0))),
            1 => Ok(Self::U(Some(1))),
            _ => Err(HRefError::BisectionIdxExceeded),
        }
    }

    pub const fn v_extened(child_idx: u8) -> Result<Self, HRefError> {
        match child_idx {
            0 => Ok(Self::V(Some(0))),
            1 => Ok(Self::V(Some(1))),
            _ => Err(HRefError::BisectionIdxExceeded),
        }
    }

    pub(crate) fn indices_and_ids(
        &self,
        id_counter: &mut usize,
    ) -> Box<dyn Iterator<Item = (usize, usize)> + '_> {
        let starting_id = *id_counter;
        match self {
            Self::T => {
                *id_counter += 4;
                Box::new((0..4).map(move |idx| (idx, starting_id + idx)))
            }
            Self::U(_) | Self::V(_) => {
                *id_counter += 2;
                Box::new((0..2).map(move |idx| (idx, starting_id + idx)))
            }
        }
    }

    pub(crate) fn loc(&self, idx: usize) -> HRefLoc {
        match self {
            Self::T => match idx {
                0 => HRefLoc::SW,
                1 => HRefLoc::SE,
                2 => HRefLoc::NW,
                3 => HRefLoc::NE,
                _ => unreachable!(),
            },
            Self::U(_) => match idx {
                0 => HRefLoc::E,
                1 => HRefLoc::W,
                _ => unreachable!(),
            },
            Self::V(_) => match idx {
                0 => HRefLoc::S,
                1 => HRefLoc::N,
                _ => unreachable!(),
            },
        }
    }
}

/// The location of a child `Elem` relative to its parent following an [HRef]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HRefLoc {
    /// South West: T(0)
    SW,
    /// South East: T(1)
    SE,
    /// North West: T(2)
    NW,
    /// North East: T(3)
    NE,
    /// West: U(0)
    W,
    /// East: U(1)
    E,
    /// South: V(0)
    S,
    /// North: V(1)
    N,
}

impl HRefLoc {
    /// Index of the child Elem relative to its parent according to the convention defined by the `Elem` Type
    pub fn index(&self) -> usize {
        match self {
            Self::SW => 0,
            Self::SE => 1,
            Self::NW => 2,
            Self::NE => 3,
            Self::W => 0,
            Self::E => 1,
            Self::S => 0,
            Self::N => 1,
        }
    }

    /// given an input range in the parents parametric space, what sup range does the child `Elem` occupy?
    pub fn sub_range(&self, [u_range, v_range]: [[f64; 2]; 2]) -> [[f64; 2]; 2] {
        [self.sub_u_range(u_range), self.sub_v_range(v_range)]
    }

    /// given an input range along the parent's `u` parametric axis, what sub-range does the child `Elem` occupy?
    pub fn sub_u_range(&self, [min_u, max_u]: [f64; 2]) -> [f64; 2] {
        let middle = (min_u + max_u) / 2.0;
        match self {
            Self::SW => [min_u, middle],
            Self::SE => [middle, max_u],
            Self::NW => [min_u, middle],
            Self::NE => [middle, max_u],
            Self::W => [min_u, middle],
            Self::E => [middle, max_u],
            Self::S => [min_u, max_u],
            Self::N => [min_u, max_u],
        }
    }

    /// given an input range along the parent's `v` parametric axis, what sub-range does the child `Elem` occupy?
    pub fn sub_v_range(&self, [min_v, max_v]: [f64; 2]) -> [f64; 2] {
        let middle = (min_v + max_v) / 2.0;
        match self {
            Self::SW => [min_v, middle],
            Self::SE => [min_v, middle],
            Self::NW => [middle, max_v],
            Self::NE => [middle, max_v],
            Self::W => [min_v, max_v],
            Self::E => [min_v, max_v],
            Self::S => [min_v, middle],
            Self::N => [middle, max_v],
        }
    }
}

#[derive(Debug)]
pub enum HRefError {
    MinEdgeLength(usize),
    ElemHasChildren(usize),
    EdgeHasChildren(usize),
    UninitializedElem(usize),
    ElemDoesntExist(usize),
    DoubleRefinement(usize),
    EdgeOnEqualPoints(usize),
    BisectionIdxExceeded,
}

impl fmt::Display for HRefError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::MinEdgeLength(edge_id) => write!(
                f,
                "h-refinement will result in Edge length below minimum; Cannot h-refine Edge {}!",
                edge_id
            ),
            Self::ElemHasChildren(elem_id) => {
                write!(f, "Elem {} already has children; Cannot h-refine!", elem_id)
            }
            Self::EdgeHasChildren(edge_id) => {
                write!(f, "Edge {} already has children; Cannot h-refine!", edge_id)
            }
            Self::UninitializedElem(elem_uninit_id) => write!(
                f,
                "ElemUninit {} was not fully initialized by the conclusion of h-refinement",
                elem_uninit_id
            ),
            Self::ElemDoesntExist(elem_id) => write!(
                f,
                "Elem {} does not exist; cannot apply h-Refinement!",
                elem_id
            ),
            Self::DoubleRefinement(elem_id) => write!(
                f,
                "Multiple h-refinements were specified for Elem {} in the same generation!",
                elem_id
            ),
            Self::EdgeOnEqualPoints(elem_id) => write!(
                f,
                "Attempt to generate a child-Edge between two identical points over Elem {}!",
                elem_id
            ),
            Self::BisectionIdxExceeded => write!(f, "Extended refinement index must be 0 or 1!"),
        }
    }
}
