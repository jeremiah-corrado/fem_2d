use super::ParaDir;
use std::fmt;

/// Description of an h-Refinement
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HRef {
    /// isotropic
    T,
    /// anisotropic about the u-direction (with the option for the subsequent v-refinement of one child Elem)
    U(Option<Bisection>),
    /// anisotropic about the v-direction (with the option for the subsequent u-refinement of one child Elem)
    V(Option<Bisection>),
}

impl HRef {
    pub fn indices_and_ids(
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

    pub fn location(&self, idx: usize) -> HRefLoc {
        match self {
            Self::T => HRefLoc::T(Quadrant::from(idx)),
            Self::U(_) | Self::V(_) => HRefLoc::UV(Bisection::from(idx)),
        }
    }
}

/// Quadrant of a child Elem following a T-Type h-Refinement (from the parent Elem's perspective)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Quadrant {
    /// south west
    SW,
    /// south east
    SE,
    /// north west
    NW,
    /// north east
    NE,
}

impl Quadrant {
    fn from(idx: usize) -> Self {
        match idx {
            0 => Self::SW,
            1 => Self::SE,
            2 => Self::NW,
            3 => Self::NE,
            _ => panic!("Quadrant index cannot exceed 3!"),
        }
    }

    pub fn index(&self) -> usize {
        match self {
            Self::SW => 0,
            Self::SE => 1,
            Self::NW => 2,
            Self::NE => 3,
        }
    }
}

/// Location of a child Elem following a U-Type or V-Type h-refinement (from the parent Elem's perspective).
/// Or the Location of a child Edge following an h-refinement.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Bisection {
    /// Bottom (V-type) or Left (U-type)
    BL,
    /// Top (V-Type) or Right (U-type)
    TR,
}

impl Bisection {
    fn from(idx: usize) -> Self {
        match idx {
            0 => Self::BL,
            1 => Self::TR,
            _ => panic!("Bisected Index must be 0 or 1!"),
        }
    }

    pub fn index(&self) -> usize {
        match self {
            Self::BL => 0,
            Self::TR => 1,
        }
    }
}

/// The location of an [Elem] relative to its parent following an h-refinement
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HRefLoc {
    T(Quadrant),
    UV(Bisection),
}

impl HRefLoc {
    pub fn index(&self) -> usize {
        match self {
            Self::T(quad) => quad.index(),
            Self::UV(bi) => bi.index(),
        }
    }
}

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
            ParaDir::U => [self.u, self.v],
            ParaDir::V => [self.v, self.u],
        }
    }
}

impl Default for HLevels {
    fn default() -> Self {
        Self { u: 0, v: 0 }
    }
}

#[derive(Debug)]
pub enum HRefError {
    MinEdgeLength(usize),
    ElemHasChildren(usize),
    EdgeHasChildren(usize),
    UninitializedElem(usize),
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
        }
    }
}
