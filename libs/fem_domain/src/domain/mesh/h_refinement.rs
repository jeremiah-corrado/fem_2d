

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
        Self {
            u, v
        }
    }

    pub fn refined(&self, refinement: HRef) -> Self {
        match refinement {
            HRef::T => Self::from(self.u + 1, self.v + 1),
            HRef::U(_) => Self::from(self.u + 1, self.v),
            HRef::V(_) => Self::from(self.u, self.v + 1),
        }
    }
}

impl Default for HLevels {
    fn default() -> Self {
        Self {
            u: 1,
            v: 1,
        }
    }
}