use super::super::ContinuityCondition;
use crate::fem_domain::domain::mesh::elem::Elem;
use std::fmt;

/// A description of a Basis Function with expansion orders, vectorial direction, and an associated geometric component
///
/// This is used to describe a basis function and facilitate the construction of degrees-of-freedom. The information in this struct is also used to construct the actual `BasisFn` during numerical integration
#[derive(Clone, Debug)]
pub struct BasisSpec {
    /// ID for this Basis Function in the global list of `BasisSpec`s
    pub id: usize,
    /// Expansion order in the u-direction
    pub i: u8,
    /// Expansion order in the v-direction
    pub j: u8,
    /// Vectorial Direction of the Basis Function
    pub dir: BasisDir,
    /// ID of the element over which the Basis Function is Defined
    pub elem_id: usize,
    /// Index of this `BasisSpec` in the `Elem`s list of `BasisSpec`s
    pub elem_idx: Option<usize>,
    /// ID of the Degree of Freedom that this Basis Function is associated with
    pub dof_id: Option<usize>,
    /// The local Geometric Component this Basis Function is associated with
    pub loc: BasisLoc,
}

impl BasisSpec {
    /// Create a new `BasisSpec` over a given [Elem]
    ///
    /// # Arguments
    /// * `id`: Location of this `BasisSpec` in the global list of `BasisSpec`s
    /// * `i`: Expansion order in the u-direction
    /// * `j`: Expansion order in the v-direction
    /// * `dir`: Vectorial Direction of the Basis Function
    /// * `elem`: The [Elem] over which this `BasisSpec` is defined
    /// * `cc`: The [ContinuityCondition] that this `BasisSpec` should conform to
    ///
    /// Upon construction the `elem_idx`, and `dof_id` are set to `None`. These are to be defined during DoF construction.
    ///
    pub fn new(
        id: usize,
        [i, j]: [u8; 2],
        dir: BasisDir,
        elem: &Elem,
        cc: ContinuityCondition,
    ) -> Self {
        // select the associated geometric component based on the continuity condition, vectorial direction, and expansion orders
        let loc = match cc {
            ContinuityCondition::HCurl => match (i, j, dir) {
                (2..=u8::MAX, 2..=u8::MAX, _) => BasisLoc::ElemBs,
                (_, 0..=1, BasisDir::U) => BasisLoc::edge_bs(elem, j),
                (0..=1, _, BasisDir::V) => BasisLoc::edge_bs(elem, i + 2),
                (_, _, BasisDir::W) => match (i < 2, j < 2) {
                    (true, false) => BasisLoc::edge_bs(elem, i + 2),
                    (false, true) => BasisLoc::edge_bs(elem, j),
                    (true, true) => BasisLoc::node_bs(elem, i + 2 * j),
                    (false, false) => BasisLoc::ElemBs,
                },
                (_, _, _) => BasisLoc::ElemBs,
            },
            _ => unimplemented!(),
        };

        Self {
            id,
            i,
            j,
            dir,
            elem_id: elem.id,
            elem_idx: None,
            dof_id: None,
            loc,
        }
    }

    /// Checks whether two edge-type BasisSpecs are compatible for matching along their shared edge
    ///
    /// panics if the basis specs are not edge-type or if they are not attached to the same edge
    pub fn matches_with_edge(&self, other: &Self) -> bool {
        match (self.loc, other.loc) {
            (BasisLoc::EdgeBs(idx_0, edge_id_0), BasisLoc::EdgeBs(idx_1, edge_id_1)) => {
                assert_eq!(
                    edge_id_0, edge_id_1,
                    "Cannot attempt to match Edge-Type BasisSpecs associated with different Edges!"
                );
                match (self.dir, other.dir) {
                    (BasisDir::U, BasisDir::U) => {
                        self.i == other.i && self.j + other.j == 1 && idx_0 + idx_1 == 1
                    }
                    (BasisDir::V, BasisDir::V) => {
                        self.j == other.j && self.i + other.i == 1 && idx_0 + idx_1 == 5
                    }
                    (BasisDir::W, BasisDir::W) => {
                        if self.i >= 2 && other.i >= 2 {
                            self.i == other.i && self.j + other.j == 1 && idx_0 + idx_1 == 1
                        } else if self.j >= 2 && other.j >= 2 {
                            self.j == other.j && self.i + other.i == 1 && idx_0 + idx_1 == 5
                        } else {
                            false
                        }
                    }
                    (_, _) => false,
                }
            }
            (_, _) => {
                panic!("Cannot test for edge-type BasisSpec match with non-edge-type BasisSpecs!")
            }
        }
    }

    // TODO: Implement Node-type matching

    /// Set the `dof_id` and `elem_idx` (the position of this BasisSpec in it's Elem's Vec<BasisSpec>)
    ///
    /// Panics if these indices have already been set
    pub fn set_dof_and_idx(&mut self, dof_id: usize, elem_idx: usize) {
        assert!(
            self.dof_id.is_none(),
            "DoF ID was already set to {} on BasisSpec {}; cannot set to {}!",
            self.dof_id.unwrap(),
            self.id,
            dof_id
        );
        assert!(
            self.elem_idx.is_none(),
            "Elem index was already set to {} on BasisSpec {}; cannot set to {}!",
            self.elem_idx.unwrap(),
            self.id,
            elem_idx
        );

        self.dof_id = Some(dof_id);
        self.elem_idx = Some(elem_idx);
    }

    /// Get a tuple of information needed to compute an integral with this BasisSpec
    pub fn integration_data(&self) -> ([usize; 2], BasisDir, usize) {
        debug_assert!(
            self.dof_id.is_some(),
            "Cannot get integration data for unmatched BasisSpec: {}",
            self.id
        );
        (
            [self.i as usize, self.j as usize],
            self.dir,
            self.dof_id.unwrap(),
        )
    }
}

impl fmt::Display for BasisSpec {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "id:{} \t [{}, {}] - {} dir \t elem id: {} \t dof_id: {}",
            self.id,
            self.i,
            self.j,
            self.dir,
            self.elem_id,
            match self.dof_id {
                Some(id) => id.to_string(),
                None => String::from("_"),
            }
        )
    }
}

impl PartialEq for BasisSpec {
    fn eq(&self, other: &Self) -> bool {
        self.i == other.i && self.j == other.j && self.loc == other.loc && self.dir == other.dir
    }
}

/// Vectorial Direction of a Basis Function in Parametric Space (or which component of the vectorial solution it belongs to)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BasisDir {
    /// X-Directed Solution
    U,
    /// Y-Directed Solution
    V,
    /// Z-Directed Solution
    W,
}

impl fmt::Display for BasisDir {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            BasisDir::U => write!(f, "U"),
            BasisDir::V => write!(f, "V"),
            BasisDir::W => write!(f, "W"),
        }
    }
}

/// A Description of a Basis Functions geometric associations within a given [Elem]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BasisLoc {
    /// The BasisSpec's location is the same as its `elem_id`
    ElemBs,
    /// (`edge_index` [0, 1, 2, or 3], `edge_id`)
    EdgeBs(u8, usize),
    /// (`node_index` [0, 1, 2, or 3], `node_id`)
    NodeBs(u8, usize),
}

impl BasisLoc {
    /// Generate an Edge-Type `BasisLoc`
    pub fn edge_bs(elem: &Elem, idx: u8) -> Self {
        Self::EdgeBs(idx, elem.edges[idx as usize])
    }

    /// Generate an Node-Type `BasisLoc`
    pub fn node_bs(elem: &Elem, idx: u8) -> Self {
        Self::NodeBs(idx, elem.nodes[idx as usize])
    }
}

#[derive(Clone, Copy, Debug)]
/// A BasisSpec's index in an Elem's `Vec<BasisSpec>`
pub struct BSAddress {
    /// Index of the `Elem` that the basis spec belongs to
    pub elem_id: usize,
    /// Index of the `BasisSpec` in the `Vec<BasisSpec>` of the `Elem`
    pub elem_idx: usize,
}

impl BSAddress {
    /// Create a new BSAddress from an Elem's ID and an index into its `Vec<BasisSpec>`
    pub fn new(elem_id: usize, elem_idx: usize) -> Self {
        Self { elem_id, elem_idx }
    }
}

impl fmt::Display for BSAddress {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}({})", self.elem_id, self.elem_idx)
    }
}
