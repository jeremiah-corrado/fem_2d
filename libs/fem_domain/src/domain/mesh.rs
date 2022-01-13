mod edge;
mod elem;
mod element;
mod h_refinement;
mod node;
mod p_refinement;
mod space;

pub use edge::Edge;
pub use elem::{Elem, ElemUninit};
pub use element::{Element, Materials};
pub use h_refinement::{HRef, HRefError, HRefLoc, Quadrant, Bisection};
pub use p_refinement::{PRef, PRefError};
pub use node::Node;
pub use space::{ParaDir, Point, M2D, V2D};

/// Minimum Edge length in parametric space. h-Refinements will fail after edges are smaller than this value.
pub const MIN_EDGE_LENGTH: f64 = 1e-6;

/// Maximum Polynomial expansion. p-Refinements will fail when Elem's expansion orders exceed this value.
pub const MAX_POLYNOMIAL_ORDER: u8 = 20;

/// Information used to Define the geometric structure and refinement state of a Domain.
pub struct Mesh<'e> {
    pub elements: Vec<Element>,
    pub elems: Vec<Elem<'e>>,
    pub nodes: Vec<Node>,
    pub edges: Vec<Edge>,
}
