
mod space;
mod elem;
mod edge;
mod node;
mod element;
mod h_refinement;

pub use space::{Point, M2D, V2D, ParaDir};
pub use elem::Elem;
pub use node::Node;
pub use edge::Edge;
pub use element::Element;
pub use h_refinement::HRef;

pub struct Mesh<'e> {
    elements: Vec<Element>,
    elems: Vec<Elem<'e>>,
    nodes: Vec<Node>,
    edges: Vec<Edge>,
}