


mod mesh;
mod dof;

pub use mesh::{Mesh, Elem, M2D, V2D, Point, ParaDir};
pub use dof::DoF;

pub struct Domain<'e> {
    mesh: Mesh<'e>,
    dofs: Vec<DoF>,
}