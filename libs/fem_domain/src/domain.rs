mod dof;
mod mesh;

pub use dof::DoF;
pub use mesh::*;

pub struct Domain<'e> {
    pub mesh: Mesh<'e>,
    pub dofs: Vec<DoF>,
}
