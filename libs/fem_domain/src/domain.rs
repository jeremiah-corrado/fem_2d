mod dof;
mod mesh;

pub use dof::DoF;
pub use mesh::*;

pub struct Domain {
    pub mesh: Mesh,
    pub dofs: Vec<DoF>,
}
