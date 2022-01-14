mod dof;
mod mesh;

pub use dof::DoF;
pub use mesh::*;

pub struct Domain {
    pub mesh: Mesh,
    pub dofs: Vec<DoF>,
}

impl Domain {
    pub fn blank() -> Self {
        Self {
            mesh: Mesh::blank(),
            dofs: Vec::new(),
        }
    }

    pub fn elems<'a>(&'a self) -> impl Iterator<Item = &'a Elem> + '_ {
        self.mesh.elems.iter()
    }

    pub fn edges<'a>(&'a self) -> impl Iterator<Item = &'a Edge> + '_ {
        self.mesh.edges.iter()
    }

    pub fn nodes<'a>(&'a self) -> impl Iterator<Item = &'a Node> + '_ {
        self.mesh.nodes.iter()
    }
}