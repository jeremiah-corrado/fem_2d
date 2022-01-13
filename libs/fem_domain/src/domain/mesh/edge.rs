
use super::{ParaDir, h_refinement::Bisection};

#[derive(Debug)]
pub struct Edge {
    id: usize,
    nodes: [usize; 2],
    children: Option<[usize; 2]>,
    parent: Option<(usize, Bisection)>,
    dir: ParaDir,
    boundary: bool,
}