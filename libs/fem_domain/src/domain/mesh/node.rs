use super::Point;
use smallvec::SmallVec;

#[derive(Debug)]
pub struct Node {
    id: usize,
    coords: Point,
    boundary: bool,
}

impl Node {
    pub fn new(id: usize, coords: Point, boundary: bool) -> Self {
        Self {
            id,
            coords,
            boundary,
        }
    }
}