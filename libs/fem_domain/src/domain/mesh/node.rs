use super::Point;

/// A point in 2D space. 
/// Pairs of points describe Edges and Groups of 4 points describe Elems.
#[derive(Debug)]
pub struct Node {
    pub id: usize,
    pub coords: Point,
    pub boundary: bool,
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
