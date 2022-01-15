use super::h_refinement::HRefError;
use super::Point;
use json::JsonValue;

/// A point in 2D space.
/// Pairs of points describe Edges and Groups of 4 points describe Elems.
#[derive(Debug)]
pub struct Node {
    pub id: usize,
    pub coords: Point,
    pub boundary: bool,
}

impl Node {
    /// construct a new Node defined by its location in real space
    pub fn new(id: usize, coords: Point, boundary: bool) -> Self {
        Self {
            id,
            coords,
            boundary,
        }
    }

    /// Produce a Json Object that describes this Node
    pub fn to_json(&self) -> JsonValue {
        object! {
            "id": self.id,
            "boundary": self.boundary,
            "point": self.coords,
        }
    }
}
