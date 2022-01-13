use super::Point;

#[derive(Debug)]
pub struct Node {
    id: usize,
    coords: Point,
    boundary: bool,
}