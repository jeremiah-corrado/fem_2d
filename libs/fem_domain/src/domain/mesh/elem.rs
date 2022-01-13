use super::{Point, V2D, M2D, Element, HRef, h_refinement::{HRefLoc, HLevels}};
use smallvec::SmallVec;

#[derive(Debug)]
pub struct Elem<'e> {
    pub id: usize,
    pub nodes: [usize; 4],
    pub edges: [usize; 4],
    pub element: &'e Element,
    children: Option<(SmallVec<[usize; 4]>, HRef)>,
    parent: Option<HRefLoc>,
    h_levels: HLevels,
    poly_orders: PolyOrders,
}

impl<'e> Elem<'e> {
    pub fn parametric_projection(&self, real: Point) -> V2D {
        self.element.parametric_projection(real)
    } 

    pub fn parametric_gradient(&self, parametric_coords: V2D) -> M2D {
        self.element.parametric_gradient(parametric_coords)
    } 
}

