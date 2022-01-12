extern crate basis;
extern crate domain;
extern crate sparse_matrix;

extern crate eigensolver;
extern crate integration;

pub use basis::{BasisFn, BasisFnSampler, KOLShapeFn, MaxOrthoShapeFn, ShapeFn};
pub use domain::{Elem, Point, M2D, V2D};
pub use integration::{CurlProduct, Integral, IntegralResult, L2InnerProduct};
pub use sparse_matrix::{AIJMatrix, SparseMatrix};

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
