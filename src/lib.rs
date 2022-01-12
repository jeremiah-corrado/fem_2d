extern crate domain;
extern crate sparse_matrix;
extern crate basis;

extern crate integration;
extern crate eigensolver;

pub use sparse_matrix::{SparseMatrix, AIJMatrix};
pub use basis::{BasisFn, ShapeFn, BasisFnSampler, KOLShapeFn, MaxOrthoShapeFn};
pub use domain::{Elem, M2D, V2D, Point};
pub use integration::{Integral, IntegralResult, CurlProduct, L2InnerProduct};

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
