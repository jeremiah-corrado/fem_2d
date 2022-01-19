extern crate basis;
extern crate fem_domain;

extern crate eigensolver;
extern crate integration;

pub use basis::{BasisFn, BasisFnSampler, KOLShapeFn, MaxOrthoShapeFn, ShapeFn};
pub use fem_domain::{Domain, Mesh, Point, M2D, V2D, DoF, HRef, PRef, HRefError, PRefError};
pub use integration::{CurlProduct, Integral, IntegralResult, L2InnerProduct, fill_matrices, fill_matrices_parallel};
pub use eigensolver::{GEP, SparseMatrix};

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
