extern crate domain;
extern crate sparse_matrix;
extern crate basis;

extern crate integration;
extern crate eigensolver;

pub use basis::{BasisFn, ShapeFn, BasisFnSampler};
pub use domain::{Elem, M2D, V2D, Point};

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
