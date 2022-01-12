#![feature(const_fn_floating_point_arithmetic)]

extern crate domain;

mod basis_fn;

pub use basis_fn::{BasisFn, BasisFnSampler, KOLShapeFn, MaxOrthoShapeFn, ShapeFn};

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
