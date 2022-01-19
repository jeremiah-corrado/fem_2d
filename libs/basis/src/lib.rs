#![feature(const_fn_floating_point_arithmetic)]

extern crate fem_domain;
mod basis_fn;

pub use basis_fn::{
    BasisFn, BasisFnSampler, KOLShapeFn, MaxOrthoShapeFn, ParBasisFnSampler, ShapeFn,
};

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
