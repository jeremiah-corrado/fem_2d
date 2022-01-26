#![cfg_attr(feature="max_ortho_basis", feature(const_fn_floating_point_arithmetic))]

extern crate fem_domain;
mod basis_fn;

pub use basis_fn::{
    BasisFn, BasisFnSampler, KOLShapeFn, ParBasisFnSampler, ShapeFn,
};

#[cfg(feature="max_ortho_basis")]
pub use basis_fn::{
    MaxOrthoShapeFn,
};

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
