extern crate basis;
extern crate domain;

mod integrals;

pub use integrals::{
    real_gauss_quad, real_gauss_quad_edge, real_gauss_quad_inner, CurlProduct, Integral,
    IntegralResult, L2InnerProduct,
};

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
