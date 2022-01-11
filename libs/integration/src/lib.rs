extern crate domain;
extern crate basis;

mod integrals;

pub use integrals::{Integral, IntegralResult, real_gauss_quad, real_gauss_quad_edge, real_gauss_quad_inner, L2InnerProduct, CurlProduct};

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
