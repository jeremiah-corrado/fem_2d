extern crate basis;
extern crate eigensolver;
extern crate fem_domain;
extern crate rayon;

mod fields;
mod integrals;
mod matrix_filling;

pub use integrals::{
    real_gauss_quad, real_gauss_quad_edge, real_gauss_quad_inner, CurlProduct, Integral,
    IntegralResult, L2InnerProduct,
};

pub use matrix_filling::{fill_matrices, fill_matrices_parallel};

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
