use domain::ParaDir;
use basis::{BasisFn, ShapeFn};

mod inner_products;
mod curl_products;
mod glq_integration;

pub use glq_integration::{real_gauss_quad, real_gauss_quad_edge, real_gauss_quad_inner};
pub use inner_products::L2InnerProduct;
pub use curl_products::CurlProduct;

/// Return type of an [Integral]
pub enum IntegralResult {
    /// Overall Integral Result
    Full(f64),
    /// By-Parts Integral Result (face, [edge 0, edge 1, edge 2, edge 3])
    ByParts(f64, [f64; 4]),
}

impl IntegralResult {
    pub fn surface(&self) -> f64 {
        match self {
            Self::Full(surface_result) => *surface_result,
            Self::ByParts(surface_result, _) => *surface_result,
        }
    }
}

/// A trait to describe an "integrator" which can compute 2D integrals over some function of two [BasisFn]'s
pub trait Integral {
    fn with_weights(u_weights: &[f64], v_weights: &[f64]) -> Self;

    /// Compute an integral between [BasisFn]'s P and Q, where P and Q both have a parametric direction ([ParaDir]) and orders i and j.
    fn integrate<SF: ShapeFn>(
        &self,
        p_dir: ParaDir, 
        q_dir: ParaDir, 
        p_orders: [usize; 2], 
        q_orders: [usize; 2], 
        p_basis: &BasisFn<SF>,
        q_basis: &BasisFn<SF>,
    ) -> IntegralResult;

    /// Compute an integral-by-parts between [BasisFn]'s P and Q, where P and Q both have a parametric direction ([ParaDir]) and orders i and j.
    /// This function may still return a "Full" integral result if the solution is known to be zero along the edges.
    fn integrate_by_parts<SF: ShapeFn>(
        &self,
        p_dir: ParaDir, 
        q_dir: ParaDir, 
        p_orders: [usize; 2], 
        q_orders: [usize; 2], 
        p_basis: &BasisFn<SF>,
        q_basis: &BasisFn<SF>,
    ) -> IntegralResult;
}