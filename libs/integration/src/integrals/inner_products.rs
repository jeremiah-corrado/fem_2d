
use super::{ParaDir, BasisFn, ShapeFn, real_gauss_quad, real_gauss_quad_inner, Integral, IntegralResult};
use domain::V2D;


/// The L2 Inner product of two Basis Functions
pub struct L2InnerProduct {
    u_weights: Vec<f64>,
    v_weights: Vec<f64>,
}

impl Integral for L2InnerProduct {
    fn with_weights(u_weights: &[f64], v_weights: &[f64]) -> Self {
        Self {
            u_weights: u_weights.to_vec(), 
            v_weights: v_weights.to_vec(),
        }
    }

    fn integrate<SF: ShapeFn>(
        &self,
        p_dir: ParaDir, 
        q_dir: ParaDir, 
        p_orders: [usize; 2], 
        q_orders: [usize; 2], 
        p_basis: &BasisFn<SF>,
        q_basis: &BasisFn<SF>,
    ) -> IntegralResult {
        IntegralResult::Full(p_basis.glq_scale()
            * q_basis.glq_scale()
            * match (p_dir, q_dir) {
                (ParaDir::U, ParaDir::U) => real_gauss_quad(&self.u_weights, &self.v_weights, |m, n| {
                    V2D::dot(p_basis.f_u(p_orders, [m, n]), q_basis.f_u(q_orders, [m, n]))
                        * p_basis.sample_scale([m, n])
                        * q_basis.sample_scale([m, n])
                }),
                (ParaDir::U, ParaDir::V) => real_gauss_quad(&self.u_weights, &self.v_weights, |m, n| {
                    V2D::dot(p_basis.f_u(p_orders, [m, n]), q_basis.f_v(q_orders, [m, n]))
                        * p_basis.sample_scale([m, n])
                        * q_basis.sample_scale([m, n])
                }),
                (ParaDir::V, ParaDir::U) => real_gauss_quad(&self.u_weights, &self.v_weights, |m, n| {
                    V2D::dot(p_basis.f_v(p_orders, [m, n]), q_basis.f_u(q_orders, [m, n]))
                        * p_basis.sample_scale([m, n])
                        * q_basis.sample_scale([m, n])
                }),
                (ParaDir::V, ParaDir::V) => real_gauss_quad(&self.u_weights, &self.v_weights, |m, n| {
                    V2D::dot(p_basis.f_v(p_orders, [m, n]), q_basis.f_v(q_orders, [m, n]))
                        * p_basis.sample_scale([m, n])
                        * q_basis.sample_scale([m, n])
                }),
            })

    }

    fn integrate_by_parts<SF: ShapeFn>(
        &self,
        p_dir: ParaDir, 
        q_dir: ParaDir, 
        p_orders: [usize; 2], 
        q_orders: [usize; 2], 
        p_basis: &BasisFn<SF>,
        q_basis: &BasisFn<SF>,
    ) -> IntegralResult {
        IntegralResult::Full(p_basis.glq_scale()
            * q_basis.glq_scale()
            * match (p_dir, q_dir) {
                (ParaDir::U, ParaDir::U) => real_gauss_quad_inner(&self.u_weights, &self.v_weights, |m, n| {
                    V2D::dot(p_basis.f_u(p_orders, [m, n]), q_basis.f_u(q_orders, [m, n]))
                        * p_basis.sample_scale([m, n])
                        * q_basis.sample_scale([m, n])
                }),
                (ParaDir::U, ParaDir::V) => real_gauss_quad_inner(&self.u_weights, &self.v_weights, |m, n| {
                    V2D::dot(p_basis.f_u(p_orders, [m, n]), q_basis.f_v(q_orders, [m, n]))
                        * p_basis.sample_scale([m, n])
                        * q_basis.sample_scale([m, n])
                }),
                (ParaDir::V, ParaDir::U) => real_gauss_quad_inner(&self.u_weights, &self.v_weights, |m, n| {
                    V2D::dot(p_basis.f_v(p_orders, [m, n]), q_basis.f_u(q_orders, [m, n]))
                        * p_basis.sample_scale([m, n])
                        * q_basis.sample_scale([m, n])
                }),
                (ParaDir::V, ParaDir::V) => real_gauss_quad_inner(&self.u_weights, &self.v_weights, |m, n| {
                    V2D::dot(p_basis.f_v(p_orders, [m, n]), q_basis.f_v(q_orders, [m, n]))
                        * p_basis.sample_scale([m, n])
                        * q_basis.sample_scale([m, n])
                }),
            })
    }
}
