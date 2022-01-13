use super::{
    real_gauss_quad, real_gauss_quad_edge, real_gauss_quad_inner, BasisFn, Integral,
    IntegralResult, ParaDir, ShapeFn,
};
use fem_domain::V2D;

/// The L2 Inner product of the Curl of two Basis Functions
pub struct CurlProduct {
    u_weights: Vec<f64>,
    v_weights: Vec<f64>,
}

impl Integral for CurlProduct {
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
        IntegralResult::Full(
            p_basis.glq_scale()
                * q_basis.glq_scale()
                * match (p_dir, q_dir) {
                    (ParaDir::U, ParaDir::U) => {
                        real_gauss_quad(&self.u_weights, &self.v_weights, |m, n| {
                            let p_curl = p_basis
                                .f_u_d1(p_orders, [m, n], q_basis.deriv_scale())
                                .dot_with(&CURL_OP);
                            let q_curl = q_basis
                                .f_u_d1(q_orders, [m, n], p_basis.deriv_scale())
                                .dot_with(&CURL_OP);

                            p_curl * q_curl
                        })
                    }
                    (ParaDir::U, ParaDir::V) => {
                        real_gauss_quad(&self.u_weights, &self.v_weights, |m, n| {
                            let p_curl = p_basis
                                .f_u_d1(p_orders, [m, n], q_basis.deriv_scale())
                                .dot_with(&CURL_OP);
                            let q_curl = q_basis
                                .f_v_d1(q_orders, [m, n], p_basis.deriv_scale())
                                .dot_with(&CURL_OP);

                            p_curl * q_curl
                        })
                    }
                    (ParaDir::V, ParaDir::U) => {
                        real_gauss_quad(&self.u_weights, &self.v_weights, |m, n| {
                            let p_curl = p_basis
                                .f_v_d1(p_orders, [m, n], q_basis.deriv_scale())
                                .dot_with(&CURL_OP);
                            let q_curl = q_basis
                                .f_u_d1(q_orders, [m, n], p_basis.deriv_scale())
                                .dot_with(&CURL_OP);

                            p_curl * q_curl
                        })
                    }
                    (ParaDir::V, ParaDir::V) => {
                        real_gauss_quad(&self.u_weights, &self.v_weights, |m, n| {
                            let p_curl = p_basis
                                .f_v_d1(p_orders, [m, n], q_basis.deriv_scale())
                                .dot_with(&CURL_OP);
                            let q_curl = q_basis
                                .f_v_d1(q_orders, [m, n], p_basis.deriv_scale())
                                .dot_with(&CURL_OP);

                            p_curl * q_curl
                        })
                    }
                },
        )
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
        let surface_term = p_basis.glq_scale()
            * q_basis.glq_scale()
            * match (p_dir, q_dir) {
                (ParaDir::U, ParaDir::U) => {
                    real_gauss_quad_inner(&self.u_weights, &self.v_weights, |m, n| {
                        let p_d2 = p_basis.f_u_d2(p_orders, [m, n], q_basis.deriv_scale());
                        let p_dd = p_basis.f_u_dd(p_orders, [m, n], q_basis.deriv_scale());

                        let p = V2D::from([p_dd[1] + p_d2[0], p_d2[1] + p_dd[0]]);
                        let q = q_basis.f_u(q_orders, [m, n]);

                        V2D::dot(p, q)
                    }) * -1.0
                }
                (ParaDir::U, ParaDir::V) => {
                    real_gauss_quad_inner(&self.u_weights, &self.v_weights, |m, n| {
                        let p_d2 = p_basis.f_u_d2(p_orders, [m, n], q_basis.deriv_scale());
                        let p_dd = p_basis.f_u_dd(p_orders, [m, n], q_basis.deriv_scale());

                        let p = V2D::from([p_dd[1] + p_d2[0], p_d2[1] + p_dd[0]]);
                        let q = q_basis.f_v(q_orders, [m, n]);

                        V2D::dot(p, q)
                    })
                }
                (ParaDir::V, ParaDir::U) => {
                    real_gauss_quad_inner(&self.u_weights, &self.v_weights, |m, n| {
                        let p_d2 = p_basis.f_v_d2(p_orders, [m, n], q_basis.deriv_scale());
                        let p_dd = p_basis.f_v_dd(p_orders, [m, n], q_basis.deriv_scale());

                        let p = V2D::from([p_dd[1] + p_d2[0], p_d2[1] + p_dd[0]]);
                        let q = q_basis.f_u(q_orders, [m, n]);

                        V2D::dot(p, q)
                    })
                }
                (ParaDir::V, ParaDir::V) => {
                    real_gauss_quad_inner(&self.u_weights, &self.v_weights, |m, n| {
                        let p_d2 = p_basis.f_v_d2(p_orders, [m, n], q_basis.deriv_scale());
                        let p_dd = p_basis.f_v_dd(p_orders, [m, n], q_basis.deriv_scale());

                        let p = V2D::from([p_dd[1] + p_d2[0], p_d2[1] + p_dd[0]]);
                        let q = q_basis.f_v(q_orders, [m, n]);

                        V2D::dot(p, q)
                    }) * -1.0
                }
            };

        let edge_terms = (0..4)
            .map(|edge_idx| {
                p_basis.edge_glq_scale(edge_idx)
                    * q_basis.edge_glq_scale(edge_idx)
                    * match (p_dir, q_dir, edge_idx) {
                        (ParaDir::U, ParaDir::U, 0 | 1) => real_gauss_quad_edge(
                            &self.u_weights,
                            &self.v_weights,
                            edge_idx,
                            |m, n| {
                                let p_curl = p_basis
                                    .f_u_d1(p_orders, [m, n], q_basis.deriv_scale())
                                    .dot_with(&CURL_OP);
                                let q = q_basis
                                    .f_u(q_orders, [m, n])
                                    .dot_with(&EDGE_UNIT_VECTORS[edge_idx]);

                                p_curl * q
                            },
                        ),
                        (ParaDir::V, ParaDir::U, 0 | 1) => real_gauss_quad_edge(
                            &self.u_weights,
                            &self.v_weights,
                            edge_idx,
                            |m, n| {
                                let p_curl = p_basis
                                    .f_v_d1(p_orders, [m, n], q_basis.deriv_scale())
                                    .dot_with(&CURL_OP);
                                let q = q_basis
                                    .f_u(q_orders, [m, n])
                                    .dot_with(&EDGE_UNIT_VECTORS[edge_idx]);

                                p_curl * q
                            },
                        ),
                        (ParaDir::U, ParaDir::V, 2 | 3) => {
                            real_gauss_quad_edge(
                                &self.u_weights,
                                &self.v_weights,
                                edge_idx,
                                |m, n| {
                                    let p_curl = p_basis
                                        .f_u_d1(p_orders, [m, n], q_basis.deriv_scale())
                                        .dot_with(&CURL_OP);
                                    let q = q_basis
                                        .f_v(q_orders, [m, n])
                                        .dot_with(&EDGE_UNIT_VECTORS[edge_idx]);

                                    p_curl * q
                                },
                            ) * -1.0
                        }
                        (ParaDir::V, ParaDir::V, 2 | 3) => {
                            real_gauss_quad_edge(
                                &self.u_weights,
                                &self.v_weights,
                                edge_idx,
                                |m, n| {
                                    let p_curl = p_basis
                                        .f_v_d1(p_orders, [m, n], q_basis.deriv_scale())
                                        .dot_with(&CURL_OP);
                                    let q = q_basis
                                        .f_v(q_orders, [m, n])
                                        .dot_with(&EDGE_UNIT_VECTORS[edge_idx]);

                                    p_curl * q
                                },
                            ) * -1.0
                        }
                        (_, _, _) => 0.0,
                    }
            })
            .collect::<Vec<f64>>()
            .try_into()
            .unwrap();

        IntegralResult::ByParts(surface_term, edge_terms)
    }
}

const CURL_OP: V2D = V2D::from([-1.0, 1.0]);

const EDGE_UNIT_VECTORS: [V2D; 4] = [
    V2D::from([-1.0, 0.0]),
    V2D::from([1.0, 0.0]),
    V2D::from([0.0, -1.0]),
    V2D::from([0.0, 1.0]),
];
