pub mod glq;
use super::{Integral, IntegralResult};
use crate::basis::{BasisFn, ShapeFn};
use crate::domain::{dof::basis_spec::BasisDir, mesh::element::Materials, mesh::space::V2D};

use glq::*;

/// <∇ × u, ∇ × ρ>
pub mod curl_curl {
    use super::*;

    /// The L2 Inner product of the Curl of two Basis Functions
    pub struct CurlCurl {
        u_weights: Vec<f64>,
        v_weights: Vec<f64>,
    }

    impl Integral for CurlCurl {
        fn with_weights(u_weights: &[f64], v_weights: &[f64]) -> Self {
            Self {
                u_weights: u_weights.to_vec(),
                v_weights: v_weights.to_vec(),
            }
        }

        fn integrate<SF: ShapeFn>(
            &self,
            p_dir: BasisDir,
            q_dir: BasisDir,
            p_orders: [usize; 2],
            q_orders: [usize; 2],
            p_basis: &BasisFn<SF>,
            q_basis: &BasisFn<SF>,
            materials: &Materials,
        ) -> IntegralResult {
            IntegralResult::Full(
                (1.0 / materials.mu_rel.re)
                    // * p_basis.glq_scale()
                    // * q_basis.glq_scale()
                    * match (p_dir, q_dir) {
                        (BasisDir::U, BasisDir::U) => {
                            real_gauss_quad(&self.u_weights, &self.v_weights, |m, n| {
                                let p_curl = p_basis
                                    .f_u_d1(p_orders, [m, n], q_basis.deriv_scale())
                                    .dot_with(&CURL_OP);
                                let q_curl = q_basis
                                    .f_u_d1(q_orders, [m, n], p_basis.deriv_scale())
                                    .dot_with(&CURL_OP);

                                p_curl * q_curl * max_uv_ratios(p_basis, q_basis, [m, n])
                            })
                        }
                        (BasisDir::U, BasisDir::V) => {
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
                        (BasisDir::V, BasisDir::U) => {
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
                        (BasisDir::V, BasisDir::V) => {
                            real_gauss_quad(&self.u_weights, &self.v_weights, |m, n| {
                                let p_curl = p_basis
                                    .f_v_d1(p_orders, [m, n], q_basis.deriv_scale())
                                    .dot_with(&CURL_OP);
                                let q_curl = q_basis
                                    .f_v_d1(q_orders, [m, n], p_basis.deriv_scale())
                                    .dot_with(&CURL_OP);

                                p_curl * q_curl * max_vu_ratios(p_basis, q_basis, [m, n])
                            })
                        }
                        (_, _) => 0.0,
                    },
            )
        }

        fn integrate_by_parts<SF: ShapeFn>(
            &self,
            p_dir: BasisDir,
            q_dir: BasisDir,
            p_orders: [usize; 2],
            q_orders: [usize; 2],
            p_basis: &BasisFn<SF>,
            q_basis: &BasisFn<SF>,
            materials: &Materials,
        ) -> IntegralResult {
            let surface_term = (1.0 / materials.mu_rel.re)
                * p_basis.glq_scale()
                * q_basis.glq_scale()
                * match (p_dir, q_dir) {
                    (BasisDir::U, BasisDir::U) => {
                        real_gauss_quad_inner(&self.u_weights, &self.v_weights, |m, n| {
                            let p_d2 = p_basis.f_u_d2(p_orders, [m, n], q_basis.deriv_scale());
                            let p_dd = p_basis.f_u_dd(p_orders, [m, n], q_basis.deriv_scale());

                            let p = V2D::from([p_dd[1] + p_d2[0], p_d2[1] + p_dd[0]]);
                            let q = q_basis.f_u(q_orders, [m, n]);

                            V2D::dot(p, q) / q_basis.glq_scale().powi(2)
                        }) * -1.0
                    }
                    (BasisDir::U, BasisDir::V) => {
                        real_gauss_quad_inner(&self.u_weights, &self.v_weights, |m, n| {
                            let p_d2 = p_basis.f_u_d2(p_orders, [m, n], q_basis.deriv_scale());
                            let p_dd = p_basis.f_u_dd(p_orders, [m, n], q_basis.deriv_scale());

                            let p = V2D::from([p_dd[1] + p_d2[0], p_d2[1] + p_dd[0]]);
                            let q = q_basis.f_v(q_orders, [m, n]);

                            V2D::dot(p, q) / q_basis.glq_scale().powi(2)
                        })
                    }
                    (BasisDir::V, BasisDir::U) => {
                        real_gauss_quad_inner(&self.u_weights, &self.v_weights, |m, n| {
                            let p_d2 = p_basis.f_v_d2(p_orders, [m, n], q_basis.deriv_scale());
                            let p_dd = p_basis.f_v_dd(p_orders, [m, n], q_basis.deriv_scale());

                            let p = V2D::from([p_dd[1] + p_d2[0], p_d2[1] + p_dd[0]]);
                            let q = q_basis.f_u(q_orders, [m, n]);

                            V2D::dot(p, q) / q_basis.glq_scale().powi(2)
                        })
                    }
                    (BasisDir::V, BasisDir::V) => {
                        real_gauss_quad_inner(&self.u_weights, &self.v_weights, |m, n| {
                            let p_d2 = p_basis.f_v_d2(p_orders, [m, n], q_basis.deriv_scale());
                            let p_dd = p_basis.f_v_dd(p_orders, [m, n], q_basis.deriv_scale());

                            let p = V2D::from([p_dd[1] + p_d2[0], p_d2[1] + p_dd[0]]);
                            let q = q_basis.f_v(q_orders, [m, n]);

                            V2D::dot(p, q) / q_basis.glq_scale().powi(2)
                        }) * -1.0
                    }
                    (_, _) => 0.0,
                };

            let edge_terms = (0..4)
                .map(|edge_idx| {
                    -1.0 * (1.0 / materials.mu_rel.re)
                        * p_basis.edge_glq_scale(edge_idx)
                        * q_basis.edge_glq_scale(edge_idx)
                        * match (p_dir, q_dir, edge_idx) {
                            (BasisDir::U, BasisDir::U, 0 | 1) => real_gauss_quad_edge(
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

                                    p_curl * q / q_basis.glq_scale() * max_uv_ratios(p_basis, q_basis, [m, n])
                                },
                            ),
                            (BasisDir::V, BasisDir::U, 0 | 1) => real_gauss_quad_edge(
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

                                    p_curl * q / q_basis.glq_scale()
                                },
                            ),
                            (BasisDir::U, BasisDir::V, 2 | 3) => {
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

                                        p_curl * q / q_basis.glq_scale()
                                    },
                                ) * -1.0
                            }
                            (BasisDir::V, BasisDir::V, 2 | 3) => {
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

                                        p_curl * q / q_basis.glq_scale() * max_vu_ratios(p_basis, q_basis, [m, n])
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

    #[inline]
    fn max_uv_ratios<SF: ShapeFn>(p_basis: &BasisFn<SF>, q_basis: &BasisFn<SF>, [m, n]: [usize; 2]) -> f64 {
        ((p_basis.dt[m][n] >= q_basis.dt[m][n]) as u8) as f64 * p_basis.uv_ratio([m, n]) + 
            ((p_basis.dt[m][n] < q_basis.dt[m][n]) as u8) as f64 * q_basis.uv_ratio([m, n])
    }

    #[inline]
    fn max_vu_ratios<SF: ShapeFn>(p_basis: &BasisFn<SF>, q_basis: &BasisFn<SF>, [m, n]: [usize; 2]) -> f64 {
        ((p_basis.dt[m][n] >= q_basis.dt[m][n]) as u8) as f64 * p_basis.vu_ratio([m, n]) + 
            ((p_basis.dt[m][n] < q_basis.dt[m][n]) as u8) as f64 * q_basis.vu_ratio([m, n])
    }

}

/// <u, ρ>
pub mod inner {
    use super::*;

    /// The L2 Inner product of two Basis Functions
    pub struct L2Inner {
        u_weights: Vec<f64>,
        v_weights: Vec<f64>,
    }

    impl Integral for L2Inner {
        fn with_weights(u_weights: &[f64], v_weights: &[f64]) -> Self {
            Self {
                u_weights: u_weights.to_vec(),
                v_weights: v_weights.to_vec(),
            }
        }

        fn integrate<SF: ShapeFn>(
            &self,
            p_dir: BasisDir,
            q_dir: BasisDir,
            p_orders: [usize; 2],
            q_orders: [usize; 2],
            p_basis: &BasisFn<SF>,
            q_basis: &BasisFn<SF>,
            materials: &Materials,
        ) -> IntegralResult {
            IntegralResult::Full(
                materials.eps_rel.re
                    * p_basis.glq_scale()
                    * q_basis.glq_scale()
                    * match (p_dir, q_dir) {
                        (BasisDir::U, BasisDir::U) => {
                            real_gauss_quad(&self.u_weights, &self.v_weights, |m, n| {
                                V2D::dot(
                                    p_basis.f_u(p_orders, [m, n]),
                                    q_basis.f_u(q_orders, [m, n]),
                                ) * partial_max(
                                    p_basis.sample_scale([m, n]),
                                    q_basis.sample_scale([m, n]),
                                )
                            })
                        }
                        (BasisDir::U, BasisDir::V) => {
                            real_gauss_quad(&self.u_weights, &self.v_weights, |m, n| {
                                V2D::dot(
                                    p_basis.f_u(p_orders, [m, n]),
                                    q_basis.f_v(q_orders, [m, n]),
                                ) * partial_max(
                                    p_basis.sample_scale([m, n]),
                                    q_basis.sample_scale([m, n]),
                                )
                            })
                        }
                        (BasisDir::V, BasisDir::U) => {
                            real_gauss_quad(&self.u_weights, &self.v_weights, |m, n| {
                                V2D::dot(
                                    p_basis.f_v(p_orders, [m, n]),
                                    q_basis.f_u(q_orders, [m, n]),
                                ) * partial_max(
                                    p_basis.sample_scale([m, n]),
                                    q_basis.sample_scale([m, n]),
                                )
                            })
                        }
                        (BasisDir::V, BasisDir::V) => {
                            real_gauss_quad(&self.u_weights, &self.v_weights, |m, n| {
                                V2D::dot(
                                    p_basis.f_v(p_orders, [m, n]),
                                    q_basis.f_v(q_orders, [m, n]),
                                ) * partial_max(
                                    p_basis.sample_scale([m, n]),
                                    q_basis.sample_scale([m, n]),
                                )
                            })
                        }
                        (_, _) => 0.0,
                    },
            )
        }

        fn integrate_by_parts<SF: ShapeFn>(
            &self,
            p_dir: BasisDir,
            q_dir: BasisDir,
            p_orders: [usize; 2],
            q_orders: [usize; 2],
            p_basis: &BasisFn<SF>,
            q_basis: &BasisFn<SF>,
            materials: &Materials,
        ) -> IntegralResult {
            IntegralResult::Full(
                materials.eps_rel.re
                    * p_basis.glq_scale()
                    * q_basis.glq_scale()
                    * match (p_dir, q_dir) {
                        (BasisDir::U, BasisDir::U) => {
                            real_gauss_quad_inner(&self.u_weights, &self.v_weights, |m, n| {
                                V2D::dot(
                                    p_basis.f_u(p_orders, [m, n]),
                                    q_basis.f_u(q_orders, [m, n]),
                                ) * partial_max(
                                    p_basis.sample_scale([m, n]),
                                    q_basis.sample_scale([m, n]),
                                )
                            })
                        }
                        (BasisDir::U, BasisDir::V) => {
                            real_gauss_quad_inner(&self.u_weights, &self.v_weights, |m, n| {
                                V2D::dot(
                                    p_basis.f_u(p_orders, [m, n]),
                                    q_basis.f_v(q_orders, [m, n]),
                                ) * partial_max(
                                    p_basis.sample_scale([m, n]),
                                    q_basis.sample_scale([m, n]),
                                )
                            })
                        }
                        (BasisDir::V, BasisDir::U) => {
                            real_gauss_quad_inner(&self.u_weights, &self.v_weights, |m, n| {
                                V2D::dot(
                                    p_basis.f_v(p_orders, [m, n]),
                                    q_basis.f_u(q_orders, [m, n]),
                                ) * partial_max(
                                    p_basis.sample_scale([m, n]),
                                    q_basis.sample_scale([m, n]),
                                )
                            })
                        }
                        (BasisDir::V, BasisDir::V) => {
                            real_gauss_quad_inner(&self.u_weights, &self.v_weights, |m, n| {
                                V2D::dot(
                                    p_basis.f_v(p_orders, [m, n]),
                                    q_basis.f_v(q_orders, [m, n]),
                                ) * partial_max(
                                    p_basis.sample_scale([m, n]),
                                    q_basis.sample_scale([m, n]),
                                )
                            })
                        }
                        (_, _) => 0.0,
                    },
            )
        }
    }

    fn partial_max(v1: f64, v2: f64) -> f64 {
        std::cmp::max_by(v1, v2, |a, b| a.partial_cmp(b).unwrap())
    }
}
