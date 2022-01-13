mod glq;
mod kol;
mod max_ortho;

use fem_domain::{Elem, Point, M2D, V2D};
use glq::{gauss_quadrature_points, scale_gauss_quad_points};
use std::marker::PhantomData;

pub use kol::KOLShapeFn;
pub use max_ortho::MaxOrthoShapeFn;

/// Hierarchical Shape Function along a single direction (defined over (-1.0, +1.0)).
/// [KOLShapeFn] and [MaxOrthoShapeFn] implement this trait.
/// Alternate Hierarchical Basis Functions can be used by implementing this trait.
pub trait ShapeFn {
    fn with(max_order: usize, points: &[f64], compute_d2: bool) -> Self;

    fn power(&self, n: usize, p: usize) -> f64;
    fn power_d1(&self, n: usize, p: usize) -> f64;
    fn power_d2(&self, n: usize, p: usize) -> f64;

    fn poly(&self, n: usize, p: usize) -> f64;
    fn poly_d1(&self, n: usize, p: usize) -> f64;
    fn poly_d2(&self, n: usize, p: usize) -> f64;
}

/// Structure used to generate [BasisFn]'s over [Elem]'s in a Domain.
/// This structure contains settings for a
pub struct BasisFnSampler<SF: ShapeFn> {
    /// Type of [ShapeFn] used in [BasisFn]'s
    shape_type: PhantomData<SF>,
    /// Maximum u-directed expansion order
    pub i_max: usize,
    /// Maximum v-directed expansion order                
    pub j_max: usize,
    /// Whether the 2nd derivatives of the [ShapeFn]'s will be computed                 
    pub compute_d2: bool,
    /// Gauss Legendre Quadrature points evaluated along u-direction (Defined from -1 to +1)         
    u_points: Vec<f64>,
    /// Gauss Legendre Quadrature points evaluated along v-direction (Defined from -1 to +1)
    v_points: Vec<f64>,
}

impl<SF: ShapeFn> BasisFnSampler<SF> {
    pub fn with(
        num_u_points: usize,
        num_v_points: usize,
        i_max: usize,
        j_max: usize,
        compute_2nd_derivs: bool,
    ) -> (Self, [Vec<f64>; 2]) {
        let (u_points, u_weights) = gauss_quadrature_points(num_u_points, compute_2nd_derivs);
        let (v_points, v_weights) = gauss_quadrature_points(num_v_points, compute_2nd_derivs);

        (
            Self {
                shape_type: PhantomData,
                i_max,
                j_max,
                compute_d2: compute_2nd_derivs,
                u_points,
                v_points,
            },
            [u_weights, v_weights],
        )
    }

    /// Generate a [BasisFn] defined over an [Elem]. Can be defined over a subset of the Element.
    pub fn sample_basis_fn(&self, elem: &Elem, sampled_space: Option<[Point; 2]>) -> BasisFn<SF> {
        BasisFn::with(
            self.i_max,
            self.j_max,
            self.compute_d2,
            &self.u_points,
            &self.v_points,
            elem,
            sampled_space,
        )
    }
}

/// Structure used to evaluate [ShapeFn]'s and their derivatives over some area
pub struct BasisFn<SF: ShapeFn> {
    /// Raw transformation matrices at each sample point. Describes transformation from real space to sampled parametric space.
    pub t: Vec<Vec<M2D>>,
    // Inverse of transformation matrices at each sample point.
    pub ti: Vec<Vec<M2D>>,
    /// Determinants of the "Sampling Jacobian" at each point.
    pub dt: Vec<Vec<f64>>,
    /// Parametric scaling factors (used to scale derivatives in parametric space as necessary)
    pub para_scale: V2D,
    u_shapes: SF,
    v_shapes: SF,
}

impl<SF: ShapeFn> BasisFn<SF> {
    pub fn with(
        i_max: usize,
        j_max: usize,
        compute_d2: bool,
        raw_u_points: &[f64],
        raw_v_points: &[f64],
        elem: &Elem,
        sampled_space: Option<[Point; 2]>,
    ) -> Self {
        let [(u_glq_scale, u_points_scaled), (v_glq_scale, v_points_scaled)] = match sampled_space {
            Some([sample_min, sample_max]) => {
                let para_min = elem.parametric_projection(sample_min);
                let para_max = elem.parametric_projection(sample_max);

                [
                    scale_gauss_quad_points(raw_u_points, para_min[0], para_max[0]),
                    scale_gauss_quad_points(raw_u_points, para_min[1], para_max[1]),
                ]
            }
            None => [(1.0, raw_u_points.to_vec()), (1.0, raw_v_points.to_vec())],
        };

        let t: Vec<Vec<M2D>> = u_points_scaled
            .iter()
            .map(|u| {
                v_points_scaled
                    .iter()
                    .map(|v| elem.parametric_gradient(V2D::from([*u, *v])))
                    .collect()
            })
            .collect();

        let ti: Vec<Vec<M2D>> = t
            .iter()
            .map(|row| row.iter().map(|v| v.inverse()).collect())
            .collect();

        let dt = if sampled_space.is_some() {
            t.iter()
                .map(|row| row.iter().map(|v| v.det()).collect())
                .collect()
        } else {
            vec![vec![1.0; raw_v_points.len()]; raw_u_points.len()]
        };

        Self {
            t,
            ti,
            dt,
            para_scale: V2D::from([v_glq_scale, u_glq_scale]),
            u_shapes: SF::with(i_max, &u_points_scaled, compute_d2),
            v_shapes: SF::with(j_max, &v_points_scaled, compute_d2),
        }
    }

    pub fn f_u(&self, [i, j]: [usize; 2], [m, n]: [usize; 2]) -> V2D {
        self.ti[m][n].u * self.u_shapes.power(i, m) * self.v_shapes.poly(j, n)
    }

    pub fn f_v(&self, [i, j]: [usize; 2], [m, n]: [usize; 2]) -> V2D {
        self.ti[m][n].v * self.u_shapes.poly(i, m) * self.v_shapes.power(j, n)
    }

    pub fn f_u_d1(&self, [i, j]: [usize; 2], [m, n]: [usize; 2], para_scale: &V2D) -> V2D {
        self.ti[m][n].u
            * V2D::from([
                self.u_shapes.power(i, m) * self.v_shapes.poly_d1(j, n),
                self.u_shapes.power_d1(i, m) * self.v_shapes.poly(j, n),
            ]) * para_scale
    }

    pub fn f_v_d1(&self, [i, j]: [usize; 2], [m, n]: [usize; 2], para_scale: &V2D) -> V2D {
        self.ti[m][n].v
            * V2D::from([
                self.u_shapes.poly(i, m) * self.v_shapes.power_d1(j, n),
                self.u_shapes.poly_d1(i, m) * self.v_shapes.power(j, n),
            ]) * para_scale
    }

    pub fn f_u_d2(&self, [i, j]: [usize; 2], [m, n]: [usize; 2], para_scale: &V2D) -> V2D {
        self.ti[m][n].u
            * V2D::from([
                self.u_shapes.power(i, m) * self.v_shapes.poly_d2(j, n),
                self.u_shapes.power_d2(i, m) * self.v_shapes.poly(j, n),
            ]) * para_scale * para_scale
    }

    pub fn f_v_d2(&self, [i, j]: [usize; 2], [m, n]: [usize; 2], para_scale: &V2D) -> V2D {
        self.ti[m][n].v
            * V2D::from([
                self.u_shapes.poly(i, m) * self.v_shapes.power_d2(j, n),
                self.u_shapes.poly_d2(i, m) * self.v_shapes.power(j, n),
            ]) * para_scale * para_scale
    }

    pub fn f_u_dd(&self, [i, j]: [usize; 2], [m, n]: [usize; 2], para_scale: &V2D) -> V2D {
        self.ti[m][n].u * self.u_shapes.power_d1(i, m) * self.v_shapes.poly_d1(j, n) * para_scale[0] * self.para_scale[1]
    }

    pub fn f_v_dd(&self, [i, j]: [usize; 2], [m, n]: [usize; 2], para_scale: &V2D) -> V2D {
        self.ti[m][n].v * self.u_shapes.poly_d1(i, m) * self.v_shapes.power_d1(j, n) * para_scale[0] * self.para_scale[1]
    }

    #[inline]
    pub fn glq_scale(&self) -> f64 {
        self.para_scale.dot_with(&V2D::from([1.0, 1.0]))
    }

    #[inline] pub fn edge_glq_scale(&self, edge_idx: usize) -> f64 {
        match edge_idx {
            0 | 1 => self.para_scale[1],
            2 | 3 => self.para_scale[0],
            _ => panic!("edge_idx must not exceed 3; cannot get glq scaling factor!"),
        }
    }

    #[inline]
    pub fn u_glq_scale(&self) -> f64 {
        self.para_scale[1]
    }

    #[inline]
    pub fn v_glq_scale(&self) -> f64 {
        self.para_scale[0]
    }

    #[inline]
    pub fn deriv_scale(&self) -> &V2D {
        &self.para_scale
    }

    #[inline]
    pub fn sample_scale(&self, [m, n]: [usize; 2]) -> f64 {
        self.dt[m][n]
    }
}
