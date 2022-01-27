use crate::basis::{BasisFn, ShapeFn};
use crate::domain::dof::basis_spec::BasisDir;

/// Specific Implementations of the `Integral` Trait
pub mod integrals;

/// Return type of an [Integral]
pub enum IntegralResult {
    /// Overall Integral Result
    Full(f64),
    /// By-Parts Integral Result (face, [edge 0, edge 1, edge 2, edge 3])
    ByParts(f64, [f64; 4]),
}

impl IntegralResult {
    /// get the full solution regardless of the variant
    /// * Full: yields the solution as is
    /// * ByPars: yields "face + edges.sum()""
    pub fn full_solution(self) -> f64 {
        match self {
            Self::Full(full) => full,
            Self::ByParts(face, edges) => face + edges.iter().sum::<f64>(),
        }
    }

    /// get the `face` and `edge` solutions, panicking if the variant is Full
    pub fn unwrap_parts(self) -> (f64, [f64; 4]) {
        match self {
            Self::Full(_) => {
                panic!("Integral solution was computed in one part; cannot get By-Parts solution!")
            }
            Self::ByParts(face, edges) => (face, edges),
        }
    }

    /// get the solution over the `face` of the integrated area regardless of the variant
    pub fn get_face(&self) -> f64 {
        match self {
            Self::Full(full) => *full,
            Self::ByParts(face, _) => *face,
        }
    }

    /// get the solution over the `edges` of the integrated area regardless of the variant
    ///
    /// returns an array of zeros for the `Full` variant
    pub fn get_edges(&self) -> [f64; 4] {
        match self {
            Self::Full(_) => [0.0; 4],
            Self::ByParts(_, edges) => *edges,
        }
    }
}

/// A trait to describe an "integrator" which can compute 2D integrals over some function of two [BasisFn]'s
pub trait Integral: Sync + Send {
    /// Construct the Integral with u and v directed Gauss-Leg-Quad weights.
    ///
    /// The weight vectors must match the dimension of the [BasisFn]s used in later calls to `integrate` or `integrate_by_parts`
    fn with_weights(u_weights: &[f64], v_weights: &[f64]) -> Self;

    /// Compute an integral between [BasisFn]'s P and Q, where P and Q both have a direction ([BasisDir]) and orders `i` and `j`.
    fn integrate<SF: ShapeFn>(
        &self,
        p_dir: BasisDir,
        q_dir: BasisDir,
        p_orders: [usize; 2],
        q_orders: [usize; 2],
        p_basis: &BasisFn<SF>,
        q_basis: &BasisFn<SF>,
    ) -> IntegralResult;

    /// Compute an integral-by-parts between [BasisFn]'s P and Q, where P and Q both have a direction ([BasisDir]) and orders `i` and `j`.
    ///
    /// This function may still return a the `Full` variant of [IntegralResult] if the solution is known to be zero along the edges.
    fn integrate_by_parts<SF: ShapeFn>(
        &self,
        p_dir: BasisDir,
        q_dir: BasisDir,
        p_orders: [usize; 2],
        q_orders: [usize; 2],
        p_basis: &BasisFn<SF>,
        q_basis: &BasisFn<SF>,
    ) -> IntegralResult;
}