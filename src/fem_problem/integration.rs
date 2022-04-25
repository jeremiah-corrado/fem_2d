use crate::fem_domain::basis::{HierCurlBasisFn, HierCurlBasisFnSpace};
use crate::fem_domain::domain::{dof::basis_spec::BasisDir, mesh::element::Materials};

/// Methods to assist in Gauss-Legendre-Quadrature integration
pub mod glq;

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
    /// Retrieve the full solution regardless of the variant
    /// * Full: yields the solution as is
    /// * ByPars: yields "face + edges.sum()"
    pub fn full_solution(self) -> f64 {
        match self {
            Self::Full(full) => full,
            Self::ByParts(face, edges) => face + edges.iter().sum::<f64>(),
        }
    }

    /// Retrieve the `face` and `edge` solutions separately, panicking if the variant is Full
    pub fn unwrap_parts(self) -> (f64, [f64; 4]) {
        match self {
            Self::Full(_) => {
                panic!("Integral solution was computed in one part; cannot get By-Parts solution!")
            }
            Self::ByParts(face, edges) => (face, edges),
        }
    }

    /// Retrieve the solution over the `face` of the integrated area regardless of the variant
    pub fn get_face(&self) -> f64 {
        match self {
            Self::Full(full) => *full,
            Self::ByParts(face, _) => *face,
        }
    }

    /// Retrieve the solution over the `edges` of the integrated area regardless of the variant. (returns an array of zeros for the `Full` variant)
    pub fn get_edges(&self) -> [f64; 4] {
        match self {
            Self::Full(_) => [0.0; 4],
            Self::ByParts(_, edges) => *edges,
        }
    }
}

// TODO: make the use of &Materials generic s.t. multiple problems can leverage identical Integrals with slight variations in material parameter invocation

/// A trait to describe an "integrator" which can compute 2D integrals over some function of two Hierarchical Curl-Conforming Basis Functions
pub trait HierCurlIntegral: Sync + Send {
    /// Assign a set of Gauss-Legendre-Quadrature weights to this integrator.
    ///
    /// The weight vectors must match the dimension of the [BasisFn]s used in later calls to `integrate` or `integrate_by_parts`
    fn with_weights(u_weights: &[f64], v_weights: &[f64]) -> Self;

    /// Compute an integral between [BasisFn]'s P and Q, where P and Q both have a direction ([BasisDir]) and orders `i` and `j`.
    fn integrate<BSpace: HierCurlBasisFnSpace>(
        &self,
        p_dir: BasisDir,
        q_dir: BasisDir,
        p_orders: [usize; 2],
        q_orders: [usize; 2],
        p_basis: &HierCurlBasisFn<BSpace>,
        q_basis: &HierCurlBasisFn<BSpace>,
        materials: &Materials,
    ) -> IntegralResult;

    /// Compute an integral-by-parts between [BasisFn]'s P and Q, where P and Q both have a direction ([BasisDir]) and orders `i` and `j`.
    ///
    /// This function may still return a the `Full` variant of [IntegralResult] if the solution is known to be zero along the edges.
    fn integrate_by_parts<BSpace: HierCurlBasisFnSpace>(
        &self,
        p_dir: BasisDir,
        q_dir: BasisDir,
        p_orders: [usize; 2],
        q_orders: [usize; 2],
        p_basis: &HierCurlBasisFn<BSpace>,
        q_basis: &HierCurlBasisFn<BSpace>,
        materials: &Materials,
    ) -> IntegralResult;
}
