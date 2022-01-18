mod aij_matrix;

pub use aij_matrix::AIJMatrix;

/// Generalized Eigenproblem 
/// 
/// Au = λBu
pub struct GEP {
    a: AIJMatrix,
    b: AIJMatrix,
}
