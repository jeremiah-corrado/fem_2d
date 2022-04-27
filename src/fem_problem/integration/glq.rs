use nalgebra::{DMatrix, SymmetricEigen};

/// 2D Gauss Legendre Quadrature integral of some function F defined over an m by n rectangular region.
/// ```
/// use fem_2d::fem_problem::integration::glq::*;
///
/// // define glq points over `(-1, 1)` in both directions
/// let (u_points, u_weights) = gauss_quadrature_points(10, false);
/// let (v_points, v_weights) = gauss_quadrature_points(10, false);
///
/// // compute the integral of (u^2 * v^2)
/// let solution = real_gauss_quad(&u_weights, &v_weights, |m, n| {
///    u_points[m].powi(2) * v_points[n].powi(2)
/// });
///
/// assert!((solution - 4.0 / 9.0).abs() < 1e-12);
///
/// ```
pub fn real_gauss_quad<F>(u_weights: &[f64], v_weights: &[f64], integrand: F) -> f64
where
    F: Fn(usize, usize) -> f64,
{
    let mut solution = 0.0;
    for (m, u_w) in u_weights.iter().enumerate() {
        let mut inner_solution = 0.0;
        for (n, v_w) in v_weights.iter().enumerate() {
            inner_solution += integrand(m, n) * v_w;
        }
        solution += inner_solution * u_w;
    }
    solution
}

/// 2D Gauss Legendre Quadrature integral of some function F defined over an m by n rectangular region.
///
/// This is the same as [real_gauss_quad] except, the outer edge (the first and last elements of 'u_weights' and 'v_weights') are ignored.
/// Intended to be used in scenarios where BasisFns are defined for By-Parts integration but the solutions on the edges are known to be zero.
///
/// ```
/// use fem_2d::fem_problem::integration::glq::*;
///
/// // define glq points over `[-1, 1]` in both directions
/// let (u_points, u_weights) = gauss_quadrature_points(10, true);
/// let (v_points, v_weights) = gauss_quadrature_points(10, true);
///
/// // compute the integral of (u^2 * v^2)
/// let solution = real_gauss_quad_inner(&u_weights, &v_weights, |m, n| {
///    u_points[m].powi(2) * v_points[n].powi(2)
/// });
///
/// assert!((solution - 4.0 / 9.0).abs() < 1e-12);
/// ```
pub fn real_gauss_quad_inner<F>(u_weights: &[f64], v_weights: &[f64], integrand: F) -> f64
where
    F: Fn(usize, usize) -> f64,
{
    let mut solution = 0.0;
    for (m, u_w) in u_weights
        .iter()
        .enumerate()
        .skip(1)
        .take(u_weights.len() - 2)
    {
        let mut inner_solution = 0.0;
        for (n, v_w) in v_weights
            .iter()
            .enumerate()
            .skip(1)
            .take(v_weights.len() - 2)
        {
            inner_solution += integrand(m, n) * v_w;
        }
        solution += inner_solution * u_w;
    }
    solution
}

/// 1D integral over some function F, which is defined along one edge of a rectangular parametric region.
///
/// ```
/// use fem_2d::fem_problem::integration::glq::*;
///
/// // define glq points over `[-1, 1]` in both directions
/// let (u_points, u_weights) = gauss_quadrature_points(10, true);
/// let (v_points, v_weights) = gauss_quadrature_points(10, true);
///
/// // compute the integral of (u^2 * v^2) along all four "edges"
/// let edge_solutions : Vec<f64> = (0..4).map(|edge_idx| {
///     real_gauss_quad_edge(&u_weights, &v_weights, edge_idx, |m, n| {
///        u_points[m].powi(2) * v_points[n].powi(2)
///    })
/// }).collect();
///
/// assert!((edge_solutions[0] - 2.0 / 3.0).abs() < 1e-12);
/// assert!((edge_solutions[1] - 2.0 / 3.0).abs() < 1e-12);
/// assert!((edge_solutions[2] - 2.0 / 3.0).abs() < 1e-12);
/// assert!((edge_solutions[3] - 2.0 / 3.0).abs() < 1e-12);
/// ```
pub fn real_gauss_quad_edge<F>(
    u_weights: &[f64],
    v_weights: &[f64],
    edge_index: usize,
    integrand: F,
) -> f64
where
    F: Fn(usize, usize) -> f64,
{
    let mut solution = 0.0;
    match edge_index {
        0 => {
            for (m, u_w) in u_weights
                .iter()
                .enumerate()
                .skip(1)
                .take(u_weights.len() - 2)
            {
                solution += integrand(m, 0) * u_w
            }
        }
        1 => {
            for (m, u_w) in u_weights
                .iter()
                .enumerate()
                .skip(1)
                .take(u_weights.len() - 2)
                .rev()
            {
                solution += integrand(m, v_weights.len() - 1) * u_w
            }
        }
        2 => {
            for (n, v_w) in v_weights
                .iter()
                .enumerate()
                .skip(1)
                .take(v_weights.len() - 2)
                .rev()
            {
                solution += integrand(0, n) * v_w
            }
        }
        3 => {
            for (n, v_w) in v_weights
                .iter()
                .enumerate()
                .skip(1)
                .take(v_weights.len() - 2)
            {
                solution += integrand(u_weights.len() - 1, n) * v_w
            }
        }
        _ => unreachable!(),
    }

    solution
}

/// Get a set of n Gauss-Legendre-Quadrature Integration points and weights
///
/// ```
/// use fem_2d::fem_problem::integration::glq::*;
///
/// // generate 10 GLQ points and weights over the range `(-1, 1)`
/// let (points, weights) = gauss_quadrature_points(10, false);
/// assert_eq!(points.len(), 10);
/// assert_eq!(weights.len(), 10);
/// assert!(points.iter().sum::<f64>().abs() < 1e-12);
///
/// // generate 10 GLQ points and weights with points defined on the edges `[-1, 1]`
/// let (points_with_end, weights_with_end) = gauss_quadrature_points(10, true);
/// assert_eq!(points_with_end.len(), 12);
/// assert_eq!(weights_with_end.len(), 12);
/// assert!((-1.0 - points_with_end[0].abs()) < 1e-12);
/// assert!((1.0 - points_with_end[11].abs()) < 1e-12);
///
/// ```
// https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss%E2%80%93Legendre_quadrature
// https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/23972/versions/22/previews/chebfun/examples/quad/html/GaussQuad.html
pub fn gauss_quadrature_points(n: usize, include_endpoints: bool) -> (Vec<f64>, Vec<f64>) {
    let betas: Vec<f64> = (1..n)
        .map(|i| 0.5 / (1.0 - (2.0 * i as f64).powi(-2)).sqrt())
        .collect();

    let polymat: DMatrix<f64> = DMatrix::from_fn(n, n, |r, c| {
        if r == c + 1 {
            betas[r - 1]
        } else if c == r + 1 {
            betas[c - 1]
        } else {
            0.0
        }
    });

    let eigen_decomp = SymmetricEigen::new(polymat);

    let mut xw: Vec<(f64, f64)> = eigen_decomp
        .eigenvalues
        .iter()
        .cloned()
        .zip(
            eigen_decomp
                .eigenvectors
                .row(0)
                .iter()
                .map(|weight| (*weight).powi(2) * 2.0),
        )
        .collect();

    xw.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());

    let (mut points, mut weights): (Vec<_>, Vec<_>) = xw.drain(0..).unzip();

    if include_endpoints {
        points.insert(0, -1.0);
        points.push(1.0);

        weights.insert(0, 1.0);
        weights.push(1.0)
    }

    (points, weights)
}

/// Scale a set of Gauss-Legendre-Quadrature Integration points to fall within a specific range
///
/// ```
/// use fem_2d::fem_problem::integration::glq::*;
/// // generate 10 GLQ points and weights over the range `[-1, 1]`
/// let (points, weights) = gauss_quadrature_points(10, true);
///
/// // scale the points to the range `(-0.75, 0.25)`
/// let (scale, points_scaled) = scale_gauss_quad_points(&points, -0.75, 0.25);
///
/// assert!((-0.75 - points_scaled[0]).abs() < 1e-12);
/// assert!((0.25 - points_scaled[11]).abs() < 1e-12);
/// assert!((0.5 - scale).abs() < 1e-12);
/// ```
pub fn scale_gauss_quad_points(points: &[f64], min: f64, max: f64) -> (f64, Vec<f64>) {
    let scale_factor = (max - min) / 2.0;
    let offset = (max + min) / 2.0;

    (
        scale_factor,
        points
            .iter()
            .map(|x| x * scale_factor + offset)
            .collect::<Vec<f64>>(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    const GLQ_ACCURACY: f64 = 1e-9;
    // test points
    const X_20: [f64; 20] = [
        -0.993128599,
        -0.963971927,
        -0.912234428,
        -0.839116972,
        -0.746331906,
        -0.636053681,
        -0.510867002,
        -0.373706089,
        -0.227785851,
        -0.076526521,
        0.076526521,
        0.227785851,
        0.373706089,
        0.510867002,
        0.636053681,
        0.746331906,
        0.839116972,
        0.912234428,
        0.963971927,
        0.993128599,
    ];
    const W_20: [f64; 20] = [
        0.017614007,
        0.04060143,
        0.062672048,
        0.083276742,
        0.10193012,
        0.118194532,
        0.131688638,
        0.142096109,
        0.149172986,
        0.152753387,
        0.152753387,
        0.149172986,
        0.142096109,
        0.131688638,
        0.118194532,
        0.10193012,
        0.083276742,
        0.062672048,
        0.04060143,
        0.017614007,
    ];

    const X_20_SCALED: [f64; 20] = [
        0.250858925,
        0.254503509,
        0.260970696,
        0.270110379,
        0.281708512,
        0.29549329,
        0.311141625,
        0.328286739,
        0.346526769,
        0.365434185,
        0.384565815,
        0.403473231,
        0.421713261,
        0.438858375,
        0.45450671,
        0.468291488,
        0.479889621,
        0.489029304,
        0.495496491,
        0.499141075,
    ];

    #[test]
    fn glq_point_generation_and_scaling() {
        let (glq_points, glq_weights) = gauss_quadrature_points(20, false);

        for (glq_ref, glq_test) in X_20.iter().zip(glq_points.iter()) {
            assert!((glq_ref - glq_test).abs() < GLQ_ACCURACY);
        }

        for (glq_w_ref, glq_w_test) in W_20.iter().zip(glq_weights.iter()) {
            assert!((glq_w_ref - glq_w_test).abs() < GLQ_ACCURACY);
        }

        let (glq_scale, glq_scaled_points) = scale_gauss_quad_points(&glq_points, 0.25, 0.5);

        assert!((glq_scale - 0.125).abs() < 1e-14);
        for (glq_s_ref, glq_s_test) in X_20_SCALED.iter().zip(glq_scaled_points.iter()) {
            assert!((glq_s_ref - glq_s_test).abs() < GLQ_ACCURACY);
        }
    }

    #[test]
    fn glq_point_generation_with_endpoints() {
        let (glq_points, glq_weights) = gauss_quadrature_points(20, true);

        assert!((glq_points.first().unwrap() + 1.0).abs() < 1e-14);
        assert!((glq_points.last().unwrap() - 1.0).abs() < 1e-14);

        assert!((glq_weights.first().unwrap() - 1.0).abs() < 1e-14);
        assert!((glq_weights.last().unwrap() - 1.0).abs() < 1e-14);

        for (glq_ref, glq_test) in X_20.iter().zip(glq_points.iter().skip(1).take(20)) {
            assert!((glq_ref - glq_test).abs() < GLQ_ACCURACY);
        }

        for (glq_w_ref, glq_w_test) in W_20.iter().zip(glq_weights.iter().skip(1).take(20)) {
            assert!((glq_w_ref - glq_w_test).abs() < GLQ_ACCURACY);
        }

        let (glq_scale, glq_scaled_points) = scale_gauss_quad_points(&glq_points, 0.25, 0.5);

        assert!((glq_scaled_points.first().unwrap() - 0.25).abs() < 1e-14);
        assert!((glq_scaled_points.last().unwrap() - 0.5).abs() < 1e-14);
        assert!((glq_scale - 0.125).abs() < 1e-14);

        for (glq_s_ref, glq_s_test) in X_20_SCALED
            .iter()
            .zip(glq_scaled_points.iter().skip(1).take(20))
        {
            assert!((glq_s_ref - glq_s_test).abs() < GLQ_ACCURACY);
        }
    }
}
