extern crate nalgebra;
use nalgebra::{DMatrix, SymmetricEigen};

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
