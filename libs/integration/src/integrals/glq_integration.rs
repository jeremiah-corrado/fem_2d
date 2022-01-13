/// 2D Gauss Legendre Quadrature integral of some function F defined over an m by n rectangular region.
/// It is assumed that u_weights.len() == m and v_weights.len() == n.
pub fn real_gauss_quad<F>(u_weights: &Vec<f64>, v_weights: &Vec<f64>, integrand: F) -> f64
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
/// This is the same as [real_gauss_quad] except, the outer edge (the first and last elements of 'u_weights' and 'v_weights') is ignored.
pub fn real_gauss_quad_inner<F>(u_weights: &Vec<f64>, v_weights: &Vec<f64>, integrand: F) -> f64
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
pub fn real_gauss_quad_edge<F>(
    u_weights: &Vec<f64>,
    v_weights: &Vec<f64>,
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
