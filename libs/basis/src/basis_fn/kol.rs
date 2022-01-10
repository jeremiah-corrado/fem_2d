use super::ShapeFn;

pub struct KOLShapeFn {
    pows: Vec<Vec<f64>>,
    pows_d1: Vec<Vec<f64>>,
    pows_d2: Vec<Vec<f64>>,
    polys: Vec<Vec<f64>>,
    polys_d1: Vec<Vec<f64>>,
    has_d2: bool,
    dim: usize,
    order: usize,
}

impl KOLShapeFn {
    fn new_with_d2(n_max: usize, points: &[f64]) -> Self {
        let mut pows = Vec::with_capacity(n_max + 1);
        let mut pows_d1 = Vec::with_capacity(n_max + 1);
        let mut pows_d2 = Vec::with_capacity(n_max + 1);

        let mut polys = Vec::with_capacity(n_max + 1);
        let mut polys_d1 = Vec::with_capacity(n_max + 1);

        let num_points = points.len();

        for n in 0..=n_max {
            let n_ = n as f64;
            match n {
                0 => {
                    polys.push(points.iter().map(|x| 1.0 - x).collect::<Vec<f64>>());
                    polys_d1.push(vec![-1.0; num_points]);

                    pows.push(vec![1.0; num_points]);
                    pows_d1.push(vec![0.0; num_points]);
                    pows_d2.push(vec![0.0; num_points]);
                }
                1 => {
                    polys.push(points.iter().map(|x| 1.0 + x).collect::<Vec<f64>>());
                    polys_d1.push(vec![1.0; num_points]);

                    pows.push(points.to_vec().clone());
                    pows_d1.push(vec![1.0; num_points]);
                    pows_d2.push(vec![0.0; num_points]);
                }
                _ => {
                    pows.push(
                        pows[n - 1]
                            .iter()
                            .zip(points.iter())
                            .map(|(pow_prev, x)| pow_prev * x)
                            .collect(),
                    );
                    pows_d1.push(pows[n - 1].iter().map(|pow_prev| n_ * pow_prev).collect());

                    if n == 2 {
                        pows_d2.push(vec![2.0; num_points]);
                    } else if n == 3 {
                        pows_d2.push(points.iter().map(|p| 6.0 * p).collect());
                    } else {
                        pows_d2.push(
                            pows[n - 2]
                                .iter()
                                .map(|pow| n_ * (n_ - 1.0) * pow)
                                .collect(),
                        );
                    }

                    if n % 2 == 0 {
                        polys.push(pows[n].iter().map(|pow| pow - 1.0).collect());
                        polys_d1.push(pows_d1[n].clone());
                    } else {
                        polys.push(
                            pows[n]
                                .iter()
                                .zip(points.iter())
                                .map(|(pow, x)| pow - x)
                                .collect(),
                        );
                        polys_d1.push(pows_d1[n].iter().map(|pow_d1| pow_d1 - 1.0).collect());
                    }
                }
            }
        }

        Self {
            pows,
            pows_d1,
            pows_d2,
            polys,
            polys_d1,
            has_d2: true,
            dim: points.len(),
            order: n_max,
        }
    }

    fn new_without_d2(n_max: usize, points: &[f64]) -> Self {
        let mut pows = Vec::with_capacity(n_max + 1);
        let mut pows_d1 = Vec::with_capacity(n_max + 1);

        let mut polys = Vec::with_capacity(n_max + 1);
        let mut polys_d1 = Vec::with_capacity(n_max + 1);

        let num_points = points.len();

        for n in 0..=n_max {
            let n_ = n as f64;
            match n {
                0 => {
                    polys.push(points.iter().map(|x| 1.0 - x).collect::<Vec<f64>>());
                    polys_d1.push(vec![-1.0; num_points]);

                    pows.push(vec![1.0; num_points]);
                    pows_d1.push(vec![0.0; num_points]);
                }
                1 => {
                    polys.push(points.iter().map(|x| 1.0 + x).collect::<Vec<f64>>());
                    polys_d1.push(vec![1.0; num_points]);

                    pows.push(points.to_vec().clone());
                    pows_d1.push(vec![1.0; num_points]);
                }
                _ => {
                    pows.push(
                        pows[n - 1]
                            .iter()
                            .zip(points.iter())
                            .map(|(pow_prev, x)| pow_prev * x)
                            .collect(),
                    );
                    pows_d1.push(pows[n - 1].iter().map(|pow_prev| n_ * pow_prev).collect());

                    if n % 2 == 0 {
                        polys.push(pows[n].iter().map(|pow| pow - 1.0).collect());
                        polys_d1.push(pows_d1[n].clone());
                    } else {
                        polys.push(
                            pows[n]
                                .iter()
                                .zip(points.iter())
                                .map(|(pow, x)| pow - x)
                                .collect(),
                        );
                        polys_d1.push(pows_d1[n].iter().map(|pow_d1| pow_d1 - 1.0).collect());
                    }
                }
            }
        }

        Self {
            pows,
            pows_d1,
            pows_d2: Vec::new(),
            polys,
            polys_d1,
            has_d2: false,
            dim: points.len(),
            order: n_max,
        }
    }
}

impl ShapeFn for KOLShapeFn {
    fn with(n_max: usize, points: &[f64], compute_2nd_deriv: bool) -> Self {
        if compute_2nd_deriv {
            Self::new_with_d2(n_max, points)
        } else {
            Self::new_without_d2(n_max, points)
        }
    }

    fn power(&self, n: usize, p: usize) -> f64 {
        self.pows[n][p]
    }

    fn power_d1(&self, n: usize, p: usize) -> f64 {
        self.pows_d1[n][p]
    }

    fn power_d2(&self, n: usize, p: usize) -> f64 {
        self.pows_d2[n][p]
    }

    fn poly(&self, n: usize, p: usize) -> f64 {
        self.polys[n][p]
    }

    fn poly_d1(&self, n: usize, p: usize) -> f64 {
        self.polys_d1[n][p]
    }

    fn poly_d2(&self, n: usize, p: usize) -> f64 {
        // coincidentally same as pows_d2
        self.pows_d2[n][p]
    }
}
