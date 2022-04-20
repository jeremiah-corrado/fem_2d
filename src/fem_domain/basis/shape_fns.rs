pub mod kol {
    use super::super::ShapeFn;

    /// A simple Curl-Conforming Hierarchical Shape Function
    #[derive(Clone, Debug)]
    pub struct KOLShapeFn {
        pows: Vec<Vec<f64>>,
        pows_d1: Vec<Vec<f64>>,
        pows_d2: Vec<Vec<f64>>,
        polys: Vec<Vec<f64>>,
        polys_d1: Vec<Vec<f64>>,
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

        fn tang(&self, n: usize, p: usize) -> f64 {
            self.pows[n][p]
        }

        fn tang_d1(&self, n: usize, p: usize) -> f64 {
            self.pows_d1[n][p]
        }

        fn tang_d2(&self, n: usize, p: usize) -> f64 {
            self.pows_d2[n][p]
        }

        fn norm(&self, n: usize, p: usize) -> f64 {
            self.polys[n][p]
        }

        fn norm_d1(&self, n: usize, p: usize) -> f64 {
            self.polys_d1[n][p]
        }

        fn norm_d2(&self, n: usize, p: usize) -> f64 {
            // coincidentally same as pows_d2
            self.pows_d2[n][p]
        }
    }
}

#[cfg(feature = "max_ortho_basis")]
mod max_ortho {
    use super::super::ShapeFn;

    //https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6470651
    const EUC_NORM_COEFFS: [f64; 12] = [
        0.968246, 2.561738, 0.838525, 4.248161, 0.816397, 5.882766, 0.808509, 1.0, 1.0, 1.0, 1.0,
        1.0,
    ];

    const Q_NUMERATORS: [&[i32]; 12] = [
        &[-1, 0, 1],
        &[0, -3, 0, 3],
        &[-1, 0, -5, 0, 6],
        &[0, -3, 0, -7, 0, 10],
        &[-1, 0, -5, 0, -9, 0, 15],
        &[0, -3, 0, -7, 0, -11, 0, 21],
        &[-1, 0, -5, 0, -9, 0, -13, 0, 28],
        &[0, -3, 0, -7, 0, -11, 0, -15, 0, 36],
        &[-1, 0, -5, 0, -9, 0, -13, 0, -17, 0, 40],
        &[0, -3, 0, -7, 0, -11, 0, -15, 0, -19, 0, 55],
        &[-1, 0, -5, 0, -9, 0, -13, 0, -17, 0, -21, 0, 66],
        &[0, -3, 0, -7, 0, -11, 0, -15, 0, -19, 0, -23, 0, 72],
    ];
    const Q_DENOMINATORS: [i32; 12] = [1, 3, 6, 10, 15, 21, 28, 36, 40, 55, 66, 72];

    const fn get_q_weight_vector<const DIM: usize>() -> [f64; DIM] {
        let mut coeffs = [0.0; DIM];
        let mut i = 0;
        let index = DIM - 3;

        while i < DIM {
            coeffs[i] = (Q_NUMERATORS[index][i] as f64) / (Q_DENOMINATORS[index] as f64);
            i += 1;
        }

        coeffs
    }

    const Q_WEIGHTS: [&[f64]; 11] = [
        &get_q_weight_vector::<3>(),
        &get_q_weight_vector::<4>(),
        &get_q_weight_vector::<5>(),
        &get_q_weight_vector::<6>(),
        &get_q_weight_vector::<7>(),
        &get_q_weight_vector::<8>(),
        &get_q_weight_vector::<9>(),
        &get_q_weight_vector::<10>(),
        &get_q_weight_vector::<11>(),
        &get_q_weight_vector::<12>(),
        &get_q_weight_vector::<13>(),
    ];

    /// An advanced Hierarchical Type Shape Function which maximizes orthogonality between polynomial orders
    ///
    /// Based on: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6470651
    #[derive(Clone, Debug)]
    pub struct MaxOrthoShapeFn {
        pub q_fn: QFunction,
        pub l_fn: LegendrePoly,
    }

    impl ShapeFn for MaxOrthoShapeFn {
        fn with(max_order: usize, points: &[f64], compute_d2: bool) -> Self {
            let l_fn = LegendrePoly::with(max_order as u8, points, compute_d2);
            Self {
                q_fn: QFunction::with(max_order as u8, points, &l_fn, compute_d2),
                l_fn,
            }
        }

        fn tang(&self, n: usize, p: usize) -> f64 {
            self.l_fn.l[n][p]
        }
        fn tang_d1(&self, n: usize, p: usize) -> f64 {
            self.l_fn.d1[n][p]
        }
        fn tang_d2(&self, n: usize, p: usize) -> f64 {
            self.l_fn.d2[n][p]
        }

        fn norm(&self, n: usize, p: usize) -> f64 {
            self.q_fn.q[n][p]
        }
        fn norm_d1(&self, n: usize, p: usize) -> f64 {
            self.q_fn.d1[n][p]
        }
        fn norm_d2(&self, n: usize, p: usize) -> f64 {
            self.q_fn.d2[n][p]
        }
    }

    #[derive(Clone, Debug)]
    pub struct QFunction {
        pub q: Vec<Vec<f64>>,
        pub d1: Vec<Vec<f64>>,
        pub d2: Vec<Vec<f64>>,
    }

    impl QFunction {
        pub fn with(
            max_n: u8,
            points: &[f64],
            leg_poly: &LegendrePoly,
            compute_2nd_derivs: bool,
        ) -> Self {
            if compute_2nd_derivs {
                Self::with_specs_with_2nd_derivs(max_n, points, leg_poly)
            } else {
                Self::with_specs_and_no_2nd_derivs(max_n, points, leg_poly)
            }
        }

        fn with_specs_and_no_2nd_derivs(
            max_n: u8,
            points: &[f64],
            leg_poly: &LegendrePoly,
        ) -> Self {
            let mut values = Vec::with_capacity(max_n as usize);
            let mut primes = Vec::with_capacity(max_n as usize);

            for i in 0..=(max_n as usize) {
                match i {
                    0 => {
                        values.push(points.iter().map(|p| 1.0 - p).collect());
                        primes.push((0..points.len()).map(|_| -1.0).collect());
                    }
                    1 => {
                        values.push(points.iter().map(|p| 1.0 + p).collect());
                        primes.push((0..points.len()).map(|_| 1.0).collect());
                    }
                    _ => {
                        values.push(
                            leg_poly.weighted_value_sum(&Q_WEIGHTS[i - 2], EUC_NORM_COEFFS[i - 2]),
                            // leg_poly.weighted_value_sum(Q_SEGMENT_WEIGHTS[i - 2], 1.0),
                        );
                        primes.push(
                            leg_poly.weighted_prime_sum(&Q_WEIGHTS[i - 2], EUC_NORM_COEFFS[i - 2]),
                            // leg_poly.weighted_prime_sum(Q_SEGMENT_WEIGHTS[i - 2], 1.0),
                        );
                    }
                }
            }

            Self {
                q: values,
                d1: primes,
                d2: Vec::new(),
            }
        }

        fn with_specs_with_2nd_derivs(max_n: u8, points: &[f64], leg_poly: &LegendrePoly) -> Self {
            let n = max_n as usize;
            let np = points.len();

            let mut values = Vec::with_capacity(n);
            let mut primes = Vec::with_capacity(n);
            let mut double_primes = Vec::with_capacity(n);

            for i in 0..=n {
                match i {
                    0 => {
                        values.push(points.iter().map(|p| 1.0 - p).collect());
                        primes.push(vec![-1.0; np]);
                        double_primes.push(vec![0.0; np]);
                    }
                    1 => {
                        values.push(points.iter().map(|p| 1.0 + p).collect());
                        primes.push(vec![1.0; np]);
                        double_primes.push(vec![0.0; np]);
                    }
                    _ => {
                        values.push(
                            leg_poly.weighted_value_sum(&Q_WEIGHTS[i - 2], EUC_NORM_COEFFS[i - 2]),
                        );
                        primes.push(
                            leg_poly.weighted_prime_sum(&Q_WEIGHTS[i - 2], EUC_NORM_COEFFS[i - 2]),
                        );
                        double_primes.push(
                            leg_poly.weighted_double_prime_sum(
                                &Q_WEIGHTS[i - 2],
                                EUC_NORM_COEFFS[i - 2],
                            ),
                        )
                    }
                }
            }

            Self {
                q: values,
                d1: primes,
                d2: double_primes,
            }
        }
    }

    #[derive(Clone, Debug)]
    pub struct LegendrePoly {
        pub l: Vec<Vec<f64>>,
        pub d1: Vec<Vec<f64>>,
        pub d2: Vec<Vec<f64>>,
        has_2nd_derivs: bool,
    }

    impl LegendrePoly {
        pub fn new() -> Self {
            Self {
                l: Vec::new(),
                d1: Vec::new(),
                d2: Vec::new(),
                has_2nd_derivs: false,
            }
        }

        pub fn with(max_n: u8, points: &[f64], compute_2nd_derivs: bool) -> Self {
            if compute_2nd_derivs {
                Self::with_specs_with_2nd_derivs(max_n, points)
            } else {
                Self::with_specs_and_no_2nd_derivs(max_n, points)
            }
        }

        fn with_specs_and_no_2nd_derivs(max_n: u8, points: &[f64]) -> Self {
            let mut values = Vec::with_capacity(max_n as usize);
            let mut primes = Vec::with_capacity(max_n as usize);

            for i in 0..=(max_n as usize) {
                values.push(Vec::with_capacity(points.len()));
                primes.push(Vec::with_capacity(points.len()));

                let i_f = i as f64;
                for (p, &point) in points.iter().enumerate() {
                    match i {
                        0 => {
                            values[i].push(1.0);
                            primes[i].push(0.0);
                        }
                        1 => {
                            values[i].push(point);
                            primes[i].push(1.0);
                        }
                        _ => {
                            let v = ((2.0 * i_f - 1.0) * point * values[i - 1][p]
                                - (i_f - 1.0) * values[i - 2][p])
                                / i_f;
                            values[i].push(v);

                            let p = i_f * values[i - 1][p] + point * primes[i - 1][p];
                            primes[i].push(p)
                        }
                    }
                }
            }

            Self {
                l: values,
                d1: primes,
                d2: Vec::new(),
                has_2nd_derivs: false,
            }
        }

        fn with_specs_with_2nd_derivs(max_n: u8, points: &[f64]) -> Self {
            let n = max_n as usize;
            let num_p = points.len();

            let mut values = Vec::with_capacity(n);
            let mut primes = Vec::with_capacity(n);
            let mut double_primes = Vec::with_capacity(n);

            // construct (1 - x^2) without singularities at the endpoints
            let d2_denom: Vec<f64> = match (
                (points[0].abs() - 1.0).abs() < 1e-15,
                (points[num_p - 1].abs() - 1.0).abs() < 1e-15,
            ) {
                (true, true) => {
                    let mut den = vec![1.0];
                    den.extend(
                        points
                            .iter()
                            .skip(1)
                            .take(num_p - 2)
                            .map(|x| 1.0 - x.powi(2)),
                    );
                    den.push(1.0);
                    den
                }
                (true, false) => {
                    let mut den = vec![1.0];
                    den.extend(
                        points
                            .iter()
                            .skip(1)
                            .take(num_p - 1)
                            .map(|x| 1.0 - x.powi(2)),
                    );
                    den
                }
                (false, true) => {
                    let mut den: Vec<_> = points
                        .iter()
                        .take(num_p - 1)
                        .map(|x| 1.0 - x.powi(2))
                        .collect();
                    den.push(1.0);
                    den
                }
                (false, false) => points.iter().map(|x| (1.0 - x.powi(2))).collect(),
            };

            for i in 0..=n {
                let i_f = i as f64;
                match i {
                    0 => {
                        values.push(vec![1.0; num_p]);
                        primes.push(vec![0.0; num_p]);

                        double_primes.push(vec![0.0; num_p]);
                    }
                    1 => {
                        values.push(Vec::from(points));
                        primes.push(vec![1.0; num_p]);

                        double_primes.push(vec![0.0; num_p]);
                    }
                    2 => {
                        values.push(
                            points
                                .iter()
                                .enumerate()
                                .map(|(p, x)| {
                                    ((2.0 * i_f - 1.0) * x * values[i - 1][p]
                                        - (i_f - 1.0) * values[i - 2][p])
                                        / i_f
                                })
                                .collect(),
                        );
                        primes.push(
                            points
                                .iter()
                                .enumerate()
                                .map(|(p, x)| i_f * values[i - 1][p] + x * primes[i - 1][p])
                                .collect(),
                        );

                        double_primes.push(vec![3.0; num_p]);
                    }
                    _ => {
                        values.push(
                            points
                                .iter()
                                .enumerate()
                                .map(|(p, x)| {
                                    ((2.0 * i_f - 1.0) * x * values[i - 1][p]
                                        - (i_f - 1.0) * values[i - 2][p])
                                        / i_f
                                })
                                .collect(),
                        );
                        primes.push(
                            points
                                .iter()
                                .enumerate()
                                .map(|(p, x)| i_f * values[i - 1][p] + x * primes[i - 1][p])
                                .collect(),
                        );

                        double_primes.push(
                            primes[i]
                                .iter()
                                .zip(values[i].iter())
                                .enumerate()
                                .map(|(p, (prime, value))| {
                                    (2.0 * points[p] * prime - (i_f * (i_f + 1.0)) * value)
                                        / d2_denom[p]
                                })
                                .collect(),
                        )
                    }
                }
            }

            Self {
                l: values,
                d1: primes,
                d2: double_primes,
                has_2nd_derivs: true,
            }
        }

        pub fn weighted_value_sum(&self, weights: &[f64], normalization_coeff: f64) -> Vec<f64> {
            let mut sum = vec![0.0; self.l[0].len()];
            for (order, &weight) in weights.iter().enumerate() {
                for (p, &value) in self.l[order].iter().enumerate() {
                    sum[p] += weight * value;
                }
            }

            for s in sum.iter_mut() {
                *s *= normalization_coeff;
            }

            sum
        }

        pub fn weighted_prime_sum(&self, weights: &[f64], normalization_coeff: f64) -> Vec<f64> {
            let mut sum = vec![0.0; self.l[0].len()];
            for (order, &weight) in weights.iter().enumerate() {
                for (p, &value) in self.d1[order].iter().enumerate() {
                    sum[p] += weight * value;
                }
            }

            for s in sum.iter_mut() {
                *s *= normalization_coeff;
            }

            sum
        }

        pub fn weighted_double_prime_sum(
            &self,
            weights: &[f64],
            normalization_coeff: f64,
        ) -> Vec<f64> {
            assert!(
                self.has_2nd_derivs,
                "2nd Derivatives not evaluated on Legendre Polynomial; cannot compute weighted sum!"
            );
            let mut sum = vec![0.0; self.l[0].len()];

            for (order, &weight) in weights.iter().enumerate() {
                for (p, &value) in self.d2[order].iter().enumerate() {
                    sum[p] += weight * value;
                }
            }

            for s in sum.iter_mut() {
                *s *= normalization_coeff;
            }

            sum
        }
    }
}
