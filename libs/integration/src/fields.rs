use basis::{BasisFn, BasisFnSampler, ShapeFn};
use fem_domain::{BasisDir, Domain};
use std::collections::{BTreeMap, HashMap};

pub struct UniformFieldSpace<'d> {
    quantities: HashMap<String, FieldQuantity>,
    points: [Vec<f64>; 2],
    densities: [usize; 2],
    domain: &'d Domain,
}

impl<'d> UniformFieldSpace<'d> {
    /// Generate a FieldSpace over a [Domain]
    pub fn new(domain: &'d Domain, densities: [usize; 2]) -> Self {
        Self {
            quantities: HashMap::new(),
            points: [
                uniform_range(-1.0, 1.0, densities[0]),
                uniform_range(-1.0, 1.0, densities[1]),
            ],
            densities,
            domain,
        }
    }

    /// Use an eigenvector and associated [ShapeFn] to compute the X and Y fields over the `Domain`
    /// 
    /// The X and Y fields will be stored as X_{vector_name} and Y_{vector_name} respectively. The Names are returned in an array in that order.
    pub fn xy_fields<SF: ShapeFn>(
        &mut self,
        vector_name: String,
        eigenvector: Vec<f64>,
    ) -> Result<[String; 2], String> {
        if eigenvector.len() != self.domain.dofs.len() {
            Err(format!(
                "NDofs != Eigenvector length ({} != {}); Cannot compute xy fields over Domain",
                self.domain.dofs.len(),
                eigenvector.len()
            ))
        } else {
            let x_q_name = format!("X_{}", vector_name);
            let y_q_name = format!("Y_{}", vector_name);

            let mut x_quantity = FieldQuantity::new(&x_q_name);
            let mut y_quantity = FieldQuantity::new(&y_q_name);

            let [i_max, j_max] = self.domain.mesh.max_expansion_orders();

            for shell_elem in self.domain.mesh.elems.iter().filter(|e| !e.has_children()) {
                let mut x_values = vec![vec![0.0; self.densities[0]]; self.densities[1]];
                let mut y_values = vec![vec![0.0; self.densities[0]]; self.densities[1]];

                for anc_elem_id in self.domain.mesh.ancestor_elems(shell_elem.id, true).unwrap().iter() {
                    let (bf, [u_idx_range, v_idx_range]): (BasisFn<SF>, _) = BasisFn::over_anc(
                        i_max as usize,
                        j_max as usize,
                        false,
                        &self.points[0],
                        &self.points[1],
                        shell_elem,
                        &self.domain.mesh.elems[*anc_elem_id],
                    );

                    for bs in self.domain.local_basis_specs(*anc_elem_id).unwrap() {
                        for (m, &m_glob) in u_idx_range.iter().enumerate() {
                            for (n, &n_glob) in v_idx_range.iter().enumerate() {
                                let u = bf.f_u([bs.i as usize, bs.j as usize], [m, n]);
                                let v = bf.f_v([bs.i as usize, bs.j as usize], [m, n]);

                                x_values[m_glob][n_glob] += u.x();
                                x_values[m_glob][n_glob] += v.x();

                                y_values[m_glob][n_glob] += u.y();
                                y_values[m_glob][n_glob] += v.y();
                            }
                        }
                    }
                }

                x_quantity.insert_elem_values(shell_elem.id, x_values);
                y_quantity.insert_elem_values(shell_elem.id, y_values);
            }

            self.quantities.insert(x_q_name.clone(), x_quantity);
            self.quantities.insert(y_q_name.clone(), y_quantity);

            Ok([x_q_name, y_q_name])
        }
    }
}

struct FieldQuantity {
    values: BTreeMap<usize, Vec<Vec<f64>>>,
    name: String,
}

impl FieldQuantity {
    pub fn new(name: &String) -> Self {
        Self {
            values: BTreeMap::new(),
            name: name.clone(),
        }
    }

    pub fn insert_elem_values(&mut self, elem_id: usize, values: Vec<Vec<f64>>) {
        if let Some(prev_entry) = self.values.insert(elem_id, values) {
            panic!("Field Quantity '{}' already had values for Elem {}; cannot assign new values!", self.name, elem_id);
        }
    }
}

fn uniform_range(min: f64, max: f64, n: usize) -> Vec<f64> {
    let step = (max - min) / ((n - 1) as f64);
    (0..n).map(|i| (i as f64) * step + min).collect()
}
