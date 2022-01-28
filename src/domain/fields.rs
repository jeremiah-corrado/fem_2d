use super::{Domain, dof::basis_spec::BasisDir, mesh::space::V2D};
use crate::basis::{BasisFn, ShapeFn};

use std::collections::{BTreeMap, HashMap};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::time::SystemTime;

// TODO: rework UniformFieldSpace and print_to_vtk functions to support curvilinear elements
// TODO: implement a constant (over x and y) density FieldSpace structure which supports field image exports

/// A collection of Field Solutions over a [Domain].
///
/// Solutions can be operated on and printed to VTK files for visualization.
pub struct UniformFieldSpace<'d> {
    quantities: HashMap<String, FieldQuantity>,
    parametric_points: [Vec<f64>; 2],
    densities: [usize; 2],
    domain: &'d Domain,
}

impl<'d> UniformFieldSpace<'d> {
    /// Generate a FieldSpace over a [Domain]
    ///
    /// * `domain`: is a reference to a Domain which must outlive this Structure. It is used for all subsequent computations
    /// * `densities`: the number of evaluation points in the u and v directions respectively. These values apply to Elems on the shell of the Domain (their ancestor Elems will have more evaluation points)
    pub fn new(domain: &'d Domain, densities: [usize; 2]) -> Self {
        Self {
            quantities: HashMap::new(),
            parametric_points: [
                uniform_range(-1.0, 1.0, densities[0]),
                uniform_range(-1.0, 1.0, densities[1]),
            ],
            densities,
            domain,
        }
    }

    // TODO: add option to include z-directed fields (after W-Dir & node-type Basis functions are implemented)

    /// Use an eigenvector and associated [ShapeFn] to compute the X and Y fields over the `Domain`
    ///
    /// The X and Y field quantities will be stored as X_{vector_name} and Y_{vector_name} respectively. The Names are returned in an array in that order.
    pub fn xy_fields<SF: ShapeFn>(
        &mut self,
        vector_name: &'static str,
        eigenvector: Vec<f64>,
    ) -> Result<[String; 2], String> {
        if eigenvector.len() != self.domain.dofs.len() {
            Err(format!(
                "NDofs != Eigenvector length ({} != {}); Cannot compute xy fields over Domain",
                self.domain.dofs.len(),
                eigenvector.len()
            ))
        } else {
            let x_q_name = format!("{}_x", vector_name);
            let y_q_name = format!("{}_y", vector_name);

            let mut x_quantity = FieldQuantity::new(&x_q_name);
            let mut y_quantity = FieldQuantity::new(&y_q_name);

            let [i_max, j_max] = self.domain.mesh.max_expansion_orders();

            for shell_elem in self.domain.mesh.elems.iter().filter(|e| !e.has_children()) {
                let mut x_values = vec![vec![0.0; self.densities[0]]; self.densities[1]];
                let mut y_values = vec![vec![0.0; self.densities[0]]; self.densities[1]];

                for anc_elem_id in self
                    .domain
                    .mesh
                    .ancestor_elems(shell_elem.id, true)
                    .unwrap()
                    .iter()
                {
                    let bf: BasisFn<SF> = BasisFn::mapped_over_desc(
                        i_max as usize,
                        j_max as usize,
                        false,
                        &self.parametric_points[0],
                        &self.parametric_points[1],
                        &self.domain.mesh.elems[*anc_elem_id],
                        Some(shell_elem),
                    );

                    for bs in self.domain.local_basis_specs(*anc_elem_id).unwrap() {
                        for m in 0..self.densities[0] {
                            for n in 0..self.densities[1] {
                                let dof_id = bs.dof_id.unwrap();
                                let value = match bs.dir {
                                    BasisDir::U => bf.f_u([bs.i as usize, bs.j as usize], [m, n]),
                                    BasisDir::V => bf.f_v([bs.i as usize, bs.j as usize], [m, n]),
                                    _ => V2D::from([0.0, 0.0]),
                                } * eigenvector[dof_id];

                                x_values[m][n] += value.x();
                                y_values[m][n] += value.y();
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

    /// create a VTK file at the designated `path` (with the file `name.vtk`) including all Field Quantities
    ///
    /// These files can be plotted using [Visit](https://wci.llnl.gov/simulation/computer-codes/visit)
    pub fn print_all_to_vtk(&self, path: impl AsRef<str>) -> std::io::Result<()> {
        let all_q_names = self.quantities.keys().cloned().collect();
        self.print_quantities_to_vkt(path, all_q_names)
    }

    /// create a VTK file at the designated `path` (with the file `name.vtk`) including a list of Field Quantities
    ///
    /// These files can be plotted using [Visit](https://wci.llnl.gov/simulation/computer-codes/visit)
    pub fn print_quantities_to_vkt(
        &self,
        path: impl AsRef<str>,
        quantity_names: Vec<String>,
    ) -> std::io::Result<()> {
        let output_file = File::create(path.as_ref())?;
        let mut writer = BufWriter::new(&output_file);

        let nx = self.densities[0];
        let ny = self.densities[1];

        // header
        writeln!(writer, "# vtk DataFile Version 3.0")?;
        writeln!(
            writer,
            "# File generated by fem_2d on: {:?}\n",
            SystemTime::now().duration_since(SystemTime::UNIX_EPOCH).unwrap()
        )?;
        writeln!(writer, "ASCII")?;
        writeln!(writer, "DATASET UNSTRUCTURED_GRID")?;

        // points
        let num_shell_elems = self
            .domain
            .mesh
            .elems
            .iter()
            .filter(|e| !e.has_children())
            .count();
        let num_points = nx * ny * num_shell_elems;
        writeln!(writer, "\nPOINTS {} double", num_points)?;
        for shell_elem in self
            .domain
            .mesh
            .elems
            .iter()
            .filter(|elem| !elem.has_children())
        {
            let diag_points = self.domain.mesh.elem_diag_points(shell_elem.id);
            for x in uniform_range(diag_points[0].x, diag_points[1].x, self.densities[0]) {
                for y in uniform_range(diag_points[0].y, diag_points[1].y, self.densities[1]) {
                    writeln!(writer, "{:.10} {:.10} 0.0", x, y)?;
                }
            }
        }

        // cells
        let num_cells = (nx - 1) * (ny - 1) * num_shell_elems;
        writeln!(writer, "\nCELLS {} {}", num_cells, 5 * num_cells)?;
        for k in 0..num_shell_elems {
            for i in 0..(nx - 1) {
                for j in 0..(ny - 1) {
                    let initial_pt = nx * i + j + (nx * ny) * k;

                    writeln!(
                        writer,
                        "4\t{}\t{}\t{}\t{}",
                        initial_pt,
                        initial_pt + 1,
                        initial_pt + nx + 1,
                        initial_pt + nx,
                    )?;
                }
            }
        }

        // cell types
        writeln!(writer, "\nCELL_TYPES {}", num_cells)?;
        for _ in 0..num_cells {
            write!(writer, " 9")?;
        }
        writeln!(writer)?;

        // field values
        writeln!(writer, "POINT_DATA {}", num_points)?;
        for q_name in quantity_names {
            match self.quantities.get(&q_name) {
                Some(field_quant) => {
                    field_quant.write_vtk_quantity(&mut writer)?;
                }
                None => println!(
                    "Field Space does not have Quantity '{}'; cannot write to VTK!",
                    q_name
                ),
            };
        }
        writeln!(writer)?;

        Ok(())
    }

    /// map an operation over a field quantity (`name`) and store the result in a new quantity (`result_name`)
    pub fn map_to_quantity<F>(
        &mut self,
        name: impl AsRef<str>,
        result_name: impl AsRef<str>,
        operator: F,
    ) -> Result<(), String>
    where
        F: Fn(&f64) -> f64 + Copy,
    {
        let q_key = String::from(name.as_ref());
        let q_new_key = String::from(result_name.as_ref());

        if !self.quantities.contains_key(&q_key) {
            Err(format!(
                "FieldSpace does not have quantity: {}; cannot apply operation!",
                q_key
            ))
        } else {
            let q_new = self
                .quantities
                .get(&q_key)
                .unwrap()
                .operation(operator, &q_new_key);
            self.quantities.insert(q_new_key, q_new);
            Ok(())
        }
    }

    /// evaluate an expression of two field quantities and store the result in a new quantity (`result_name`)
    pub fn expression_2arg<F>(
        &mut self,
        operand_names: [impl AsRef<str>; 2],
        result_name: impl AsRef<str>,
        expression: F,
    ) -> Result<(), String>
    where
        F: Fn(f64, f64) -> f64,
    {
        let op_a = String::from(operand_names[0].as_ref());
        let op_b = String::from(operand_names[1].as_ref());
        let q_new_key = String::from(result_name.as_ref());

        if !self.quantities.contains_key(&op_a) || !self.quantities.contains_key(&op_b) {
            Err(format!(
                "FieldSpace does not have quantities {} and {}; cannot apply operation!",
                op_a, op_b
            ))
        } else {
            let mut q_new = FieldQuantity::new(&q_new_key);
            let q_a = self.quantities.get(&op_a).unwrap();
            let q_b = self.quantities.get(&op_b).unwrap();

            for ((shell_elem_id, elem_values_a), elem_values_b) in
                q_a.values.iter().zip(q_b.values.values())
            {
                let mut result_values = vec![vec![0.0; self.densities[0]]; self.densities[1]];
                for m in 0..self.densities[0] {
                    for n in 0..self.densities[1] {
                        result_values[m][n] = expression(elem_values_a[m][n], elem_values_b[m][n]);
                    }
                }

                q_new.insert_elem_values(*shell_elem_id, result_values);
            }

            self.quantities.insert(q_new_key, q_new);
            Ok(())
        }
    }

    // TODO: implement 3arg, Narg, and convolution.
}

struct FieldQuantity {
    pub values: BTreeMap<usize, Vec<Vec<f64>>>,
    name: String,
}

impl FieldQuantity {
    pub fn new(name: &str) -> Self {
        Self {
            values: BTreeMap::new(),
            name: name.to_string(),
        }
    }

    pub fn insert_elem_values(&mut self, elem_id: usize, values: Vec<Vec<f64>>) {
        if self.values.insert(elem_id, values).is_some() {
            panic!(
                "Field Quantity '{}' already had values for Elem {}; cannot assign new values!",
                self.name, elem_id
            );
        }
    }

    pub fn write_vtk_quantity(&self, writer: &mut BufWriter<&File>) -> std::io::Result<()> {
        writeln!(
            writer,
            "SCALARS {} double 1 \nLOOKUP_TABLE default",
            self.name
        )?;

        for (_, shell_elem_values) in self.values.iter() {
            for shell_row_values in shell_elem_values {
                for value in shell_row_values {
                    write!(writer, "{:.15} ", value)?;
                }
            }
        }

        Ok(())
    }

    pub fn operation<F>(&self, operator: F, new_name: &str) -> Self
    where
        F: Fn(&f64) -> f64 + Copy,
    {
        Self {
            values: self
                .values
                .iter()
                .map(|(elem_id, elem_values)| {
                    (
                        *elem_id,
                        elem_values
                            .iter()
                            .map(|col| col.iter().map(operator).collect())
                            .collect(),
                    )
                })
                .collect(),
            name: new_name.to_string(),
        }
    }
}

fn uniform_range(min: f64, max: f64, n: usize) -> Vec<f64> {
    let step = (max - min) / ((n - 1) as f64);
    (0..n).map(|i| (i as f64) * step + min).collect()
}
