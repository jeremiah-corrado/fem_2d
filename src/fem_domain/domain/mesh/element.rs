use super::super::space::{ParaDir, Point, M2D, V2D};
use json::{object, JsonValue};
use num_complex::Complex64;
use std::fmt;

/// The `Element`s are the basic geometric unit of the Mesh in Real Space.
///
/// Elements are responsible for:
/// * Keeping a mapping between Real and Parametric Space in their region of the Mesh (curvilinear Elements are not fully implemented yet)
/// * Keeping track of the material parameters in their portion of the Mesh
///
/// JSON mesh files describe the `Element`s in the domain; not the `Elem`s
/// Upon `Mesh` construction, each `Element` has one associated `Elem`, but more can be added through h-Refinements
#[derive(Debug)]
pub struct Element {
    pub id: usize,
    pub points: [Point; 4],
    pub materials: Materials,
}

impl Element {
    /// Create a new element defined by its coordinates in real space and its material properties
    pub fn new(id: usize, points: [Point; 4], materials: Materials) -> Self {
        Self {
            id,
            points,
            materials,
        }
    }

    // TODO: update this method to support curvilinear Elements
    /// Get the mapping between Real and Parametric Space in the Element
    pub fn parametric_mapping(
        &self,
        _: V2D,
        [[u_min, u_max], [v_min, v_max]]: [[f64; 2]; 2],
    ) -> M2D {
        let real_x_min = map_range(u_min, -1.0, 1.0, self.points[0].x, self.points[3].x);
        let real_x_max = map_range(u_max, -1.0, 1.0, self.points[0].x, self.points[3].x);

        let real_y_min = map_range(v_min, -1.0, 1.0, self.points[0].y, self.points[3].y);
        let real_y_max = map_range(v_max, -1.0, 1.0, self.points[0].y, self.points[3].y);

        // println!("{} \t min (x: {:.5} y: {:.5})  max (x: {:.5} y: {:.5})", self.id, real_x_min, real_y_min, real_x_max, real_y_max);

        let dx_du = (real_x_max - real_x_min) / 2.0;
        let dy_dv = (real_y_max - real_y_min) / 2.0;

        M2D::from([dx_du, 0.0], [0.0, dy_dv])
    }

    // TODO: update this method to support curvilinear Elements
    /// Get the ordering of two points within the Element
    ///
    /// Points closer to the origin (0.0, 0.0) are smaller than points further from the origin
    pub fn order_points(&self, p0: &Point, p1: &Point) -> std::cmp::Ordering {
        match p0.orientation_with(p1) {
            ParaDir::U => p0.x_order(p1),
            ParaDir::V => p0.y_order(p1),
        }
    }

    /// Produce a Json Object that describes this Element
    #[cfg(feature = "json_export")]
    pub fn to_json(&self) -> JsonValue {
        object! {
            "id": self.id,
            "eps_rel": self.materials.eps_rel.re,
            "mu_rel": self.materials.mu_rel.re,
            "eps_rel_im": self.materials.eps_rel.im,
            "mu_rel_im": self.materials.mu_rel.im,
        }
    }
}

fn map_range(val: f64, in_min: f64, in_max: f64, out_min: f64, out_max: f64) -> f64 {
    (val - in_min) * (out_max - out_min) / (in_max - in_min) + out_min
}

/// Complex valued material parameters
#[derive(Clone, Debug)]
pub struct Materials {
    /// Relative Permittivity (ε_r)
    pub eps_rel: Complex64,
    /// Relative permeability (μ_r)
    pub mu_rel: Complex64,
}

impl Materials {
    pub fn from_array(properties: [f64; 4]) -> Self {
        Self {
            eps_rel: Complex64::new(properties[0], properties[1]),
            mu_rel: Complex64::new(properties[2], properties[3]),
        }
    }
}

impl Default for Materials {
    fn default() -> Self {
        Self {
            eps_rel: Complex64::from(1.0),
            mu_rel: Complex64::from(1.0),
        }
    }
}

impl fmt::Display for Materials {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "(ε_re: {}, μ_re: {})", self.eps_rel, self.mu_rel)
    }
}
