use super::{ParaDir, Point, M2D, V2D};
use json::JsonValue;
use num_complex::Complex64;
use std::fmt;

#[derive(Debug)]
/// Basic geometric unit of the FEM Domain.
/// Describes the geometric structure and material properties of a rectangle in real space.
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

    pub fn parametric_projection(&self, real: &Point) -> V2D {
        assert!(
            real.x < self.points[3].x && real.x > self.points[0].x,
            "Real Point is outside Elem {}'s X Bounds; Cannot project onto Parametric Space!",
            self.id
        );
        assert!(
            real.y < self.points[3].y && real.y > self.points[0].y,
            "Real Point is outside Elem {}'s Y Bounds; Cannot project onto Parametric Space!",
            self.id
        );

        V2D::from([
            map_range(real.x, self.points[0].x, self.points[3].x, -1.0, 1.0),
            map_range(real.y, self.points[0].y, self.points[3].y, -1.0, 1.0),
        ])
    }

    pub fn parametric_gradient(&self, _: V2D) -> M2D {
        let dx_du = (self.points[3].x - self.points[0].x) / 2.0;
        let dy_dv = (self.points[3].y - self.points[0].y) / 2.0;

        M2D::from([dx_du, 0.0], [0.0, dy_dv])
    }

    pub fn order_points(&self, p0: &Point, p1: &Point) -> std::cmp::Ordering {
        match p0.orientation_with(p1) {
            ParaDir::U => p0.x_order(p1),
            ParaDir::V => p0.y_order(p1),
        }
    }

    /// Produce a Json Object that describes this Element
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
