use super::{Point, V2D, M2D};
use std::fmt;
use num_complex::Complex64;

#[derive(Debug)]
pub struct Element {
    id: usize,
    points: [Point; 4],
}

impl Element {
    pub fn parametric_projection(&self, real: Point) -> V2D {
        assert!(real.x < self.points[3].x && real.x > self.points[0].x, "Real Point is outside Elem {}'s X Bounds; Cannot project to Parametric Space!", self.id);
        assert!(real.y < self.points[3].y && real.y > self.points[0].y, "Real Point is outside Elem {}'s Y Bounds; Cannot project to Parametric Space!", self.id);

        V2D::from(
            [
                map_range(real.x, self.points[0].x, self.points[3].x, -1.0, 1.0),
                map_range(real.y, self.points[0].y, self.points[3].y, -1.0, 1.0),
            ]
        )
    } 

    pub fn parametric_gradient(&self, _: V2D) -> M2D {
        let dx_du = (self.points[3].x - self.points[0].x) / 2.0;
        let dy_dv = (self.points[3].y - self.points[0].y) / 2.0;


        M2D::from(
            [dx_du, 0.0],
            [0.0, dy_dv]
        )
    } 
}

fn map_range(val: f64, in_min: f64, in_max: f64, out_min: f64, out_max: f64) -> f64 {
    (val - in_min) * (out_max - out_min) / (in_max - in_min) + out_min
}


#[derive(Clone, Debug)]
pub struct Materials {
    pub eps_rel: Complex64,
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