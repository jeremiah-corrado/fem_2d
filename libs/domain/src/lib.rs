
mod vector_space_2d;
pub use vector_space_2d::{V2D, M2D};

#[derive(Clone, Copy, Debug)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

impl Point {
    pub const fn at(x: f64, y: f64) -> Self {
        Self {x, y}
    }

    pub const fn from([x, y]: [f64; 2]) -> Self {
        Self {x, y}
    }
}

impl Default for Point {
    fn default() -> Self {
        Self { x: 0.0, y: 0.0 }
    }
}

// not correct!
impl PartialEq for Point {
    fn eq(&self, other: &Self) -> bool {
        true
    }
}

pub struct Elem {
    id: usize,
    p_min: Point,
    p_max: Point,
}

impl Elem {
    pub fn parametric_projection(&self, real: Point) -> V2D {
        assert!(real.x < self.p_max.x && real.x > self.p_min.x, "Real Point is outside Elem {}'s X Bounds; Cannot project to Parametric Space!", self.id);
        assert!(real.y < self.p_max.y && real.y > self.p_min.y, "Real Point is outside Elem {}'s Y Bounds; Cannot project to Parametric Space!", self.id);

        V2D::from(
            [
                map_range(real.x, self.p_min.x, self.p_max.x, -1.0, 1.0),
                map_range(real.y, self.p_min.y, self.p_max.y, -1.0, 1.0),
            ]
        )
    } 

    pub fn parametric_gradient(&self, _: V2D) -> M2D {
        let dx_du = (self.p_max.x - self.p_min.x) / 2.0;
        let dy_dv = (self.p_max.y - self.p_min.y) / 2.0;


        M2D::from(
            [dx_du, 0.0],
            [0.0, dy_dv]
        )
    } 
}

fn map_range(val: f64, in_min: f64, in_max: f64, out_min: f64, out_max: f64) -> f64 {
    (val - in_min) * (out_max - out_min) / (in_max - in_min) + out_min
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
