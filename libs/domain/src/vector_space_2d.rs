use std::fmt;
use std::ops::{Add, Div, Index, Mul};

#[derive(Clone, Copy, Debug)]
pub struct V2D {
    inner: [f64; 2],
}

impl V2D {
    pub const fn from([x, y]: [f64; 2]) -> Self {
        Self { inner: [x, y] }
    }

    pub fn dot_with(&self, other: &Self) -> f64 {
        self[0] * other[0] + self[1] * other[1]
    }

    pub fn dot(a: Self, b: Self) -> f64 {
        a[0] * b[0] + a[1] * b[1]
    }

}

impl Default for V2D {
    fn default() -> Self {
        Self {
            inner: [0.0; 2]
        }
    }
}

impl Index<usize> for V2D {
    type Output = f64;
    fn index(&self, index: usize) -> &Self::Output {
        &self.inner[index]
    }
}

impl Add for V2D {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            inner: [self[0] + other[0], self[1] + other[1]],
        }
    }
}

impl Div<f64> for V2D {
    type Output = Self;
    fn div(self, divisor: f64) -> Self {
        Self {
            inner: [self[0] / divisor, self[1] / divisor],
        }
    }
}

impl Div<Self> for V2D {
    type Output = Self;
    fn div(self, divisor: Self) -> Self {
        Self {
            inner: [self[0] / divisor[0], self[1] / divisor[1]],
        }
    }
}

impl Mul<f64> for V2D {
    type Output = Self;
    fn mul(self, coefficient: f64) -> Self {
        Self {
            inner: [self[0] * coefficient, self[1] * coefficient],
        }
    }
}

impl Mul<Self> for V2D {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self {
            inner: [self[0] * other[0], self[1] * other[1]],
        }
    }
}

impl Mul<&Self> for V2D {
    type Output = Self;
    fn mul(self, other: &Self) -> Self {
        Self {
            inner:
            [
                self[0] * other[0],
                self[1] * other[1],
            ]
        }
    }
}

/*
    | [x1, y1] |
    | [x2, y2] |
*/

#[derive(Clone, Copy, Debug)]
pub struct M2D {
    pub u: V2D,
    pub v: V2D,
}

impl M2D {
    pub const fn from(r0: [f64; 2], r1: [f64; 2]) -> Self {
        Self {
            u: V2D::from(r0),
            v: V2D::from(r1),
        }
    }

    #[inline]
    pub fn det(&self) -> f64 {
        self.u[0] * self.v[1] - self.u[1] * self.v[0]
    }

    pub fn inverse(&self) -> Self {
        Self {
            u: V2D::from([self.v[1], -1.0 * self.u[1]]),
            v: V2D::from([-1.0 * self.v[0], self.u[0]]),
        } / self.det()
    }

    pub fn transpose(&self) -> Self {
        Self {
            u: V2D::from([self.u[0], self.v[0]]),
            v: V2D::from([self.u[1], self.v[1]]),
        }
    }
}

impl Div<f64> for M2D {
    type Output = Self;
    fn div(self, divisor: f64) -> Self {
        Self {
            u: self.u / divisor,
            v: self.v / divisor,
        }
    }
}

impl Mul<Self> for M2D {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self {
            u: V2D::from([self.u.dot_with(&other.u), self.u.dot_with(&other.v)]),
            v: V2D::from([self.v.dot_with(&other.u), self.v.dot_with(&other.v)]),
        }
    }
}

impl Mul<V2D> for M2D {
    type Output = V2D;
    fn mul(self, v: V2D) -> V2D {
        V2D::from([self.u.dot_with(&v), self.v.dot_with(&v)])
    }
}

impl fmt::Display for M2D {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "u: [{:.5}, {:.5}]  v: [{:.5}, {:.5}]",
            self.u[0], self.u[1], self.v[0], self.v[1]
        )
    }
}
