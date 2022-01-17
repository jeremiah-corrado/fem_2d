use std::cmp::Ordering;
use std::fmt;
use std::hash::{Hash, Hasher};
use std::ops::{Add, Div, Index, Mul};

use std::f64::consts::FRAC_PI_4;

use json::JsonValue;

#[derive(Clone, Copy, Debug)]
/// 2D vector in Parametric Space
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
        Self { inner: [0.0; 2] }
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
            inner: [self[0] * other[0], self[1] * other[1]],
        }
    }
}

/*
    | [x1, y1] |
    | [x2, y2] |
*/

#[derive(Clone, Copy, Debug)]
/// 2 by 2 Matrix. Used to represent transformations in/into Parametric space
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

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
/// Parametric Coordinate Directions
pub enum ParaDir {
    U,
    V,
}

impl Into<JsonValue> for ParaDir {
    fn into(self) -> JsonValue {
        match self {
            Self::U => JsonValue::from("U-Dir"),
            Self::V => JsonValue::from("V-Dir"),
        }
    }
}

const POINT_UNIQUENESS_ACCURACY: f64 = 1e-12;

#[derive(Clone, Copy, Debug)]
/// Point in 2D Space
pub struct Point {
    pub x: f64,
    pub y: f64,
    x_cmp: FloatRep,
    y_cmp: FloatRep,
}

impl Point {
    pub fn new(x: f64, y: f64) -> Self {
        Self {
            x,
            y,
            x_cmp: FloatRep::from(x),
            y_cmp: FloatRep::from(y),
        }
    }

    pub fn between(a: &Self, b: &Self) -> Self {
        Self::new((a.x + b.x) / 2.0, (a.y + b.y) / 2.0)
    }

    pub fn from([x, y]: [f64; 2]) -> Self {
        Self {
            x,
            y,
            x_cmp: FloatRep::from(x),
            y_cmp: FloatRep::from(y),
        }
    }

    /// The orientation of the "edge" composed of these two Points
    /// panics if the points have the same location
    pub fn orientation_with(&self, other: &Self) -> ParaDir {
        assert!(
            self.x_cmp != other.x_cmp || self.y_cmp != other.y_cmp,
            "Cannot compute the orientation between two Points at the same location: [{} == {}]!",
            self,
            other,
        );

        let dx = (other.x - self.x).abs();
        let dy = (other.y - self.y).abs();
        let theta = (dy / dx).atan();

        if theta < FRAC_PI_4 {
            ParaDir::U
        } else {
            ParaDir::V
        }
    }

    pub fn dist(&self, other: &Self) -> f64 {
        let dx = (other.x - self.x).abs();
        let dy = (other.y - self.y).abs();

        (dx.powi(2) + dy.powi(2)).sqrt()
    }

    pub fn x_order(&self, other: &Self) -> Ordering {
        self.x_cmp.partial_cmp(&other.x_cmp).unwrap()
    }

    pub fn y_order(&self, other: &Self) -> Ordering {
        self.y_cmp.partial_cmp(&other.y_cmp).unwrap()
    }
}

impl Default for Point {
    fn default() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            x_cmp: FloatRep::from(0.0),
            y_cmp: FloatRep::from(0.0),
        }
    }
}

impl Add for Point {
    type Output = Self;

    fn add(self, other: Self) -> Self {
        Self::new(self.x + other.x, self.y + other.y)
    }
}

impl Div<f64> for Point {
    type Output = Self;

    fn div(self, divis: f64) -> Self {
        Self::new(self.x / divis, self.y / divis)
    }
}

impl Hash for Point {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.x_cmp.hash(state);
        self.y_cmp.hash(state);
    }
}

impl PartialEq for Point {
    fn eq(&self, other: &Self) -> bool {
        self.x_cmp.eq(&other.x_cmp) && self.y_cmp.eq(&other.y_cmp)
    }
}

impl Into<JsonValue> for Point {
    fn into(self) -> JsonValue {
        object! {
            "x": self.x,
            "y": self.y,
        }
    }
}

impl fmt::Display for Point {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "(x: {:.10}, y: {:.10})",
            self.x,
            self.y,
        )
    }
}

#[derive(Hash, PartialEq, Eq, Clone, Copy, Debug)]
struct FloatRep {
    sign: bool,
    bits: u64,
}

impl FloatRep {
    pub fn from(value: f64) -> Self {
        let integer_part = value.abs().trunc();
        let fractional_rounded =
            (value.abs().fract() / POINT_UNIQUENESS_ACCURACY).round() * POINT_UNIQUENESS_ACCURACY;
        let total_rounded = integer_part + fractional_rounded;

        Self {
            sign: value.is_sign_positive(),
            bits: total_rounded.to_bits(),
        }
    }
}

impl Ord for FloatRep {
    fn cmp(&self, other: &Self) -> Ordering {
        match (self.sign, other.sign) {
            (true, true) => self.bits.cmp(&other.bits),
            (false, true) => Ordering::Less,
            (true, false) => Ordering::Greater,
            (false, false) => self.bits.cmp(&other.bits).reverse(),
        }
    }
}

impl PartialOrd for FloatRep {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(match (self.sign, other.sign) {
            (true, true) => self.bits.cmp(&other.bits),
            (false, true) => Ordering::Less,
            (true, false) => Ordering::Greater,
            (false, false) => self.bits.cmp(&other.bits).reverse(),
        })
    }
}
