
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PolyOrders {
    i: u8,
    j: u8,
}

impl PolyOrders {
    fn from(i: u8, j: u8) -> Self {
        Self { i, j }
    }

    fn refined(&self, refinement: PRef) -> Result<Self, ()> {
        Ok(Self {
            i: refinement.refine_i(self.i)?,
            j: refinement.refine_j(self.j)?,
        })
    }
}

impl Default for PolyOrders {
    fn default() -> Self {
        Self {
            i: 1,
            j: 1,
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum PRefInt {
    Increment(u8),
    Decrement(u8),
    None,
}

impl PRefInt {
    fn refine(&self, n: u8) -> Result<u8, ()> {
        match self {
            Self::Increment(delta) => Ok(n + *delta),
            Self::Decrement(delta) => {
                if *delta >= n {
                    Err(())
                } else {
                    Ok(n - *delta)
                }
            }
            Self::None => Ok(n),
        }
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PRef {
    di: PRefInt,
    dj: PRefInt,
}

impl PRef {
    pub const fn from(i: i8, j: i8) -> Self {
        Self {
            di: match i {
                0 => PRefInt::None,
                d if d > 0 => PRefInt::Increment(d as u8),
                d if d < 0 => PRefInt::Decrement((-1 * d) as u8),
                _ => unreachable!()
            },
            dj: match j {
                0 => PRefInt::None,
                d if d > 0 => PRefInt::Increment(d as u8),
                d if d < 0 => PRefInt::Decrement((-1 * d) as u8),
                _ => unreachable!()
            }
        }
    }

    fn refine_i(&self, i_current: u8) -> Result<u8, ()> {
        self.di.refine(i_current)
    }

    fn refine_j(&self, j_current: u8) -> Result<u8, ()> {
        self.dj.refine(j_current)
    }

}