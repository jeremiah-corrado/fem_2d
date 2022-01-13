
extern crate smallvec;
extern crate num_complex;

mod domain;

pub use domain::{Elem, Point, ParaDir, M2D, V2D};

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
