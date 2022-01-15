#[macro_use]
extern crate json;
extern crate num_complex;
#[macro_use]
extern crate smallvec;

mod domain;

pub use domain::*;

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        let result = 2 + 2;
        assert_eq!(result, 4);
    }
}
